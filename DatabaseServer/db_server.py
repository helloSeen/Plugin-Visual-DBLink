from flask import Flask, request, make_response
import json
import docker
import threading
import os
import time
from pathlib import Path
import requests as http_req

app = Flask(__name__)
# Configure databse id and genbank fragment
db = 'nt.00'
nid = '1'
HOME = '/home/ec2-user/'
# in seconds
tout = 600


def run_docker(qid, remote_ip):
    """ This function runs a blast search in a docker container, returns the top 10 results
        with score, query coverage, and percent identity
    """
    fasta = "/blast/queries/{}.fsa".format(qid)
    results = "/blast/results/{}.out".format(qid)
    docker_client = docker.from_env()
    # Command to run, it limits results to:
    # Accession ID, Score, Query Coverage, and Identity Percentage
    cmnd = "blastn -query {} -db {} -out {} -outfmt \"6 sacc score qcovhsp pident\"".format(
        fasta, db, results)
    # Mounts local directories inside docker container
    volume_dict = {os.path.join(HOME, 'blastdb'): {'bind': '/blast/blastdb', 'mode': 'ro'},
                   os.path.join(HOME, 'blastdb_custom'): {'bind': '/blast/blastdb_custom', 'mode': 'ro'},
                   os.path.join(HOME, 'queries'): {'bind': '/blast/queries', 'mode': 'ro'},
                   os.path.join(HOME, 'results'): {'bind': '/blast/results', 'mode': 'rw'}
                   }

    container = docker_client.containers.run(
        image='ncbi/blast', command=cmnd, volumes=volume_dict, detach=True)
    
    timeout_reached = True
    # Check status of docker execution every 1 second
    for i in range(tout):
        container.reload()
        print(container.status)
        if container.status == 'exited':
            timeout_reached = False
            break
        time.sleep(1)
    container.remove(force=True)

    # Build json response to server
    local_r = os.path.join(HOME, "results", "{}.out".format(qid))
    local_f = os.path.join(HOME, "queries", "{}.fsa".format(qid))
    r_dict = {}
    r_dict['qid'] = qid
    r_dict['nid'] = nid
    r_dict['results'] = []
    if not timeout_reached:
        # If results exist
        with open(local_r, 'r') as rfile:
            for i in range(10):
                metrics = {}
                line = rfile.readline()[:-1]
                values = line.split("\t")
                if values[0]:
                    metrics['accession'] = values[0]
                    metrics['score'] = int(values[1])
                    metrics['per_cov'] = float(values[2])
                    metrics['per_id'] = float(values[3])
                    r_dict['results'].append(metrics)

    # Delete files when complete
    os.remove(local_r)
    os.remove(local_f)

    #Send response
    url = "http://{}:80/node_data/{}".format(remote_ip, qid)
    http_req.post(url, json=r_dict)


@app.route('/status')
def healthy():
    return "nice and healthy", 200


@app.route('/api/request/<qid>', methods=['POST', 'GET'])
def process_request(qid):
    ''' 
    Receives query requests from server, spins up thread and runs
    docker command in thread
    '''
    print("Got request for: {}".format(qid))
    content = request.data.decode('UTF-8')
    fasta = os.path.join(HOME, "queries", "{}.fsa".format(qid))
    with open(fasta, "w+") as fast_f:
        fast_f.write(content)
    remote_ip = request.remote_addr
    thread = threading.Thread(target=run_docker, args=(qid,remote_ip))
    thread.start()
    return "ok", 200


def main():
    app.run(host='0.0.0.0',port=80, debug=True)


if __name__ == "__main__":
    main()
