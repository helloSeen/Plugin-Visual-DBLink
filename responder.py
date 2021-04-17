from flask import Flask, request, make_response
import json
import docker
import concurrent.futures
import os
import time

app = Flask(__name__)
db = 'nt.00'
nid = '1'
HOME = '$HOME/blast/'
#HOME = "C:\\Users\\Sean Nemtzow\\Documents\\blast\\"
# in seconds
tout = 60


def run_docker(qid):
    """ This function runs a blast search in a docker container, returns the top 10 results
        with score, query coverage, and percent identity
    """
    fasta = "/blast/queries/{}.fsa".format(qid)
    results = "/blast/results/{}.out".format(qid)
    docker_client = docker.from_env()
    cmnd = "blastn -query {} -db {} -out {} -outfmt \"6 sacc score qcovhsp pident\"".format(
        fasta, db, results)
    volume_dict = {os.path.join(HOME, 'blastdb'): {'bind': '/blast/blastdb', 'mode': 'ro'},
                   os.path.join(HOME, 'blastdb_custom'): {'bind': '/blast/blastdb_custom', 'mode': 'ro'},
                   os.path.join(HOME, 'queries'): {'bind': '/blast/queries', 'mode': 'ro'},
                   os.path.join(HOME, 'results'): {'bind': '/blast/results', 'mode': 'rw'}
                   }

    container = docker_client.containers.run(
        image='ncbi/blast', command=cmnd, volumes=volume_dict, detach=True)
    
    timeout_reached = True
    for i in range(tout):
        if container.status == 'exited':
            timeout_reached = False
            break
        time.sleep(1)
    container.remove(force=True)

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
                metrics['accession'] = values[0]
                metrics['score'] = int(values[1])
                metrics['per_cov'] = float(values[2])
                metrics['per_id'] = float(values[3])
                r_dict['results'].append(metrics)

    # Delete files when complete
    os.remove(local_r)
    os.remove(local_f)
    return r_dict


@app.route('/status')
def healthy():
    return "nice and healthy", 200


@app.route('/api/request/<qid>', methods=['POST', 'GET'])
def process_request(qid):
    content = request.data.decode('UTF-8')
    fasta = os.path.join(HOME, "queries", "{}.fsa".format(qid))
    with open(fasta, "w+") as fast_f:
        fast_f.write(content)
    with concurrent.futures.ThreadPoolExecutor() as executor:
        future = executor.submit(run_docker, qid)
        metrics = future.result()
        return metrics


def main():
    app.run(host='0.0.0.0',port=80, debug=True)


if __name__ == "__main__":
    main()
