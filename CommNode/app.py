from flask import Flask, request
import requests as http_req
import hashlib
import copy
import json
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

db_nodes = [
    
"http://18.220.28.88/"

]

node_count = len(db_nodes)

max_act_prot = 5

ready_results = {}


# plugin to comm node: POST http://<com node ip>:5000/plugin_request
# w/ JSON Data: {'sequence':<FASTA FILE>}

# db node to comm node: POST http://<com node ip>:5000/node_data?qid=<seq_hash>
# w/ JSON Data (10 obj array):
#               [{'Accession #': n, 'PCT':0.7, 'Q_PCT':0.85, 'Score':238},
#                {'Accession #': n, 'PCT':0.7, 'Q_PCT':0.85, 'Score':238},... x10 ]


def get_top_ten_results(results_list, qid):
    # print(results_list)
    combined_results_list = []
    for result in results_list:
        combined_results_list.extend(result)
    sorted_results_list = sorted(combined_results_list, key=lambda res: res['score'], reverse=True)
    if len(sorted_results_list) > 10:
        sorted_results_list = sorted_results_list[0:10]

    return {"qid":qid,"data":sorted_results_list}


class QueryTracker:

    def __init__(self, max_act_proc):
        self.query_process_list = [-1]*max_act_proc
        ## Make queue also store the sequence data ... otherwise cannot process it ...
        self.query_process_queue = []
        self.query_process_results = [[] for i in range(max_act_proc)]


    def exists(self, qid, check_list=True, check_queue=True):
        # Checks if query exists in queue or list
        if (check_list and (qid in self.query_process_list)) or (check_queue and any(qid == entry['qid'] for entry in self.query_process_queue )):
            return True
        else:
            return False
    def queue_len(self):
        return len(self.query_process_queue)

    def new(self, qid, sequence):
        # Returns -1 if duplicate
        if self.exists(qid):
            return -1
        # Returns 1 if stored in active process list
        for i,ele in enumerate(self.query_process_list):
            if self.query_process_list[i] == -1:
                self.query_process_list[i] = qid
                return 1
        # Returns 0 if stored in process queue
        self.query_process_queue.append({"qid": qid, "sequence": sequence})
        return 0

    def store_results(self,qid,result):
        # Stores the results in the corresponding location given the qid
        ind = self.query_process_list.index(qid)
        self.query_process_results[ind].append(result)
        

    def all_results_received(self,qid):
        # Checks if all results were received for a given qid
        ind = self.query_process_list.index(qid)
        if len(self.query_process_results[ind]) == node_count:
            return True
        return False

    def pop_results(self, qid):
        # Returns the results for a given qid and frees its space
        ind = self.query_process_list.index(qid)
        results = copy.deepcopy(self.query_process_results[ind])
        self.query_process_list[ind] = -1
        self.query_process_results[ind] = []
        return results
    
    def status(self, qid):
        if not self.exists(qid):
            return -2
        print(self.query_process_list)
        if qid in self.query_process_list:
            ind = self.query_process_list.index(qid)
            print(ind)
            return len(self.query_process_results[ind])
        for item in self.query_process_queue:
            if qid in item['qid']:
                return -1
        return -2

    def insert_proc_from_queue(self):
        # Moves a process from the queue into the list of active processes
        if self.queue_len() == 0:
            return False
        for i,ele in enumerate(self.query_process_list):
            if self.query_process_list[i] == -1:
                popped_data = self.query_process_queue.pop(0)
                self.query_process_list[i] = popped_data['qid']
                return popped_data
        return False

qtrack = QueryTracker(max_act_prot)

@app.route('/status')
def index():
    print("Got status")
    return 'The server is running', 200

@app.route('/plugin_request', methods=['GET', 'POST'])
def plugin_request():
    # Process a plugin request
    if request.method == 'POST':
        # Get the data being sent from the plugin
        sequence = request.data.decode('UTF-8')
        # Store sequence as a md5 hash
        seq_hash = hashlib.md5(sequence.encode()).hexdigest()
        print("Got query from plugin "+seq_hash)
        # Add sequence hash to query tracker
        status = qtrack.new(seq_hash, sequence)
        # Check if duplicate request
        if(status == -1):
            return "Error: Duplicate request", 400
        # Check if request enqueued
        if(status == 0):
            return "Job Enqueued: "+seq_hash,250
        # Check if request stored in active process list
        if(status == 1):
            for db_node in db_nodes:
                print("Sending to "+db_node)
                url_post = db_node+"api/request/"+seq_hash
                header = {'Content-Type':'text/plain'}
                response = http_req.post(url_post, sequence, headers=header)
            return seq_hash, 200


@app.route('/node_data/<qid>', methods=['GET','POST'])
def node_data(qid):
    if request.method == 'POST':
        # Process node data
        seq_hash = qid
        received_data = request.get_json()
        node_id = received_data['nid']
        print("Received results from "+ str(node_id))
        results = received_data['results']
        qtrack.store_results(seq_hash, results)
        if qtrack.all_results_received(qid):
            results_list = qtrack.pop_results(qid)
            # Process results, then store with qid, ready for javascript
            data_ready = get_top_ten_results(results_list, qid)
            ready_results[qid] = data_ready
            print("Results ready to be read")
            # Pop process from queue if there are waiting queries
            new_id = qtrack.insert_proc_from_queue()
            if new_id:
                # Send request to process new query
                for db_node in db_nodes:
                    url_post = db_node+"api/request/"+new_id['qid']
                    header = {'Content-Type':'text/plain'}
                    response = http_req.post(url_post, new_id['sequence'], headers=header)

            return "sent",200
        
        else:
            return "waiting", 250 


@app.route('/plugin_poll/<qid>')
def plugin_poll(qid):
    # Check if processed results ready, if so, return genbank data with query number, otherwise, return just query number
    # and number of nodes waiting to hear back
    if qid in ready_results:
        payload = ready_results.pop(qid)
        # Check if needs to be turned into string
        return json.dumps(payload), 200
    else:
        status = qtrack.status(qid)
        if status == -2:
            return "Error: Query not found in queue", 400
        if status == -1:
            return "Query still in queue", 250
        else:
            return f"{status} out of {node_count} processes finished", 250


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=80)
