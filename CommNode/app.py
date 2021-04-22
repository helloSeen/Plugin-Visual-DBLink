from flask import Flask, request
import requests as http_req
import hashlib
import copy
import json
import re
from flask_cors import CORS

app = Flask(__name__)
CORS(app)


api_key = "f138e8c165aa0e7a283a3d7a72aca89c3908"

__node_list = [
    
"http://3.143.3.150/"

]

db_nodes = copy.deepcopy(__node_list)

node_count = len(db_nodes)

max_act_prot = 5

ready_results = {}


# plugin to comm node: POST http://<com node ip>:5000/plugin_request
# w/ JSON Data: {'sequence':<FASTA FILE>}

# db node to comm node: POST http://<com node ip>:5000/node_data?qid=<seq_hash>
# w/ JSON Data (10 obj array):
#               [{'Accession #': n, 'PCT':0.7, 'Q_PCT':0.85, 'Score':238},
#                {'Accession #': n, 'PCT':0.7, 'Q_PCT':0.85, 'Score':238},... x10 ]


def get_info_from_accession_ids(results_list, api_key_string=None):
    # Stores results list for modification later
    payload = copy.deepcopy(results_list)
    
    # Gets list of accession IDs
    accession_id_list = []
    for doc in payload['results']:
        accession_id_list.append(doc['accession'])
    
    # Joins accession ID list into one comma-separated string
    sep = ","
    accession_ids=sep.join(accession_id_list)
    
    # Creates URL for GET request
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    # Nucleotide DB
    database_name = "nuccore"
    # Request GenBank File
    rettype = "gb"
    retmode = "text"

    fetch_url = f"efetch.fcgi?db={database_name}&id={accession_ids}&rettype={rettype}&retmode={retmode}"

    api_key=""
    if api_key_string is not None:
        api_key=f"&api_key={api_key_string}"

    resp = http_req.get(base_url+fetch_url+api_key)
    
    combined_gb_files = resp.content.decode()
    # Splits into list of strings with file data, removes trailing newline from list
    gb_file_list = combined_gb_files.split("//\n")
    if '\n' in gb_file_list:
        gb_file_list.remove('\n')
                          
    for i,gb_file in enumerate(gb_file_list):
        
        info_dict = {}
        
        locus_str=""
        def_str=""
        source_str=""
        ref_str=""
        auth_str=""
        cons_str=""
        title_str=""
        journ_str=""
        pubmed_str=""

        # Parses document for information    
        pattern_locus = re.compile(r'(?<=^LOCUS\s{7})([\S]+?)(?=\s+)', re.MULTILINE)
        pattern_def = re.compile(r'(?<=^DEFINITION\s{2})([\s\S]+?)(?=\n[A-Z]{2,})', re.MULTILINE)
        pattern_source = re.compile(r'(?<=^SOURCE\s{6})([\s\S]+?)(?=\n\s*[A-Z]{2,})', re.MULTILINE)
        pattern_ref = re.compile(r'(?<=REFERENCE\s{3}1\s{2})(.*)(?=\n)$', re.MULTILINE)
        pattern_auth = re.compile(r'(?<=\s{2}AUTHORS\s{3})([\s\S]+?)(?=\s{3}[A-Z]{2,})', re.MULTILINE)
        pattern_cons = re.compile(r'(?<=\s{2}CONSRTM\s{3})([\s\S]+?)(?=\s{3}[A-Z]{2,})', re.MULTILINE)
        pattern_title = re.compile(r'(?<=\s{2}TITLE\s{5})([\s\S]+?)(?=\s{3}[A-Z]{2,})', re.MULTILINE)
        pattern_journ = re.compile(r'(?<=\s{2}JOURNAL\s{3})([\s\S]+?)(?=\n\s*[A-Z]{2,})', re.MULTILINE)
        pattern_pubmed = re.compile(r'(?<=PUBMED\s{3})([\w\d]+)(?=\n)$', re.MULTILINE)
        
        locus_regex = re.search(pattern_locus, gb_file)
        def_regex = re.search(pattern_def, gb_file)
        source_regex = re.search(pattern_source, gb_file)
        ref_regex = re.search(pattern_ref, gb_file)
        auth_regex = re.search(pattern_auth, gb_file)
        cons_regex = re.search(pattern_cons, gb_file)
        title_regex = re.search(pattern_title, gb_file)
        journ_regex = re.search(pattern_journ, gb_file)
        pubmed_regex = re.search(pattern_pubmed, gb_file)
        
        # Constructs dictionary with found parameters, removes trailing spaces
        info_dict['Link']="https://www.ncbi.nlm.nih.gov/nuccore/"+accession_id_list[i]

        if locus_regex:
            locus_str = locus_regex.group(0).replace("\n           ", "")
            info_dict['Locus'] = locus_str
        
        if def_regex:
            def_str = def_regex.group(0).replace("\n           ", "")
            info_dict['Definition'] = def_str
            
        if source_regex:
            source_str = source_regex.group(0).replace("\n           ", "")
            info_dict['Source'] = source_str
        
        if ref_regex:
            ref_str = ref_regex.group(0).replace("\n           ", "")
            info_dict['Reference'] = ref_str
            
        if auth_regex:
            auth_str = auth_regex.group(0).replace("\n           ", "")
            info_dict['Authors'] = auth_str
            
        if cons_regex:
            cons_str = cons_regex.group(0).replace("\n           ", "")
            info_dict['Consortium'] = cons_str    
            
        if title_regex:
            title_str = title_regex.group(0).replace("\n           ", "")
            info_dict['Publication Title'] = title_str
        
        if journ_regex:
            journ_str = journ_regex.group(0).replace("\n           ", "")
            info_dict['Journal'] = journ_str
            
        if pubmed_regex:
            pubmed_str = pubmed_regex.group(0).replace("\n           ", "")
            info_dict['PubMed ID'] = pubmed_str
        
        # Adds "data" field to payload, then stores the document data there
        if accession_id_list[i] == payload["results"][i]["accession"]:
            payload["results"][i]["data"] = info_dict
        else:
            print("Error!")
        
    return payload


def get_top_ten_results(results_list, qid):
    # print(results_list)
    combined_results_list = []
    for result in results_list:
        combined_results_list.extend(result)
    sorted_results_list = sorted(combined_results_list, key=lambda res: res['score'], reverse=True)
    if len(sorted_results_list) > 10:
        sorted_results_list = sorted_results_list[0:10]

    return {"qid":qid,"results":sorted_results_list}


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
    return {"status":'The server is running'}, 200

@app.route('/plugin_request', methods=['GET', 'POST'])
def plugin_request():
    # Process a plugin request
    if request.method == 'POST':
        for db_node in db_nodes:
            try:
                response = http_req.get(db_node+"status", timeout=10)
            except (http_req.ReadTimeout, http_req.ConnectionError):
                db_nodes.remove(db_node)
                continue
    if not db_nodes:
        return {"status":"No database nodes active"},500

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
        return {"qid":seq_hash}, 200


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
            sorted_results = get_top_ten_results(results_list, qid)
            data_ready = get_info_from_accession_ids(sorted_results, api_key)
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

            return {"status":"sent"},200
        
        else:
            return {"status":"waiting"}, 250 


@app.route('/plugin_poll/<qid>')
def plugin_poll(qid):
    # Check if processed results ready, if so, return genbank data with query number, otherwise, return just query number
    # and number of nodes waiting to hear back
    if qid in ready_results:
        payload = ready_results.pop(qid)
        payload["State"] = "Done"
        # Check if needs to be turned into string
        return json.dumps(payload), 200
    else:
        status = qtrack.status(qid)
        if status == -2:
            return {"State":"Error: Query not found in queue"}, 400
        if status == -1:
            return {"State": "Query still in queue"}, 250
        else:
            return {"State": f"{status} out of {node_count} processes finished"}, 250


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=80)
