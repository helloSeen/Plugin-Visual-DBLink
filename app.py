from flask import Flask, request, abort, jsonify
import os
import requests as http_req
import hashlib
import copy
import json
from flask_cors import CORS

app = Flask(__name__)

CORS(app)

commNode_url="http://18.190.154.161/"

@app.route("/status")
def status():
    return("The Visualisation Test Plugin Flask Server is up and running")

@app.route("/evaluate", methods=["POST", "GET"])
def evaluate():
    data = request.get_json(force=True)
    rdf_type = data['type']
    
    ########## REPLACE THIS SECTION WITH OWN RUN CODE #################
    #uses rdf types
    accepted_types = {'Component'}
    
    acceptable = rdf_type in accepted_types
    
    # #to ensure it shows up on all pages
    # acceptable = True
    ################## END SECTION ####################################
    
    if acceptable:
        return f'The type sent ({rdf_type}) is an accepted type', 200
    else:
        return f'The type sent ({rdf_type}) is NOT an accepted type', 415

@app.route("/run", methods=["POST"])
def run():
    data = request.get_json(force=True)
    
    top_level_url = data['top_level']
    complete_sbol = data['complete_sbol']
    instance_url = data['instanceUrl']
    size = data['size']
    rdf_type = data['type']
    shallow_sbol = data['shallow_sbol']
    
    url = complete_sbol.replace('/sbol','')
    
    resp=http_req.get(url+r'/fasta', timeout=10)

    fasta_file = resp.content
    header = {'Content-Type':'text/plain'}
    response = http_req.post(commNode_url+"plugin_request", fasta_file, headers=header, timeout=10)
    qid = response.json['qid']
    cwd = os.getcwd()
    filename = os.path.join(cwd, "index2.html")
    
    try:
        with open(filename, 'r') as htmlfile:
            result = htmlfile.read()
            
        #put in the url, uri, and instance given by synbiohub
        result = result.replace("COMM_NODE_IP", commNode_url)
        result = result.replace("QUERY_ID", qid)
        # result = result.replace("INSTANCE_REPLACE", instance_url)
           
        # result = result.replace("REQUEST_REPLACE", str(data))
            
        return result
    except Exception as e:
        print(e)
        abort(400)
