from flask import Flask, request, abort, jsonify
import os
import requests as http_req
import hashlib
import copy
import json
import traceback
from flask_cors import CORS

app = Flask(__name__)

CORS(app)

# Reads the CommIP.txt for the Communication Server IP address
with open("CommIP.txt", "r") as f:
    comm_ip = f.read()
comm_ip=comm_ip.strip()
commNode_url = f"http://{comm_ip}/"

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
    # Gets current page data from SynBioHub
    data = request.get_json(force=True)
    
    top_level_url = data['top_level']
    complete_sbol = data['complete_sbol']
    instance_url = data['instanceUrl']
    size = data['size']
    rdf_type = data['type']
    shallow_sbol = data['shallow_sbol']
    
    url = complete_sbol.replace('/sbol','')
    
    # Try Except prevents errors from affecting main page
    try:
        cwd = os.getcwd()
        # Gets the fasta file for the current page
        resp=http_req.get(url+r'/fasta', timeout=10)
        fasta_file = resp.content
        # Sends fasta file data to Communication Server, and stores response
        header = {'Content-Type':'text/plain'}
        response = http_req.post(commNode_url+"plugin_request", fasta_file, headers=header, timeout=10)
        resp_content = response.json()

        # If the file was sent successfully, render the plugin's html
        # page and send the data to the page
        if resp_content["status"] == "success":
            qid = response.json()['qid']
            filename = os.path.join(cwd, "html/index.html")
            with open(filename, 'r') as htmlfile:
                result = htmlfile.read()
            result = result.replace("COMM_NODE_IP", commNode_url)
            result = result.replace("QUERY_ID", qid)
        # Otherwise, load the error page
        else:
            with open(os.path.join(cwd, "html/error.html"), 'r') as htmlfile:
                result = htmlfile.read()
                print(resp_content["status"])
        return result
    # If an exception occurs, load the error page and print to
    # the server console the exception's traceback
    except Exception as e:
        with open(os.path.join(cwd, "html/error.html"), 'r') as htmlfile:
            result = htmlfile.read()
        print(traceback.format_exc())
        return result, 299



if __name__ == "__main__":
    app.run(host='0.0.0.0', port=5050)