from flask import Flask, request
import requests as http_req
import hashlib
import copy
import json
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

ready_results = []
status = -1
count  = 0

@app.route('/plugin_poll/<qid>')
def plugin_poll(qid):
    # Check if processed results ready, if so, return genbank data with query number, otherwise, return just query number
    # and number of nodes waiting to hear back
    global count
    global status
    global ready_results
    count += 1
    if qid in ready_results:
        ready_results = []
        with open("./tst2.json") as f:
            payload = json.load(f)
        payload["State"] = "Done"
        # Check if needs to be turned into string
        return json.dumps(payload), 200
    else:
        if status == -2:
            return {"State":"Error: Query not found in queue"}, 400
        if status == -1:
            if count == 5:
                status = 0
            return {"State": "Query still in queue"}, 250
        else:
            if count == 7:
                ready_results.append(qid)
                status = -1
                count  = 0
            return {"State": f"{count} out of 15 processes finished"}, 250


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)