from flask import Flask, jsonify, request
from flask_cors import CORS
import requests

app = Flask(__name__)
CORS(app)

# --- YOUR HOMEWORK LOGIC HERE ---
class ProteinAPI:
    def __init__(self):
        self.base_url = "https://files.rcsb.org/download"

    def get_structure(self, pdb_id):
        url = f"{self.base_url}/{pdb_id}.pdb"
        try:
            r = requests.get(url, timeout=10)
            if r.status_code == 200:
                return self.parse_backbone(r.text)
            return []
        except:
            return []

    def parse_backbone(self, text):
        backbone = []
        last_chain = None
        for line in text.splitlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    chain_id = line[21]
                    # If chain changes, add a "break" marker
                    if last_chain is not None and chain_id != last_chain:
                        backbone.append(None)
                    
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    backbone.append((x, y, z))
                    last_chain = chain_id
                except:
                    continue
        return backbone

api = ProteinAPI()

# --- THE API ENDPOINT ---
@app.route('/api/molecule/<pdb_id>')
def get_molecule(pdb_id):
    # Get the 3D points using your logic
    points = api.get_structure(pdb_id)
    
    # Calculate the center (for the frontend to center the camera)
    valid_points = [p for p in points if p is not None]
    if not valid_points:
        return jsonify({"error": "No data found"}), 404
        
    avg_x = sum(p[0] for p in valid_points) / len(valid_points)
    avg_y = sum(p[1] for p in valid_points) / len(valid_points)
    avg_z = sum(p[2] for p in valid_points) / len(valid_points)

    return jsonify({
        "points": points,
        "center": [avg_x, avg_y, avg_z]
    })

if __name__ == '__main__':
    app.run(debug=True, port=5000)
