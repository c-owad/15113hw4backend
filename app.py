from flask import Flask, jsonify
from flask_cors import CORS
import requests

app = Flask(__name__)
CORS(app)

# --- CONFIGURATION ---
# We use this to center the molecule
def get_center(points):
    valid = [p for p in points if p is not None]
    if not valid: return (0,0,0)
    avg_x = sum(p['x'] for p in valid) / len(valid)
    avg_y = sum(p['y'] for p in valid) / len(valid)
    avg_z = sum(p['z'] for p in valid) / len(valid)
    return (avg_x, avg_y, avg_z)

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
            if line.startswith("ATOM") and line[12:16].strip() == "CA": # Alpha Carbon only
                try:
                    chain_id = line[21]
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    
                    # If chain changes, we insert a None to "break" the ribbon drawing
                    if last_chain is not None and chain_id != last_chain:
                        backbone.append(None)

                    # We now return a DICTIONARY with chain info, not just a tuple
                    backbone.append({"x": x, "y": y, "z": z, "chain": chain_id})
                    last_chain = chain_id
                except:
                    continue
        return backbone

api = ProteinAPI()

@app.route('/api/molecule/<pdb_id>')
def get_molecule(pdb_id):
    points = api.get_structure(pdb_id)
    center = get_center(points)
    
    # Return the data
    return jsonify({
        "points": points,
        "center": center
    })

if __name__ == '__main__':
    app.run(debug=True, port=5000)
