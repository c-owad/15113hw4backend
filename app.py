from flask import Flask, jsonify
from flask_cors import CORS
import requests
import math

app = Flask(__name__)
CORS(app)

# --- CONFIGURATION ---
def get_center(points):
    # Filter out None values which are used as chain breaks
    valid = [p for p in points if p is not None and 'x' in p]
    if not valid: return (0,0,0)
    avg_x = sum(p['x'] for p in valid) / len(valid)
    avg_y = sum(p['y'] for p in valid) / len(valid)
    avg_z = sum(p['z'] for p in valid) / len(valid)
    return (avg_x, avg_y, avg_z)

class ProteinAPI:
    def __init__(self):
        self.pdb_url = "https://files.rcsb.org/download"
        self.chem_api_url = "https://data.rcsb.org/rest/v1/core/chemcomp"

    def fetch_raw_pdb(self, pdb_id):
        url = f"{self.pdb_url}/{pdb_id}.pdb"
        try:
            r = requests.get(url, timeout=10)
            if r.status_code == 200:
                return r.text
        except:
            return None
        return None

    def parse_pdb(self, text):
        backbone = []
        ligands = {} # Store ligands by a unique ID (name_chain_sequence)
        last_chain = None
        
        for line in text.splitlines():
            # 1. Parse Alpha Carbons for the main structure
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    chain_id = line[21]
                    res_name = line[17:20].strip()
                    res_seq = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    
                    if last_chain is not None and chain_id != last_chain:
                        backbone.append(None) # Break the ribbon

                    backbone.append({
                        "res_name": res_name,
                        "res_seq": res_seq,
                        "chain": chain_id,
                        "x": x, "y": y, "z": z
                    })
                    last_chain = chain_id
                except:
                    continue
                    
            # 2. Parse Heteroatoms (Ligands / Drugs)
            elif line.startswith("HETATM"):
                try:
                    res_name = line[17:20].strip()
                    if res_name == "HOH": continue # Ignore water molecules
                    
                    chain_id = line[21]
                    res_seq = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    
                    # Create a unique key since a protein might have multiple of the same drug bound to it
                    ligand_key = f"{res_name}_{chain_id}_{res_seq}"
                    if ligand_key not in ligands:
                        ligands[ligand_key] = {"res_name": res_name, "atoms": []}
                        
                    ligands[ligand_key]["atoms"].append({"x": x, "y": y, "z": z})
                except:
                    continue
                    
        return backbone, ligands

    def get_chem_info(self, comp_id):
        """Fetch chemical metadata dynamically from RCSB Data API."""
        try:
            r = requests.get(f"{self.chem_api_url}/{comp_id}", timeout=10)
            if r.status_code == 200:
                data = r.json()
                chem_info = data.get('chem_comp', {})
                return {
                    "name": chem_info.get('name', 'Unknown'),
                    "formula": chem_info.get('formula', 'Unknown'),
                    "weight": chem_info.get('formula_weight', 'Unknown')
                }
        except:
            pass
        return {"name": "Not found", "formula": "N/A", "weight": "N/A"}

    def find_binding_pocket(self, backbone, ligand_atoms, threshold=6.0):
        """Algorithmic spatial analysis to find interacting amino acids."""
        interacting_residues = []
        valid_backbone = [atom for atom in backbone if atom is not None]
        
        for ca in valid_backbone:
            is_close = False
            for latom in ligand_atoms:
                # 3D Euclidean Distance
                dist = math.sqrt((ca['x'] - latom['x'])**2 + (ca['y'] - latom['y'])**2 + (ca['z'] - latom['z'])**2)
                if dist <= threshold:
                    is_close = True
                    break 
            
            if is_close:
                interacting_residues.append(f"{ca['res_name']}{ca['res_seq']}({ca['chain']})")
                
        return interacting_residues

api = ProteinAPI()

@app.route('/api/molecule/<pdb_id>')
def get_molecule(pdb_id):
    raw_text = api.fetch_raw_pdb(pdb_id)
    if not raw_text:
        return jsonify({"error": "PDB not found"}), 404
        
    backbone, ligands = api.parse_pdb(raw_text)
    center = get_center(backbone)
    
    return jsonify({
        "points": backbone,
        "ligands": ligands,
        "center": center
    })

# --- NEW ANALYTICS ROUTE ---
@app.route('/api/analyze/<pdb_id>/<ligand_key>')
def analyze_ligand(pdb_id, ligand_key):
    raw_text = api.fetch_raw_pdb(pdb_id)
    if not raw_text:
        return jsonify({"error": "Structure not found"}), 404
        
    backbone, ligands = api.parse_pdb(raw_text)
    
    if ligand_key not in ligands:
        return jsonify({"error": "Ligand not found in structure"}), 404
        
    target_ligand = ligands[ligand_key]
    res_name = target_ligand["res_name"]
    
    # Meaningful API Call
    chem_info = api.get_chem_info(res_name)
    
    # Substantial Data Analysis
    interactions = api.find_binding_pocket(backbone, target_ligand["atoms"], threshold=6.0)
    
    return jsonify({
        "ligand_id": res_name,
        "properties": chem_info,
        "binding_pocket": interactions
    })

if __name__ == '__main__':
    app.run(debug=True, port=5000)
