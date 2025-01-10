from re import match
from flask import Flask, render_template, send_file, jsonify
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
import pandas as pd
import numpy as np
import io

app = Flask(__name__)

@app.route("/", methods=["GET"])
def root():
    return render_template("base.html")

def generateFingerprint(smile):
    m = Chem.MolFromSmiles(smile)
    gen = rdFingerprintGenerator.GetRDKitFPGenerator()
    return gen.GetFingerprint(m)

def findSimilarity(smile1, smile2):
    return DataStructs.TanimotoSimilarity(smile1, smile2)

@app.route("/findSimilarStructures", methods=["GET"])
def findSimilarStructures(query="Cc1ccccc1", limit=0.5):
    df = pd.read_csv('data.csv', sep=",")
    smileList = df['smiles']
    queryFingerprint = generateFingerprint(query)
    matches = []
    for i in smileList:
        if findSimilarity(queryFingerprint, generateFingerprint(i)) > limit:
            matches.append(i)
    print(matches)
    return jsonify(matches)

query = "C1=CC=C(C=C1)C=O"
limit = 0.5
df = pd.read_csv('data.csv', sep=",")

smileList = df['smiles']
# print(smileList)
queryFingerprint = generateFingerprint(query)
print(np.array(queryFingerprint))
matches = []
print(findSimilarity(queryFingerprint, generateFingerprint("C1=CC=C(C(=C1)C=O)C=O")))

for i in smileList:
    if findSimilarity(queryFingerprint, generateFingerprint(i)) > limit:
        matches.append(i)
print(matches)

@app.route("/drawMolecule", methods=["GET"])
def drawMolecule():
    m = Chem.MolFromSmiles("Cc1ccccc1")
    img = Draw.MolToImage(m)
    img_io = io.BytesIO()
    img.save(img_io, "PNG")
    img_io.seek(0)
    return send_file(img_io, mimetype="image/png")
