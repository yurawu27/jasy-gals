from re import match
from flask import Flask, render_template, send_file, jsonify
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
import pandas as pd
import numpy as np
import io

app = Flask(__name__)


@app.route("/", methods=["GET"])
def root():
    return render_template("base.html")

def findSimilarity(fp1, fp2):
    return DataStructs.TanimotoSimilarity(fp1, fp2)

@app.route("/findSimilarStructures", methods=["GET"])
def findSimilar(query="C1=CC=C(C=C1)C=O", limit=0.85):
    df = pd.read_csv("data.csv", sep=",")
    smileList = df["smiles"]

    queryMolecule = Chem.MolFromSmiles(query)
    queryScaffold = MurckoScaffold.GetScaffoldForMol(queryMolecule)
    queryFingerprint = FingerprintMols.FingerprintMol(queryScaffold)
    matches = []
    for smile in smileList:
        try:
            mol = Chem.MolFromSmiles(smile)
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            if (
                findSimilarity(
                    queryFingerprint,
                    FingerprintMols.FingerprintMol(scaffold),
                )
                >= limit
            ):
                matches.append(smile)
        except:
            print("Invalid")

    return matches

print(len(findSimilar()))

@app.route("/drawMolecule", methods=["GET"])
def drawMolecule():
    m = Chem.MolFromSmiles("Cc1ccccc1")
    img = Draw.MolToImage(m)
    img_io = io.BytesIO()
    img.save(img_io, "PNG")
    img_io.seek(0)
    return send_file(img_io, mimetype="image/png")
