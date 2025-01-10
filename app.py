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


@app.route("/findSimilarStructures", methods=["GET"])
def findSimilarity(fp1, fp2):
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def findSimilarCompounds(query, limit=0.50):
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
            similarityScaffold = findSimilarity(
                queryFingerprint, FingerprintMols.FingerprintMol(scaffold)
            )

            if similarityScaffold == 1.0:
                gen = rdFingerprintGenerator.GetRDKitFPGenerator()
                similaritySmile = findSimilarity(
                    gen.GetFingerprint(queryMolecule), gen.GetFingerprint(mol)
                )
                if similaritySmile >= limit:
                    matches.append((smile, similaritySmile))
        except:
            continue
    matches.sort(key=lambda x: x[1], reverse=True)
    return matches

@app.route("/drawMolecule", methods=["GET"])
def drawMolecule(query):
    m = Chem.MolFromSmiles(query)
    img = Draw.MolToImage(m)
    img_io = io.BytesIO()
    img.save(img_io, "PNG")
    img_io.seek(0)
    return send_file(img_io, mimetype="image/png")
