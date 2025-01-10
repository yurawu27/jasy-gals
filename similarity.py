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


def smiles_to_svg(smiles: str, width: int = 400, height: int = 400) -> bytes:
    """
    makes an SVG image of a molecule
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError("Invalid SMILES")

    Chem.rdCoordGen.AddCoords(mol)
    drawer = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
    # set drawing options on drawer.getOptions()
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText().encode()
