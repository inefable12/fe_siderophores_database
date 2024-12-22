import streamlit as st
from rdkit import Chem
from rdkit.Chem.rdmolfiles import MolFromMol2File
from stmol import showmol
from rdkit.Chem import SanitizeMol

# Título de la aplicación
st.title("Visualizador de Moléculas 3D")
st.markdown("Autor: Jesus Alvarado-Huayhuaz")

from rdkit import Chem
from rdkit.Chem import rdmolfiles

# Ruta al archivo MOL2
mol2_file = "https://raw.githubusercontent.com/inefable12/fe_siderophores_database/main/complexes/Fe_sideroforos_mol21.mol2"

# Leer el archivo MOL2
mol = rdmolfiles.MolFromMol2File(mol2_file)

if mol:
    # Convertir a SMILES
    smiles = Chem.MolToSmiles(mol)
    print(f"SMILES: {smiles}")
else:
    print("No se pudo leer el archivo MOL2.")
