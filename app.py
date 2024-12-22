import streamlit as st
from rdkit import Chem
from rdkit.Chem.rdmolfiles import MolFromMol2File
from stmol import showmol

# Título de la aplicación
st.title("Visualizador de Moléculas 3D")
st.markdown("Autor: Jesus Alvarado-Huayhuaz")

# Ruta del archivo mol2
url = "https://raw.githubusercontent.com/inefable12/fe_siderophores_database/main/complexes/Fe_sideroforos_mol21.mol2"
mol = MolFromMol2File(url, sanitize=True)
mol_block = Chem.MolToMolBlock(mol)            
showmol(mol_block, width=800, height=400)
