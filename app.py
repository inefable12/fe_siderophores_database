import streamlit as st
import py3Dmol
import requests
from rdkit import Chem
from rdkit.Chem.rdmolfiles import MolToMolBlock, MolFromMol2Block

# Título de la aplicación
st.title("Complejos Fe-sideróforos 3D")
st.markdown("Autor: Jesus Alvarado-Huayhuaz")

# URL del archivo .mol2 en GitHub
github_url = "https://raw.githubusercontent.com/inefable12/fe_siderophores_database/main/complexes/Fe_sideroforos_mol21.mol2"

# Función para cargar el archivo desde GitHub
def fetch_mol2_file(url):
    try:
        response = requests.get(url)
        response.raise_for_status()  # Asegurarse de que no hay errores de descarga
        return response.text
    except Exception as e:
        st.error(f"Error al descargar el archivo: {e}")
        return None

# Función para visualizar molécula con py3Dmol
def visualize_molecule(mol_block):
    viewer = py3Dmol.view(width=800, height=400)
    viewer.addModel(mol_block, "mol")
    viewer.setStyle({"stick": {}})
    viewer.zoomTo()
    return viewer

# Descargar el archivo .mol2
st.info("Cargando la molécula desde GitHub...")
mol2_content = fetch_mol2_file(github_url)

if mol2_content:
    # Procesar el archivo con RDKit
    mol = MolFromMol2Block(mol2_content, sanitize=True)
    if mol:
        # Convertir a formato MOL block para py3Dmol
        mol_block = MolToMolBlock(mol)
        # Mostrar la molécula en 3D
        st.subheader("Visualización 3D:")
        viewer = visualize_molecule(mol_block)
        viewer.show()
    else:
        st.error("No se pudo procesar el archivo .mol2. Verifica su contenido.")
