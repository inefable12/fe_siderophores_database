import streamlit as st
from rdkit import Chem
from rdkit.Chem.rdmolfiles import MolFromMol2File
from stmol import showmol  # Py3Dmol wrapper for Streamlit

# Título de la aplicación
st.title("Visualizador de Moléculas 3D")
st.markdown("Esta aplicación permite visualizar un archivo `.mol2` desde un repositorio en GitHub.")

# Ruta del archivo mol2
github_url = "https://raw.githubusercontent.com/<usuario>/<repositorio>/<rama>/<ruta_del_archivo>.mol2"

# Función para cargar y visualizar el archivo mol2
def visualize_molecule_from_url(url):
    try:
        # Descargar el archivo
        mol = MolFromMol2File(url, sanitize=True)
        if mol:
            # Convertir la molécula a formato MOL block
            mol_block = Chem.MolToMolBlock(mol)
            # Mostrar la molécula en 3D
            st.subheader("Visualización 3D:")
            showmol(mol_block, width=800, height=400)
        else:
            st.error("No se pudo procesar el archivo. Verifica que sea un archivo `.mol2` válido.")
    except Exception as e:
        st.error(f"Ocurrió un error al
