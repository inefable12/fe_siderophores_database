import streamlit as st
import py3Dmol
import requests
from rdkit import Chem
from rdkit.Chem.rdmolfiles import MolToMolBlock, MolFromMol2Block

# Título de la aplicación
st.title("Complejos Fe-sideróforos 3D")
st.markdown("Autor: Jesus Alvarado-Huayhuaz")

# Ruta base del repositorio en GitHub
base_url = "https://raw.githubusercontent.com/inefable12/fe_siderophores_database/main/complexes/"

# Generar lista de nombres de archivos
mol2_files = [f"Fe_sideroforos_mol{i}.mol2" for i in range(1, 201)]

# Slider para seleccionar el archivo
selected_index = st.slider("Selecciona una molécula", 1, len(mol2_files))
selected_file = mol2_files[selected_index - 1]

# URL del archivo seleccionado
selected_url = f"{base_url}{selected_file}"

# Función para cargar el contenido del archivo directamente desde GitHub
def fetch_mol2_content(url):
    try:
        response = requests.get(url)
        response.raise_for_status()  # Asegúrate de que no hay errores
        return response.text  # Devolver el contenido como texto
    except Exception as e:
        st.error(f"Error al leer el archivo desde GitHub: {e}")
        return None

# Función para visualizar molécula con py3Dmol
def visualize_molecule(mol_block):
    viewer = py3Dmol.view(width=800, height=400)
    viewer.addModel(mol_block, "mol")
    viewer.setStyle({"stick": {}})
    viewer.zoomTo()
    return viewer

# Leer el archivo .mol2 desde GitHub
st.info(f"Cargando la molécula: {selected_file}")
mol2_content = fetch_mol2_content(selected_url)

if mol2_content:
    # Procesar el contenido del archivo con RDKit
    mol = MolFromMol2Block(mol2_content, sanitize=True)
    if mol:
        # Convertir a formato MOL block para py3Dmol
        mol_block = MolToMolBlock(mol)
        # Mostrar la molécula en 3D
        st.subheader(f"Visualización 3D de: {selected_file}")
        viewer = visualize_molecule(mol_block)
        viewer.show()
    else:
        st.error(f"No se pudo procesar el archivo {selected_file}. Verifica su contenido.")
else:
    st.error("No se pudo cargar el archivo seleccionado.")
