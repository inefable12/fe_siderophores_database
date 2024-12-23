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
mol2_files = [f"Fe_sideroforos_mol2{i}.mol2" for i in range(1, 331)]

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
        showmol(viewer, height=400, width=500) #j1
    else:
        st.error(f"No se pudo procesar el archivo {selected_file}. Verifica su contenido.")
else:
    st.error("No se pudo cargar el archivo seleccionado.")

show_3d(smiles)
######################################
######################################
######################################
import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem import rdmolfiles
from io import StringIO

# Función para leer el archivo MOL2 desde GitHub
def fetch_mol2_from_github(url):
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        return None

# Función para convertir MOL2 a SMILES
def mol2_to_smiles(mol2_data):
    mol = rdmolfiles.MolFromMol2Block(mol2_data)
    if mol:
        return Chem.MolToSmiles(mol)
    else:
        return None

# Título de la aplicación
st.title('Conversión de MOL2 a SMILES')

# Ingreso de la URL del archivo MOL2 desde GitHub
##url = st.text_input('Ingresa la URL del archivo MOL2 en GitHub:', 'https://github.com/usuario/repositorio/archivo.mol2')

# Botón para convertir el archivo MOL2 a SMILES
if st.button('Convertir a SMILES'):
    mol2_data = fetch_mol2_from_github(mol2_content)
    
    if mol2_data:
        smiles = mol2_to_smiles(mol2_data)
        
        if smiles:
            st.write(f"El SMILES generado es: {smiles}")
        else:
            st.error("No se pudo convertir el archivo MOL2 a SMILES.")
    else:
        st.error("No se pudo descargar el archivo MOL2 desde GitHub.")
