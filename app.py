import os
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolfiles import MolFromMol2File
from stmol import showmol  # Py3Dmol wrapper for Streamlit

# Título de la aplicación
st.title("Visualizador de Moléculas 3D")
st.markdown("Selecciona un archivo `.mol2` de la carpeta 'complexes' para visualizar su estructura en 3D.")

# Ruta de la carpeta de moléculas
mol2_folder = "complexes"

# Verificar si la carpeta existe
if not os.path.exists(mol2_folder):
    st.error(f"La carpeta '{mol2_folder}' no existe. Asegúrate de que esté en el mismo nivel que este script.")
else:
    # Listar los archivos mol2
    mol2_files = [f for f in os.listdir(mol2_folder) if f.endswith(".mol2")]

    if len(mol2_files) == 0:
        st.warning("No se encontraron archivos `.mol2` en la carpeta 'complexes'.")
    else:
        # Seleccionar un archivo
        selected_file = st.selectbox("Selecciona un archivo:", mol2_files)

        # Cargar y procesar la molécula
        mol_path = os.path.join(mol2_folder, selected_file)
        try:
            mol = MolFromMol2File(mol_path, sanitize=True)

            if mol:
                # Generar coordenadas 3D si no están disponibles
                if not mol.GetConformer().Is3D():
                    AllChem.EmbedMolecule(mol)

                # Convertir a formato para py3Dmol
                mol_block = Chem.MolToMolBlock(mol)
                
                # Visualizar molécula
                st.subheader("Visualización 3D de la molécula:")
                showmol(mol_block, width=800, height=400)

            else:
                st.error(f"No se pudo procesar el archivo '{selected_file}'. Verifica su formato.")
        except Exception as e:
            st.error(f"Ocurrió un error al cargar el archivo: {e}")
