import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
import py3Dmol

st.set_page_config(layout="wide")

st.title("2D Cycle Representation")

# Get data from previous pages
fa_data = st.session_state.get('fa_data', {})

# Display all steps
steps = {
    "Activation": "CCCCCCCCCCCCCCCC(=O)O>>CCCCCCCCCCCCCCC(=O)SCC",
    "First oxidation": "CCCCCCCCCCCCCCC(=O)SCC>>CCCCCCCCCCCCC=CC(=O)SCC",
    # Add all steps with example SMILES
}

selected_step = st.selectbox("Select a step:", list(steps.keys()))

# Display selected step
if selected_step:
    mol = Chem.MolFromSmiles(steps[selected_step].split(">>")[0])
    img = Draw.MolToImage(mol)
    st.image(img, caption=f"Reactant: {selected_step}", use_column_width=True)
    
    # Add mechanism details, ATP calculations, etc.
    
# Navigation buttons
col1, col2 = st.columns(2)
with col1:
    if st.button("Back to 3D Spiral"):
        st.switch_page("pages/2_🌀_3D_Spiral.py")
with col2:
    if st.button("Back to Home"):
        st.switch_page("app.py")