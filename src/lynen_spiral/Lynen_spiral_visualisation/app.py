import streamlit as st
from enhanced_fatty_acid import EnhancedFattyAcidMetabolism
from rdkit import Chem
from rdkit.Chem import Draw
import sys
import os
import rdkit
print(f"RDKit version: {rdkit.__version__}")

st.title("Fatty Acid Metabolism Visualizer")

# Correct fatty acid library with specific double bond positions
FATTY_ACIDS = {
    "Palmitic acid (16:0)": "CCCCCCCCCCCCCCCC(=O)O",
    "Stearic acid (18:0)": "CCCCCCCCCCCCCCCCCC(=O)O",
    "Oleic acid (18:1(9))": "CCCCCCCC=CCCCCCCCC(=O)O",  # Δ9
    "Linoleic acid (18:2(9,12))": "CCCCCC=CCC=CCCCCCCC(=O)O",  # Δ9,12
    "α-Linolenic acid (18:3(9,12,15))": "CCC=CCC=CCC=CCCCCCCC(=O)O",  # Δ9,12,15
    "Arachidonic acid (20:4(5,8,11,14))": "CCCC=CCC=CCC=CCC=CCCCCC(=O)O"  # Δ5,8,11,14
}

# Fatty acid selection
option = st.radio("Select fatty acid source:", 
                 ("Choose from library", "Customize your own"))

if option == "Choose from library":
    selected_fa = st.selectbox("Select fatty acid:", list(FATTY_ACIDS.keys()))
    smiles = FATTY_ACIDS[selected_fa]
    notation = selected_fa.split("(")[1].replace(")", "") if "(" in selected_fa else None
else:
    # Custom fatty acid builder with proper double bond positioning
    st.subheader("Custom Fatty Acid Builder")
    col1, col2 = st.columns(2)
    with col1:
        length = st.slider("Carbon chain length:", 4, 24, 16)
    with col2:
        unsaturations = st.slider("Number of double bonds:", 0, min(6, length//2), 0)
    
    # Default positions (modify if you want user to specify)
    if unsaturations > 0:
        # Common positions for natural fatty acids
        if length == 18 and unsaturations == 2:
            positions = "9,12"  # Linoleic acid pattern
        elif length == 18 and unsaturations == 3:
            positions = "9,12,15"  # α-Linolenic acid pattern
        else:
            # Default to spaced double bonds
            positions = ",".join(str(9 + 3*i) for i in range(unsaturations))
        notation = f"C{length}:{unsaturations}({positions})"
    else:
        notation = f"C{length}:0"
    
    # Create instance to convert notation to SMILES
    try:
        fa_processor = EnhancedFattyAcidMetabolism(notation)
        smiles = Chem.MolToSmiles(fa_processor.molecule)
    except Exception as e:
        st.error(f"Error creating custom fatty acid: {str(e)}")
        st.stop()

# Visual representation
if st.button("Show fatty acid structure"):
    try:
        # Create EnhancedFattyAcidMetabolism instance
        fa_metabolism = EnhancedFattyAcidMetabolism(smiles)
        mol = fa_metabolism.molecule
        
        # Generate visualization
        img = Draw.MolToImage(mol, size=(600, 300), kekulize=True)
        st.image(img, use_container_width=True)
        
        # Show molecular information
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("**Molecular Formula:**")
            st.code(Chem.rdMolDescriptors.CalcMolFormula(mol))
        with col2:
            st.markdown("**SMILES Notation:**")
            st.code(smiles)
            
            # Systematic name display
            try:
                systematic_name = fa_metabolism.get_systematic_name()
                st.markdown("**Systematic Name:**")
                st.code(systematic_name)
            except Exception as e:
                st.warning(f"Couldn't generate systematic name: {str(e)}")
        
        # Store data for next pages
        st.session_state['fa_data'] = {
            'smiles': smiles,
            'notation': notation,
            'molecule': mol,
            'processor': fa_metabolism,
            'name': selected_fa if option == "Choose from library" else f"Custom {notation}"
        }
        
        if st.button("Proceed to 3D Spiral Visualization"):
            st.switch_page("pages/2_🌀_3D_Spiral.py")
            
    except Exception as e:
        st.error(f"Error displaying fatty acid: {str(e)}")