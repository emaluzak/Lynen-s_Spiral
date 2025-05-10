import streamlit as st
from enhanced_fatty_acid import EnhancedFattyAcidMetabolism
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import sys
import rdkit
import os
import plotly.graph_objects as go
import numpy as np
print(f"RDKit version: {rdkit.__version__}")

st.title("Fatty Acid Metabolism Visualizer")

# Correct fatty acid library with specific double bond positions
FATTY_ACIDS = {
    "Palmitic acid (16:0)": "CCCCCCCCCCCCCCCC(=O)O",
    "Heptadecanoic acid (17:0)": "CCCCCCCCCCCCCCCCC(=O)O",
    "15-Pentadecenoic acid (C15:1(10))": "CCCCCCCCCCC=CCCC(=O)O",
    "Stearic acid (18:0)": "CCCCCCCCCCCCCCCCCC(=O)O",
    "Oleic acid (18:1(9))": "CCCCCCCC\C=C/CCCCCCCC(O)=O",  # Δ9
    "Linoleic acid (18:2(9,12))": "CCCCC/C=C\C/C=C\CCCCCCCC(=O)O",  # Δ9,12
    "α-Linolenic acid (18:3(9,12,15))": "CC/C=C\C/C=C\C/C=C\CCCCCCCC(=O)O",  # Δ9,12,15
    "Arachidonic acid (20:4(5,8,11,14))": "CCCCC/C=C\C/C=C\C/C=C\C/C=C\CCCC(=O)O"  # Δ5,8,11,14
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
        
        # ====================== 2D Visualization ======================
        st.subheader("2D Structure")
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
            if hasattr(fa_metabolism, 'get_systematic_name'):
                st.markdown("**Systematic Name:**")
                st.code(fa_metabolism.get_systematic_name())

        # In your visualization section (after the 2D view):
        # ====================== 3D Visualization ======================
        st.subheader("3D Structure")

        # Generate 3D coordinates
        mol_with_h = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_with_h)
        AllChem.MMFFOptimizeMolecule(mol_with_h)

        # Extract coordinates
        conf = mol_with_h.GetConformer()
        coords = conf.GetPositions()

        # Create Plotly figure
        fig = go.Figure()

        # Atom styling
        atom_colors = {
            'C': '#808080',  # Gray
            'O': '#FF0000',  # Red
            'H': '#D3D3D3',  # Light gray
            }
        atom_sizes = {
            'C': 12,
            'O': 14,
            'H': 6
        }

        # Add atoms
        for i, atom in enumerate(mol_with_h.GetAtoms()):
            symbol = atom.GetSymbol()
            fig.add_trace(go.Scatter3d(
                x=[coords[i][0]],
                y=[coords[i][1]],
                z=[coords[i][2]],
                mode='markers',
                marker=dict(
                    size=atom_sizes.get(symbol, 10),
                    color=atom_colors.get(symbol, '#00FFFF'),
                    opacity=0.9,
                    line=dict(width=0.5, color='black')
                ),
                name=symbol,
                hoverinfo='text',
                hovertext=f"{symbol} (Atom {i+1})"
            ))

        # Bond styling - DIFFERENT STYLES FOR EACH BOND TYPE
        bond_styles = {
            Chem.BondType.SINGLE: dict(width=5, color='#A0A0A0'),  # Medium gray
            Chem.BondType.DOUBLE: dict(width=3, color='#FFA500'),   # Orange
            Chem.BondType.TRIPLE: dict(width=2, color='#FF0000'),   # Red
            Chem.BondType.AROMATIC: dict(width=4, color='#00FF00')  # Green
        }

        # Add bonds with different styles
        for bond in mol_with_h.GetBonds():
            start = bond.GetBeginAtomIdx()
            end = bond.GetEndAtomIdx()
            bond_type = bond.GetBondType()
    
            # For DOUBLE BONDS - draw two parallel lines
            if bond_type == Chem.BondType.DOUBLE:
                # Calculate perpendicular direction
                vec = coords[end] - coords[start]
                perp = np.cross(vec, [0, 0, 1])  # Arbitrary perpendicular vector
                perp = perp / np.linalg.norm(perp) * 0.1  # Normalize and scale
        
                # Draw two parallel lines
                for offset in [-perp, perp]:
                    fig.add_trace(go.Scatter3d(
                        x=[coords[start][0] + offset[0], coords[end][0] + offset[0]],
                        y=[coords[start][1] + offset[1], coords[end][1] + offset[1]],
                        z=[coords[start][2] + offset[2], coords[end][2] + offset[2]],
                        mode='lines',
                        line=dict(
                            width=3,
                            color='#FFA500'  # Orange
                        ),
                        hoverinfo='none',
                        showlegend=False
                    ))
            else:
                # Single, triple or aromatic bonds
                style = bond_styles.get(bond_type, bond_styles[Chem.BondType.SINGLE])
                fig.add_trace(go.Scatter3d(
                    x=[coords[start][0], coords[end][0]],
                    y=[coords[start][1], coords[end][1]],
                    z=[coords[start][2], coords[end][2]],
                    mode='lines',
                    line=style,
                    hoverinfo='none',
                    showlegend=False
                ))

        # Configure layout
        fig.update_layout(
            scene=dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False),
                aspectmode='data',
                camera=dict(
                    up=dict(x=0, y=0, z=1),
                    center=dict(x=0, y=0, z=0),
                    eye=dict(x=1.5, y=1.5, z=0.8)
                )
            ),
            margin=dict(l=0, r=0, b=0, t=30),
            height=600,
            showlegend=True,
            legend=dict(
                orientation="h",
                itemsizing='constant',
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1
            )
        )

        # Add bond type legend
        fig.add_annotation(
            x=0.05, y=0.95,
            xref="paper", yref="paper",
            text="<b>Bond Types:</b><br>"
                "<span style='color:#A0A0A0'>—— Single</span><br>"
                "<span style='color:#FFA500'>══ Double</span><br>"
                "<span style='color:#FF0000'>≡ Triple</span>",
            showarrow=False,
            bgcolor="white",
            bordercolor="black",
            borderwidth=1
        )

        st.plotly_chart(fig, use_container_width=True)
        
        # Store data for next pages
        st.session_state['fa_data'] = {
            'smiles': smiles,
            'notation': notation,
            'molecule': mol,
            'processor': fa_metabolism,
            'name': selected_fa if option == "Choose from library" else f"Custom {notation}"
        }

        #Navigation button
        if st.button("Proceed to Detailed 3D Spiral Visualization", key="go_to_3d"):
            # Force clear any existing components
            st.session_state['fa_data'] = {
                'smiles': smiles,
                'notation': notation,
                'molecule': mol,
                'processor': fa_metabolism,
                'name': selected_fa if option == "Choose from library" else f"Custom {notation}"
            }

            from streamlit.runtime.scriptrunner import RerunData, RerunException
            from streamlit.source_util import get_pages
            def switch_page(page_name: str):
                pages = get_pages("app.py")  # Gets all pages from main
                for page_hash, config in pages.items():
                    if config["page_name"] == page_name:
                        raise RerunException(
                            RerunData(
                                page_script_hash=page_hash,
                                page_name=page_name,
                            )
                        )
            switch_page("2_🌀_3D_Spiral")  # Exact page name as shown in sidebar
    except Exception as e:
        st.error(f"Error displaying fatty acid: {str(e)}")

