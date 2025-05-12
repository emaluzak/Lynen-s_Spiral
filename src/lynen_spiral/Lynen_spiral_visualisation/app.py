import streamlit as st
from enhanced_fatty_acid import EnhancedFattyAcidMetabolism
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import sys
import rdkit
import os
import plotly.graph_objects as go
import numpy as np
import re
print(f"RDKit version: {rdkit.__version__}")

def mol_to_svg_image(mol, width=600, height=300):
        drawer = Draw.MolDraw2DSVG(width, height)
        drawer.drawOptions().bondLineWidth = 2.0
        drawer.drawOptions().atomLabelFontSize = 18
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace("svg:", "")  # Streamlit needs clean XML
        return svg

def smiles_to_delta(smiles):
    """
    Convert a SMILES string to its delta nomenclature.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    # Count the number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')

    # Extract all double bonds in the molecule
    double_bonds = []
    num_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            # Get the atoms involved in the bond (start and end)
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            
            # Only consider bonds between carbon atoms (fatty acids)
            if begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'C':
                num_double_bonds += 1
                double_bonds.append(bond)

    # Sort double bonds by the position of the first atom in the bond
    double_bonds = sorted(double_bonds, key=lambda x: min(x.GetBeginAtomIdx(), x.GetEndAtomIdx()))

    # Generate delta nomenclature
    delta_notation = []
    configuration = []
    for bond in double_bonds:
        start = bond.GetBeginAtomIdx()
        # The position of the double bond is the index of the start atom + 1 (1-based)
        position = start + 1
        # Append the delta notation for each double bond
        delta_notation.append(f"Δ{position}")

        stereo = bond.GetStereo()
        if stereo == Chem.BondStereo.STEREOZ:
            conf = "cis"
        elif stereo == Chem.BondStereo.STEREOE:
            conf = "trans"
        else:
            conf = ""
        configuration.append(conf)
    

    # Return the delta nomenclature as a string
    return f"{num_carbons}:{num_double_bonds}" + ("(" + ",".join(configuration) + "-" if configuration else "") + (",".join(delta_notation) + ")" if delta_notation else "")
    
def display_2d_structure(mol, smiles):
    st.subheader("2D Structure (SVG)")
    svg = mol_to_svg_image(mol)
    st.image(svg, caption=f"SMILES: {Chem.MolToSmiles(mol, canonical=True)}", use_container_width=True)

    col1, col2, col3 = st.columns(3)
    with col1:
        st.markdown("**Molecular Formula:**")
        st.code(Chem.rdMolDescriptors.CalcMolFormula(mol))
    with col2:
        st.markdown("**SMILES Notation:**")
        st.code(smiles)
    with col3:
        st.markdown("**Delta-Nomenclature:**")
        st.code(smiles_to_delta(smiles))

def display_3d_structure(mol):
    st.subheader("3D Structure")
    mol_with_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_with_h)
    AllChem.MMFFOptimizeMolecule(mol_with_h)
    conf = mol_with_h.GetConformer()
    coords = conf.GetPositions()

    fig = go.Figure()

    atom_colors = {'C': '#808080', 'O': '#FF0000', 'H': '#D3D3D3'}
    atom_sizes = {'C': 12, 'O': 14, 'H': 6}

    for i, atom in enumerate(mol_with_h.GetAtoms()):
        symbol = atom.GetSymbol()
        fig.add_trace(go.Scatter3d(
            x=[coords[i][0]], y=[coords[i][1]], z=[coords[i][2]],
            mode='markers',
            marker=dict(size=atom_sizes.get(symbol, 10),
                        color=atom_colors.get(symbol, '#00FFFF'),
                        opacity=0.9,
                        line=dict(width=0.5, color='black')),
            name=symbol,
            hoverinfo='text',
            hovertext=f"{symbol} (Atom {i+1})"
        ))

    bond_styles = {
        Chem.BondType.SINGLE: dict(width=5, color='#A0A0A0'),
        Chem.BondType.DOUBLE: dict(width=3, color='#FFA500'),
        Chem.BondType.TRIPLE: dict(width=2, color='#FF0000'),
        Chem.BondType.AROMATIC: dict(width=4, color='#00FF00')
    }

    for bond in mol_with_h.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bond_type = bond.GetBondType()
        if bond_type == Chem.BondType.DOUBLE:
            vec = coords[end] - coords[start]
            perp = np.cross(vec, [0, 0, 1])
            perp = perp / np.linalg.norm(perp) * 0.1
            for offset in [-perp, perp]:
                fig.add_trace(go.Scatter3d(
                    x=[coords[start][0] + offset[0], coords[end][0] + offset[0]],
                    y=[coords[start][1] + offset[1], coords[end][1] + offset[1]],
                    z=[coords[start][2] + offset[2], coords[end][2] + offset[2]],
                    mode='lines',
                    line=dict(width=3, color='#FFA500'),
                    hoverinfo='none',
                    showlegend=False
                ))
        else:
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

    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False),
            aspectmode='data'
        ),
        margin=dict(l=0, r=0, b=0, t=30),
        height=600
    )
    st.plotly_chart(fig, use_container_width=True)


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

    # Validate the input format and other constraints
    def validate_double_bond_input(input_str, num_carbons):
        invalid_reasons = []
    
        # Check if input contains anything but numbers, commas, or spaces
        if not re.match(r'^[0-9,\s]*$', input_str):
            invalid_reasons.append("Invalid character detected. Only numbers, commas, and spaces are allowed.")
    
        # Check if the input contains any negative signs
        if "-" in input_str:
            invalid_reasons.append("Position cannot be negative.")
    
        # Parse positions
        double_bond_positions = [int(p.strip()) for p in input_str.split(",") if p.strip().isdigit()]
    
        # Check if positions are within valid range
        for pos in double_bond_positions:
            if pos < 1:
                invalid_reasons.append(f"Position {pos} is invalid — it cannot be less than 1.")
            if pos > num_carbons:
                invalid_reasons.append(f"Position {pos} is invalid — it cannot be greater than {num_carbons}.")
    
        # 1. No double bond at position 1 (carboxylic acid carbon)
        if 1 in double_bond_positions:
            invalid_reasons.append("Position 1 is invalid — it's the carboxylic acid carbon.")
    
        # 2. No adjacent double bonds
        sorted_positions = sorted(double_bond_positions)
        for i in range(1, len(sorted_positions)):
            if sorted_positions[i] == sorted_positions[i-1] or sorted_positions[i] == sorted_positions[i-1] + 1:
                invalid_reasons.append(f"Double bonds at positions {sorted_positions[i-1]} and {sorted_positions[i]} are adjacent or duplicate.")
    
        # 3. No double bond at the terminal carbon (last bond is between C[n-1] and C[n])
        if num_carbons in double_bond_positions:
            invalid_reasons.append(f"Position {num_carbons} is invalid — it's the terminal methyl group.")
    
        return invalid_reasons, double_bond_positions

    # Step 1: Choose number of carbon atoms
    num_carbons = st.slider("Number of carbon atoms (including carboxylic acid carbon at position 1)", 4, 24, 16)

    # Step 2: Enter double bond positions
    positions_input = st.text_input(
        "Enter double bond positions (Δ notation, e.g., 9,12 for Δ9,12)",
        value=""
    )

    # Validate the input string
    invalid_reasons, double_bond_positions = validate_double_bond_input(positions_input, num_carbons)

    # Display warnings if there are any issues
    if invalid_reasons:
        st.warning("⚠️ Invalid double bond positions:\n\n" + "\n- " + "\n- ".join(invalid_reasons))
        double_bond_positions = []  # Clear the list to prevent invalid positions

    # Step 3: Configuration
    if not invalid_reasons:
        cis_trans_config = {}
        if double_bond_positions:
            st.subheader("Double Bond Configuration")
            for pos in double_bond_positions:
                cis_trans_config[pos] = st.selectbox(
                    f"Bond at position Δ{pos} (between C{pos}–C{pos+1})",
                    ["cis", "trans"],
                    key=f"conf_{pos}"
                )
        print(cis_trans_config)

        # Build SMILES from positions and config
        def build_fatty_acid_smiles(length, double_bonds, config):
            """
            Build a SMILES string for a fatty acid based on the number of carbons, double bond positions, and configuration.
            Cis: F/C=C\\F or F\\C=C/F       F's are on the same side of the double bond
            Trans: F/C=C/F or F\\C=C\\F     F's are on opposite sides of the double bond

            Therefore, the cis or trans configuration is based on whether there's a slash or backslash in the SMILES string for
            previous double bonds.

            Arguments:
            length (int): Number of carbon atoms in the fatty acid chain.
            double_bonds (list): List of positions for double bonds.
            config (dict): Dictionary with double bond positions as keys and their configurations ("cis" or "trans")

            Returns:
            str: SMILES representation of the fatty acid.
            """
            if length < 2:
                raise ValueError("Chain must be at least 2 carbons long.")
    
            smiles = "OC(=O)"

            for i in range(1, length):  # C2 to Cn
                print(smiles)
                if (i) in double_bonds:  

                    # Check if i-2 position carbon has a double bond (important for SMILES generation)
                    
                    if i != (1,2) and i-2 in double_bonds:

                        if smiles[-2] == "\\":

                            if config.get(i) == "cis":
                                smiles = smiles[:-1] + "C=C/"
                            else:
                                smiles = smiles[:-1] + "C=C\\"
                                
                        if smiles[-2] == "/":
                            if config.get(i) == "cis":
                                smiles = smiles[:-1] + "C=C\\"
                            else:
                                smiles = smiles[:-1] + "C=C/"

                    elif config.get(i) == "cis":
                        smiles = smiles[:-1] + "/C=C\\"

                    elif config.get(i) == "trans":
                        smiles = smiles[:-1] + "/C=C/"
                    
                else:
                    smiles += "C"

            # Remove slashes at the end if they exist (Causes error)
            if smiles[-1] in ["/", "\\"]:
                smiles = smiles[:-1] 
                
            return smiles
        
        smiles = build_fatty_acid_smiles(num_carbons, double_bond_positions, cis_trans_config)
        notation = f"{num_carbons}:{len(double_bond_positions)}" + (
            f"({','.join(map(str, double_bond_positions))})" if double_bond_positions else ""
        )

if st.button("Visualize Fatty Acid"):
    if 'smiles' in locals() and smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Invalid SMILES string.")
        
            display_2d_structure(mol, smiles)
            display_3d_structure(mol)

            # Store session state
            fa_metabolism = EnhancedFattyAcidMetabolism(smiles)
            st.session_state['fa_data'] = {
                'smiles': smiles,
                'notation': notation,
                'molecule': mol,
                'processor': fa_metabolism,
                'name': selected_fa if option == "Choose from library" else f"Custom {notation}"
            }

            # Update session state for the next step (this controls page navigation)
            st.session_state['page'] = '2_🌀_3D_Spiral'

            # Display a button for going to the 3D Spiral page
            st.button("Proceed to Detailed 3D Spiral Visualization", on_click=lambda: switch_page("2_🌀_3D_Spiral"))

        except Exception as e:
            st.error(f"Visualization error: {e}")
    else:
        st.warning("⚠️ Please select or generate a valid fatty acid structure before visualizing.")

# Page switcher function using session state
def switch_page(page_name: str):
    st.session_state['page'] = page_name
    st.rerun()  # Force rerun to switch content

