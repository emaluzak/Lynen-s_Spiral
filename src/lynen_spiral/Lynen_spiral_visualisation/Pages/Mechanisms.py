import streamlit as st
import io
import sys
import os
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image

st.set_page_config(layout="wide")
st.title("Detailed Analysis")

sys.path.append(os.path.abspath(r"C:\Users\Lenovo\git\Lynen-s_Spiral\src\lynen_spiral\Lynen_spiral_visualisation"))
from enhanced_fatty_acid import EnhancedFattyAcidMetabolism

# ------------------ Useful Functions ------------------

def mol_to_svg_image(mol, width=600, height=300):
        
    """
    Generate an SVG image of a molecule using RDKit's 2D drawing tools.
    Higher resolution and larger size for better visualization compared 
    to RDKit's MolToImage function.

    Parameters:
        mol (rdkit.Chem.Mol): The RDKit molecule object to be drawn.
        width (int, optional): Width of the SVG image in pixels. Defaults to 600.
        height (int, optional): Height of the SVG image in pixels. Defaults to 300.

    Returns:
        str: An SVG string representation of the molecule, with cleaned XML tags
             for compatibility with Streamlit rendering.
    """

    drawer = Draw.MolDraw2DSVG(width, height)
    drawer.drawOptions().bondLineWidth = 2.0
    drawer.drawOptions().atomLabelFontSize = 18
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg

def _display_cycle_steps(cycle):

    """
    Display detailed information and molecular images for each step within a given β-oxidation cycle.

    This helper function uses Streamlit to present:
    - Reactant and product molecular structures with CoA group visualization.
    - Step-specific descriptive information including the reaction type, and expected ATP yield.
    - SMARTS pattern used to identify the reaction.
    - The chemical transformation from reactant to product in SMILES notation.

    Parameters:
        cycle (dict): A dictionary representing a single cycle, expected to have a 'steps' key 
                      containing a list of step dictionaries. Each step dictionary should include:
                        - 'step' (str): Name/type of the biochemical step (e.g., "Dehydrogenation").
                        - 'input' (str): SMILES string of the reactant molecule.
                        - 'output' (str): SMILES string of the product molecule.
                        - 'smarts' (str): SMARTS pattern string describing the chemical transformation.

    Note:
        This function depends on the mol_to_svg_image() function for molecular rendering.
    """

    for step in cycle['steps']:
        col1, col2 = st.columns([1, 2])

        with col1:

            st.markdown("**Reactant**")

            # Adding CoA group to the reactant SMILES

            reactant_smiles = step['input']
            reactant_smiles = reactant_smiles + "*"
            reactant_mol = Chem.MolFromSmiles(reactant_smiles)
            for atom in reactant_mol.GetAtoms():
                if atom.GetSymbol() == "*":
                    atom.SetProp("atomLabel", "CoA")
                    reactant_svg = mol_to_svg_image(reactant_mol)
            st.image(reactant_svg, use_container_width=True)

            st.markdown("**Product**")
            
            # Adding CoA group to the product SMILES

            product_smiles = step['output']
            product_smiles = product_smiles + "*"
            product_mol = Chem.MolFromSmiles(product_smiles)
            for atom in product_mol.GetAtoms():
                if atom.GetSymbol() == "*":
                    atom.SetProp("atomLabel", "CoA")
                    product_svg = mol_to_svg_image(product_mol)
            st.image(product_svg, use_container_width=True)

        with col2:
            st.markdown(f"**{step['step']}**")
            
            if step['step'] == "Dehydrogenation":
                st.markdown("- **Action:** Introduces a trans double bond between α and β carbons (Δ²-enoyl-CoA) via FAD-dependent oxidation")
                st.latex(r"Acyl-CoA + FAD \rightarrow Unsaturatedacyl-CoA + FADH_2")
                expected_atp = 1.5
            elif step['step'] == "Isomerase‑assisted dehydrogenation":
                st.markdown("- Occurs when a double bond is at Δ³ (between β and γ carbons)")
                st.markdown("- **Action:** Shifts the double bond to Δ² (between α and β carbons)")
                expected_atp = 0
            elif step['step'] == "Hydration":
                st.markdown("- **Action:** Adds H₂O across the double bond to form β-hydroxyacyl-CoA")
                st.latex(r"Unsaturatedacyl-CoA + H_2O \rightarrow \beta-hydroxyacyl-CoA")
                expected_atp = 0
            elif step['step'] == "Oxidation":
                st.markdown("- **Action:** Oxidizes β-hydroxy group to a keto group forming β-oxoacyl-CoA via NAD⁺-dependent oxidation")
                st.latex(r"\beta-hydroxyacyl-CoA + NAD^+ \rightarrow \beta-oxoacyl-CoA + NADH + H^+")
                expected_atp = 2.5
            elif step['step'] == "Thiolysis":
                st.markdown("- **Action:** Cleaves β-oxoacyl-CoA using CoA-SH to release acetyl-CoA and a shortened acyl-CoA")
                st.latex(r"\beta-oxoacyl-CoA + HS-CoA \rightarrow Acetyl-CoA + Acyl-CoA(-2C)")
                expected_atp = 10
            
            st.markdown(f"**SMARTS Pattern:** `{step['smarts']}`")
            
            # Show the expected ATP yield for this step
            if step['step'] in ["Dehydrogenation", "Oxidation", "Thiolysis"]:
                st.markdown(f"**ATP Yield:** +{expected_atp:.1f}")
            
            # Show the reaction transformation
            st.markdown("**Reaction:**")
            st.code(f"{step['input']} → {step['output']}")
        
        st.divider()


# ------------------ Streamlit Page ------------------

# Get data from previous pages
fa_data = st.session_state.get('fa_data', {})
smiles = fa_data.get('smiles', None)
if not smiles:
    st.error("No FA data found, please go back to Home.")
    st.stop()
    
# Check if fatty acid has changed or processor doesn't exist
if ('processor' not in st.session_state) or ('previous_fa' not in st.session_state) or (st.session_state.previous_fa != fa_data):
    st.session_state.previous_fa = fa_data
    st.session_state.processor = EnhancedFattyAcidMetabolism(smiles)
    
    with st.spinner("Running beta oxidation..."):
        try:
            # Clear any previous internal state (optional if constructor does this)
            st.session_state.processor.cycle_log = []
            st.session_state.processor.reaction_steps = []
            st.session_state.processor.reaction_results = []
            st.session_state.processor.reaction_descriptions = []

            st.session_state.oxidation_results = st.session_state.processor.run_complete_oxidation()
            
            # Ensure cycle_log exists
            if not hasattr(st.session_state.processor, 'cycle_log'):
                st.session_state.processor.cycle_log = []

        except Exception as e:
            st.error(f"Error running oxidation: {str(e)}")
            st.stop()

# Now safely assign local processor variable from session state
processor = st.session_state.processor
results = st.session_state.oxidation_results

# Display ATP calculations
st.subheader("ATP Yield Calculation")
try:
    atp_data = processor.calculate_atp_yield()
    
    # Get basic values
    activation_cost = atp_data.get('activation_cost', 2)
    base_fadh2_atp = atp_data.get('FADH2_ATP', 0)
    base_nadh_atp = atp_data.get('NADH_ATP', 0)
    base_acetyl_coa_atp = atp_data.get('acetyl_CoA_ATP', 0)
    
    # Calculate expected values based on chain length
    chain_length = processor.chain_length
    is_odd_chain = (chain_length % 2) != 0
    expected_cycles = (chain_length // 2) - 1
    
    # Adjust for double bonds (each bypasses one FADH2 production)
    double_bonds = processor.double_bonds
    adjusted_fadh2 = max(0, (expected_cycles - double_bonds) * 1.5)
    adjusted_nadh = expected_cycles * 2.5
    
    # Calculate expected acetyl-CoA count
    if is_odd_chain:
        expected_acetyl_coa = expected_cycles  # Last cycle produces propionyl-CoA
    else:
        expected_acetyl_coa = expected_cycles + 1  # Plus the final acetyl-CoA
    
    # Manual corrections if automatic calculation is off
    if abs(base_fadh2_atp - adjusted_fadh2) > 0.1:
        fadh2_atp = adjusted_fadh2
    else:
        fadh2_atp = base_fadh2_atp
    
    if abs(base_nadh_atp - adjusted_nadh) > 0.1:
        nadh_atp = adjusted_nadh
    else:
        nadh_atp = base_nadh_atp
    
    # Acetyl-CoA verification
    actual_acetyl_coa = results.get('final_products', {}).get('acetyl_coa_count', 0)
    if actual_acetyl_coa != expected_acetyl_coa:
        acetyl_coa_atp = expected_acetyl_coa * 10
    else:
        acetyl_coa_atp = base_acetyl_coa_atp

    # Calculate total ATP with proper accounting
    propionyl_coa = results.get('final_products', {}).get('propionyl_coa', False) and is_odd_chain
    
    if propionyl_coa:
        propionyl_contribution = 15
        total_atp = fadh2_atp + nadh_atp + acetyl_coa_atp + propionyl_contribution - activation_cost
    else:
        # Even-chain - add final acetyl-CoA if not already counted
        if actual_acetyl_coa == expected_cycles:  # If missing final acetyl-CoA
            acetyl_coa_atp += 10
        total_atp = fadh2_atp + nadh_atp + acetyl_coa_atp - activation_cost

    # Display metrics
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total ATP Yield", f"{total_atp:.1f}")
    col2.metric("Activation Cost", f"{activation_cost:.1f}", delta="-2 ATP")
    col3.metric("FADH₂ Yield", f"{fadh2_atp:.1f}", 
               delta=f"{expected_cycles - double_bonds} × 1.5 ATP")
    col4.metric("NADH Yield", f"{nadh_atp:.1f}", 
               delta=f"{expected_cycles} × 2.5 ATP")

    # Detailed breakdown
    breakdown_text = f"""
    **ATP Calculation Breakdown:**  
    • FADH₂: {fadh2_atp:.1f} ATP ({expected_cycles - double_bonds} × 1.5)  
    • NADH: {nadh_atp:.1f} ATP ({expected_cycles} × 2.5)  
    • Acetyl-CoA: {acetyl_coa_atp:.1f} ATP ({expected_acetyl_coa} × 10)  
    """
    
    if propionyl_coa:
        breakdown_text += f"• Propionyl-CoA: +15.0 ATP (odd-chain metabolism)  \n"
    elif is_odd_chain:
        breakdown_text += f"• Propionyl-CoA: Included in acetyl-CoA count  \n"
    
    breakdown_text += f"• Activation cost: -{activation_cost:.1f} ATP  \n\n"
    breakdown_text += f"**Total: {total_atp:.1f} ATP**"
    
    st.markdown(breakdown_text)
    
    # Show warnings if any adjustments were made
    if double_bonds > 0:
        st.info(f"Note: {double_bonds} double bond(s) detected. Each double bond reduces FADH₂ yield by 1.5 ATP.")

except Exception as e:
    st.error(f"Error displaying ATP calculations: {str(e)}")

# Add cycle selection
cycle_options = []
if hasattr(processor, 'cycle_log') and processor.cycle_log:
    cycle_options = [f"Cycle {i+1}" for i in range(len(processor.cycle_log))]
    cycle_options.append("All Cycles")
    selected_cycle = st.selectbox("Select a cycle to view:", cycle_options, index=len(cycle_options)-1)
else:
    st.warning("No oxidation cycles were recorded. Showing basic reaction steps instead.")
    selected_cycle = "Basic Steps"

# Display reaction steps
st.subheader("Beta Oxidation Steps")

try:
    if selected_cycle == "All Cycles" and hasattr(processor, 'cycle_log') and processor.cycle_log:
        # Display all cycles
        for cycle_num, cycle in enumerate(processor.cycle_log, 1):
            st.markdown(f"### Cycle {cycle_num}")
            _display_cycle_steps(cycle)
            
    elif selected_cycle != "Basic Steps" and hasattr(processor, 'cycle_log') and processor.cycle_log:
        # Display single cycle
        cycle_num = int(selected_cycle.split()[1]) - 1
        cycle = processor.cycle_log[cycle_num]
        st.markdown(f"### Cycle {cycle_num + 1}")
        _display_cycle_steps(cycle)
        
    else:
        # Fallback to showing basic reaction steps if no cycle_log is available
        if hasattr(processor, 'reaction_steps') and processor.reaction_steps:
            for i, (step, desc, mol) in enumerate(zip(
                processor.reaction_steps,
                processor.reaction_descriptions,
                processor.reaction_results
            )):
                col1, col2 = st.columns([1, 3])
                
                with col1:
                    st.markdown("**Structure**")
                    st.image(Draw.MolToImage(mol), use_container_width=True)
                    
                with col2:
                    st.markdown(f"**Step {i+1}: {step}**")
                    st.markdown(desc)
                    st.code(f"SMILES: {Chem.MolToSmiles(mol)}")
                
                st.divider()
        else:
            st.warning("No reaction step data available")

except Exception as e:
    st.error(f"Error displaying oxidation steps: {str(e)}")

# Final products summary
st.subheader("Final Products")
col1, col2 = st.columns(2)
with col1:
    st.markdown(f"**Acetyl-CoA produced:** {results.get('final_products', {}).get('acetyl_coa_count', 0)}")
with col2:
    if results.get('final_products', {}).get('propionyl_coa', False):
        st.markdown("**Propionyl-CoA produced** (odd-chain fatty acid)")
    else:
        st.markdown("**No Propionyl-CoA produced** (even-chain fatty acid)")

# Navigation buttons
col1, col2 = st.columns(2)
with col1:
    if st.button("Back to Lynen Spiral"):
        st.switch_page("pages/Lynen_Spiral.py")
with col2:
    if st.button("Back to Home Page"):
        st.switch_page("Home.py")