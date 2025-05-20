import streamlit as st
import io
import sys
import os
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image

st.set_page_config(layout="wide")
st.title("2D Cycle Representation")

sys.path.append(os.path.abspath(r"C:\Users\Lenovo\git\Lynen-s_Spiral\src\lynen_spiral\Lynen_spiral_visualisation"))
from enhanced_fatty_acid import EnhancedFattyAcidMetabolism

def draw_mol_with_coa(smiles):
    """
    Generates a 2D depiction of a molecule from a SMILES string.
    Replaces the terminal sulfur at the end of the fatty acid by a "S-CoA" label.
    The CoA label was omitted in the original code, as it was unrecognized by RDKit and lead to errors.
    It is now added as a label to the sulfur atom in the 2D rendering to better simulate the beta-oxidation process.

    Arguments:
        smiles (str): The SMILES string representing the molecule.

    Returns:
        PIL.Image.Image: An image of the rendered molecule, or None if the SMILES string is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Generate 2D coordinates
    AllChem.Compute2DCoords(mol)

    # Create the drawer and options
    drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)
    opts = drawer.drawOptions()

    # Set global atom label color to black
    opts.atomLabelColour = (0, 0, 0)  # Black text

    # Customize sulfur color to yellow
    opts.elemDict = {'S': (1.0, 1.0, 0.0)}  # Yellow sulfur

    # Modify the atom labels
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "S" and atom.GetDegree() == 1:
            # Set the label for sulfur to include "S-CoA"
            atom.SetProp("atomLabel", "S-CoA")

    # Draw the molecule
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    # Convert the drawing to an image
    return Image.open(io.BytesIO(drawer.GetDrawingText()))

def _display_cycle_steps(cycle):
    """Helper function to display steps of a single cycle."""
    for step in cycle['steps']:
        col1, col2 = st.columns([1, 3])

        with col1:

            st.markdown("**Reactant**")
            st.image(draw_mol_with_coa(step['input']), use_container_width=True)

            st.markdown("**Product**")
            st.image(draw_mol_with_coa(step['output']), use_container_width=True)

        with col2:
            # Display step information
            st.markdown(f"**{step['step']}**")
            
            if step['step'] == "Dehydrogenation":
                st.markdown("- FAD → FADH₂ (1.5 ATP)")
                st.markdown("- Catalyzed by Acyl-CoA dehydrogenase")
                expected_atp = 1.5
            elif step['step'] == "Hydration":
                st.markdown("- H₂O added across double bond")
                st.markdown("- Catalyzed by Enoyl-CoA hydratase")
                expected_atp = 0
            elif step['step'] == "Oxidation":
                st.markdown("- NAD⁺ → NADH (2.5 ATP)")
                st.markdown("- Catalyzed by 3-Hydroxyacyl-CoA dehydrogenase")
                expected_atp = 2.5
            elif step['step'] == "Thiolysis":
                st.markdown("- Cleavage by thiolase")
                st.markdown("- Produces Acetyl-CoA (10 ATP via TCA cycle)")
                expected_atp = 10
            
            st.markdown(f"**SMARTS Pattern:** `{step['smarts']}`")
            
            # Show the expected ATP yield for this step
            if step['step'] in ["Dehydrogenation", "Oxidation", "Thiolysis"]:
                st.markdown(f"**ATP Yield:** +{expected_atp:.1f}")
            
            # Show the reaction transformation
            st.markdown("**Reaction:**")
            st.code(f"{step['input']} → {step['output']}")
        
        st.divider()

# Get data from previous pages
fa_data = st.session_state.get('fa_data', {})
processor = fa_data.get('processor', None)

if not processor:
    st.warning("No fatty acid data found. Please go back to the home page and start again.")
    if st.button("Back to Home"):
        st.switch_page("app.py")
    st.stop()

# Display the current fatty acid being analyzed
current_fa = fa_data.get('name', 'Unknown Fatty Acid')
st.subheader(f"Analyzing: {current_fa}")

# Clear previous oxidation results if the fatty acid has changed
if 'previous_fa' not in st.session_state:
    st.session_state.previous_fa = current_fa
elif st.session_state.previous_fa != current_fa:
    st.session_state.pop('oxidation_results', None)
    st.session_state.previous_fa = current_fa

# Run complete oxidation if not already done or if fatty acid changed
if 'oxidation_results' not in st.session_state:
    with st.spinner("Running beta oxidation..."):
        try:
            # Reset processor state
            processor.cycle_log = []
            processor.reaction_steps = []
            processor.reaction_results = []
            processor.reaction_descriptions = []
            
            st.session_state.oxidation_results = processor.run_complete_oxidation()
            
            # Ensure cycle_log exists
            if not hasattr(processor, 'cycle_log'):
                processor.cycle_log = []
                
        except Exception as e:
            st.error(f"Error running oxidation: {str(e)}")
            st.stop()

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
    col3.metric("FADH2 Yield", f"{fadh2_atp:.1f}", 
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
    if st.button("Back to 3D Spiral"):
        st.switch_page("pages/3D_Lynen_Spiral.py")
with col2:
    if st.button("Back to Home"):
        st.switch_page("Home.py")