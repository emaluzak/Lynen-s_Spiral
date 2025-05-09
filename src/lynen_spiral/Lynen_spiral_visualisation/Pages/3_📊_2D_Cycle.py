import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
import sys
import os

st.set_page_config(layout="wide")
st.title("2D Cycle Representation")

sys.path.append(os.path.abspath(r"C:\Users\Lenovo\git\Lynen-s_Spiral\src\lynen_spiral\Lynen_spiral_visualisation"))
from enhanced_fatty_acid import EnhancedFattyAcidMetabolism

def _display_cycle_steps(cycle):
    """Helper function to display steps of a single cycle."""
    for step in cycle['steps']:
        col1, col2 = st.columns([1, 3])
        
        with col1:
            # Display reactant and product structures
            reactant = Chem.MolFromSmiles(step['input'])
            product = Chem.MolFromSmiles(step['output'])
            
            st.markdown("**Reactant**")
            st.image(Draw.MolToImage(reactant), use_container_width=True)
            
            st.markdown("**Product**")
            st.image(Draw.MolToImage(product), use_container_width=True)
            
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

# Run complete oxidation if not already done
if 'oxidation_results' not in st.session_state:
    with st.spinner("Running beta oxidation..."):
        try:
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
    fadh2_atp = atp_data.get('FADH2_ATP', 0)
    nadh_atp = atp_data.get('NADH_ATP', 0)
    acetyl_coa_atp = atp_data.get('acetyl_CoA_ATP', 0)
    
    # Calculate total ATP with proper accounting for final products
    if results.get('final_products', {}).get('propionyl_coa', False):
        # Odd-chain fatty acid - add propionyl-CoA contribution (15 ATP)
        propionyl_contribution = 15
        total_atp = fadh2_atp + nadh_atp + acetyl_coa_atp + propionyl_contribution - activation_cost
        acetyl_contribution = 0
    else:
        # Even-chain fatty acid - add final acetyl-CoA contribution (10 ATP)
        acetyl_contribution = 10
        total_atp = fadh2_atp + nadh_atp + acetyl_coa_atp + acetyl_contribution - activation_cost
        propionyl_contribution = 0

    # Display
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total ATP Yield", f"{total_atp:.1f}")
    col2.metric("Activation Cost", f"{activation_cost:.1f}", delta="-2 ATP")
    col3.metric("FADH2 Yield", f"{fadh2_atp:.1f}", delta="1.5 ATP each")
    col4.metric("NADH Yield", f"{nadh_atp:.1f}", delta="2.5 ATP each")

    # Detailed breakdown with proper bullet points
    breakdown_text = f"""
    **ATP Calculation Breakdown:**
    • FADH₂: {fadh2_atp:.1f} ATP (1.5 per FADH₂)  
    • NADH: {nadh_atp:.1f} ATP (2.5 per NADH)  
    • Acetyl-CoA: {acetyl_coa_atp + (acetyl_contribution if not propionyl_contribution else 0):.1f} ATP (10 per acetyl-CoA)  
    """
    
    if propionyl_contribution > 0:
        breakdown_text += f"• Propionyl-CoA: +{propionyl_contribution:.1f} ATP (special metabolism)  \n"
    
    breakdown_text += f"• Activation cost: -{activation_cost:.1f} ATP  \n\n"
    breakdown_text += f"**Total: {total_atp:.1f} ATP**"
    
    st.markdown(breakdown_text)

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
        st.switch_page("pages/2_🌀_3D_Spiral.py")
with col2:
    if st.button("Back to Home"):
        st.switch_page("app.py")