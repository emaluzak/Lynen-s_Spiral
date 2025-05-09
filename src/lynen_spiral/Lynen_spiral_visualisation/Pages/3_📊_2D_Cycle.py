import streamlit as st 
from rdkit import Chem
from rdkit.Chem import Draw
import sys
import os

st.set_page_config(layout="wide")
st.title("2D Cycle Representation")

sys.path.append(os.path.abspath(r"C:\Users\Lenovo\git\Lynen-s_Spiral\src\lynen_spiral\Lynen_spiral_visualisation"))
from enhanced_fatty_acid import EnhancedFattyAcidMetabolism

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
        st.session_state.oxidation_results = processor.run_complete_oxidation()

results = st.session_state.oxidation_results

# Display ATP calculations
st.subheader("ATP Yield Calculation")
atp_data = processor.calculate_atp_yield()

col1, col2, col3, col4 = st.columns(4)
col1.metric("Total ATP Yield", f"{atp_data['total_ATP']:.1f}")
col2.metric("Activation Cost", f"{atp_data['breakdown']['activation_cost']:.1f}", delta="-2 ATP")
col3.metric("FADH2 Yield", f"{atp_data['breakdown']['fadh2_yield']:.1f}", delta="1.5 ATP each")
col4.metric("NADH Yield", f"{atp_data['breakdown']['nadh_yield']:.1f}", delta="2.5 ATP each")

st.markdown("""
**ATP Calculation Breakdown:**
- Each FADH₂ produces ~1.5 ATP (through electron transport chain)
- Each NADH produces ~2.5 ATP (through electron transport chain)
- Each acetyl-CoA produces ~10 ATP (through TCA cycle)
- Activation costs 2 ATP (ATP → AMP + 2Pi)
""")

# Add cycle selection at the top
if processor.cycle_log:
    cycle_options = [f"Cycle {i+1}" for i in range(len(processor.cycle_log))]
    cycle_options.append("All Cycles")
    selected_cycle = st.selectbox("Select a cycle to view:", cycle_options, index=len(cycle_options)-1)
else:
    st.warning("No oxidation cycles found")
    selected_cycle = "All Cycles"

# Display selected cycle(s)
st.subheader("Beta Oxidation Steps")

if selected_cycle == "All Cycles":
    cycles_to_display = processor.cycle_log
else:
    cycle_num = int(selected_cycle.split()[1]) - 1
    cycles_to_display = [processor.cycle_log[cycle_num]]

for cycle_num, cycle in enumerate(cycles_to_display, 1 if selected_cycle == "All Cycles" else int(selected_cycle.split()[1])):
    st.markdown(f"### Cycle {cycle_num}")
    
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
            elif step['step'] == "Hydration":
                st.markdown("- H₂O added across double bond")
                st.markdown("- Catalyzed by Enoyl-CoA hydratase")
            elif step['step'] == "Oxidation":
                st.markdown("- NAD⁺ → NADH (2.5 ATP)")
                st.markdown("- Catalyzed by 3-Hydroxyacyl-CoA dehydrogenase")
            elif step['step'] == "Thiolysis":
                st.markdown("- Cleavage by thiolase")
                st.markdown("- Produces Acetyl-CoA")
            
            st.markdown(f"**SMARTS Pattern:** `{step['smarts']}`")
            
            if step['atp_yield'] > 0:
                st.markdown(f"**ATP Yield:** +{step['atp_yield']:.1f}")
            
            # Show the reaction transformation
            st.markdown("**Reaction:**")
            st.code(f"{step['input']} → {step['output']}")
        
        st.divider()

# Final products summary
st.subheader("Final Products")
col1, col2 = st.columns(2)
with col1:
    st.markdown(f"**Acetyl-CoA produced:** {results['final_products']['acetyl_coa_count']}")
with col2:
    if results['final_products']['propionyl_coa']:
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
