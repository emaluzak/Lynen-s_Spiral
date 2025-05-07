import streamlit as st
import numpy as np
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from enhanced_fatty_acid import EnhancedFattyAcidMetabolism

# --- Page Setup ---
st.set_page_config(layout="wide")
st.title("🌀 3D Lynen Spiral Visualization")

# --- Check Session State ---
if 'fa_data' not in st.session_state:
    st.error("❌ No molecule data found. Please start from the Home page.")
    if st.button("← Return to Home"):
        st.switch_page("app.py")
    st.stop()

# --- Load Data ---
fa_data = st.session_state['fa_data']
mol = fa_data['molecule']

# --- Run β-Oxidation ---
@st.cache_data
def run_beta_oxidation(_mol):
    metabolism = EnhancedFattyAcidMetabolism(Chem.MolToSmiles(_mol))
    metabolism.run_complete_oxidation()
    return metabolism.prepare_data_for_visualization()

with st.spinner("Simulating β-oxidation cycles..."):
    viz_data = run_beta_oxidation(mol)
    cycles = viz_data["cycles"]
    total_cycles = len(cycles)
    total_atp = viz_data["metadata"]["total_atp_yield"]

# --- 3D Spiral Visualization ---
def create_lynens_spiral(cycles):
    fig = go.Figure()
    
    # Spiral parameters
    cycles_per_revolution = 1.25  # Controls spiral tightness
    height_per_cycle = 0.5        # Vertical spacing
    base_radius = 2.0             # Starting radius
    
    # Color scheme for different steps
    step_colors = {
        "Activation": "#FF6B6B",
        "Dehydrogenation": "#4ECDC4",
        "Hydration": "#45B7D1",
        "Oxidation": "#FFA07A",
        "Thiolysis": "#98D8C8"
    }
    
    # Create coordinates
    x, y, z, labels, colors = [], [], [], [], []
    
    for cycle_idx, cycle in enumerate(cycles):
        steps = cycle["steps"]
        for step_idx, step in enumerate(steps):
            # Spiral coordinates
            angle = 2 * np.pi * (cycle_idx/cycles_per_revolution + step_idx/(len(steps)*cycles_per_revolution))
            radius = base_radius * (0.9 ** cycle_idx)  # Slightly decrease radius each cycle
            
            x.append(radius * np.cos(angle))
            y.append(radius * np.sin(angle))
            z.append(cycle_idx * height_per_cycle)
            
            # Labels and colors
            step_name = step["step"]
            labels.append(f"<b>Cycle {cycle_idx+1}</b><br>{step_name}<br>ATP: {step['atp_yield']}")
            colors.append(step_colors.get(step_name.split()[0], "#AAAAAA"))  # Default gray
    
    # Create spiral trace
    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers+lines+text',
        marker=dict(
            size=8,
            color=colors,
            opacity=0.9,
            line=dict(width=1, color='black')
        ),
        line=dict(
            color='rgba(150,150,150,0.4)',
            width=4
        ),
        text=labels,
        hoverinfo='text',
        name='β-Oxidation Pathway'
    ))
    
    # Add cycle labels
    for cycle_idx in range(total_cycles):
        angle = 2 * np.pi * (cycle_idx/cycles_per_revolution)
        radius = base_radius * (0.9 ** cycle_idx)
        fig.add_trace(go.Scatter3d(
            x=[radius * np.cos(angle)],
            y=[radius * np.sin(angle)],
            z=[cycle_idx * height_per_cycle],
            mode='text',
            text=[f"<b>Cycle {cycle_idx+1}</b>"],
            textfont=dict(size=14),
            showlegend=False
        ))
    
    # Layout configuration
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(title='<b>Cycle Progression</b>'),
            camera=dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=1.5, y=1.5, z=0.6)
            ),
            aspectratio=dict(x=1, y=1, z=0.7)
        ),
        title=f"<b>Lynen Spiral: {fa_data['name']}</b><br>Total ATP Yield: {total_atp}",
        margin=dict(l=0, r=0, b=0, t=100),
        height=800,
        showlegend=False
    )
    return fig
# --- Main Display ---
col1, col2 = st.columns([3, 1])

with col1:
    st.plotly_chart(create_lynens_spiral(cycles), use_container_width=True)

with col2:
    st.subheader("Cycle Navigator")
    selected_cycle = st.selectbox(
        "Select cycle",
        range(1, total_cycles+1),
        format_func=lambda x: f"Cycle {x}"
    )
    
    cycle_data = cycles[selected_cycle-1]
    st.markdown(f"**Steps in Cycle {selected_cycle}:**")
    for step in cycle_data["steps"]:
        st.markdown(f"- {step['step']} (ATP: {step['atp_yield']})")
    
    st.markdown("---")
    st.subheader("Molecular View")
    step_to_show = st.select_slider(
        "Select step",
        options=list(range(len(cycle_data["steps"]))),
        format_func=lambda x: f"Step {x+1}"
    )
    
    step_smiles = cycle_data["steps"][step_to_show]["output"]
    mol_img = Draw.MolToImage(Chem.MolFromSmiles(step_smiles), size=(300, 300))
    st.image(mol_img, caption=f"Step {step_to_show+1} Structure")

# --- Navigation ---
st.markdown("---")
col1, col2 = st.columns(2)
with col1:
    if st.button("← Back to Home", use_container_width=True):
        st.switch_page("app.py")
with col2:
    if st.button("View 2D Cycle Details →", use_container_width=True):
        st.switch_page("pages/3_📊_2D_Cycle.py")