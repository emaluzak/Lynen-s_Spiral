import streamlit as st
import numpy as np
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import sys
from io import BytesIO
from enhanced_fatty_acid import EnhancedFattyAcidMetabolism


# --- Setup ---
st.set_page_config(layout="wide", page_title="Lynen Spiral Visualizer")
st.title("🌀 3D Lynen Spiral: β-Oxidation Pathway")

# --- Load Fatty Acid Data ---
if 'fa_data' not in st.session_state:
    st.warning("Please start from the home page")
    st.stop()

fa_data = st.session_state['fa_data']
mol = fa_data['molecule']

# --- Run β-Oxidation ---
@st.cache_data
def run_beta_oxidation(mol):
    metabolism = EnhancedFattyAcidMetabolism(Chem.MolToSmiles(mol))
    metabolism.run_complete_oxidation()
    return metabolism.prepare_data_for_visualization()

viz_data = run_beta_oxidation(mol)
cycles = viz_data["cycles"]

# --- 3D Spiral Visualization ---
def create_lynens_spiral(cycles):
    fig = go.Figure()
    
    # Spiral parameters
    cycles_per_revolution = 1.5  # How many cycles per full spiral turn
    height_per_cycle = 0.5       # Vertical spacing between cycles
    radius = 2.0                 # Base radius of spiral
    
    # Create spiral coordinates
    all_x, all_y, all_z = [], [], []
    step_labels = []
    colors = []
    cycle_labels = []
    
    for cycle_idx, cycle in enumerate(cycles):
        steps = cycle["steps"]
        n_steps = len(steps)
        
        for step_idx, step in enumerate(steps):
            # Spiral coordinates
            angle = 2 * np.pi * (cycle_idx / cycles_per_revolution + step_idx / (n_steps * cycles_per_revolution))
            r = radius * (1 - 0.1 * cycle_idx/len(cycles))  # Slightly decrease radius
            
            x = r * np.cos(angle)
            y = r * np.sin(angle)
            z = cycle_idx * height_per_cycle
            
            all_x.append(x)
            all_y.append(y)
            all_z.append(z)
            step_labels.append(f"Cycle {cycle_idx+1}, Step {step_idx+1}: {step['step']}")
            colors.append(step_idx)  # Color by step position
            cycle_labels.append(f"Cycle {cycle_idx+1}")
    
    # Create the spiral trace
    fig.add_trace(go.Scatter3d(
        x=all_x,
        y=all_y,
        z=all_z,
        mode='markers+lines+text',
        marker=dict(
            size=6,
            color=colors,
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title='Step in Cycle')
        ),
        line=dict(
            width=4,
            color='rgba(100,100,100,0.4)'
        ),
        text=step_labels,
        hoverinfo='text',
        name='β-Oxidation Pathway'
    ))
    
    # Add cycle labels
    for cycle_idx in range(len(cycles)):
        angle = 2 * np.pi * (cycle_idx / cycles_per_revolution)
        r = radius * (1 - 0.1 * cycle_idx/len(cycles))
        x = r * np.cos(angle)
        y = r * np.sin(angle)
        z = cycle_idx * height_per_cycle
        
        fig.add_trace(go.Scatter3d(
            x=[x],
            y=[y],
            z=[z],
            mode='text',
            text=[f"<b>Cycle {cycle_idx+1}</b>"],
            textfont=dict(size=14, color='black'),
            hoverinfo='none',
            showlegend=False
        ))
    
    # Layout configuration
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(title='Cycle Progression'),
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=0.8)
            ),
            aspectratio=dict(x=1, y=1, z=0.7)
        ),
        title=f"3D Lynen Spiral: {fa_data['name']} (β-Oxidation)",
        margin=dict(l=0, r=0, b=0, t=40),
        height=800,
        showlegend=False
    )
    
    return fig

# --- Display the Visualization ---
st.plotly_chart(create_lynens_spiral(cycles), use_container_width=True)

# --- Cycle Information ---
st.subheader("Cycle Information")
selected_cycle = st.selectbox("Select cycle to view details", 
                             range(1, len(cycles)+1), 
                             format_func=lambda x: f"Cycle {x}")

cycle_data = cycles[selected_cycle-1]
col1, col2 = st.columns(2)

with col1:
    st.markdown("### Steps in Cycle")
    for step_idx, step in enumerate(cycle_data["steps"]):
        st.markdown(f"**{step_idx+1}. {step['step']}**")
        st.code(f"SMARTS: {step['smarts']}\nATP: {step['atp_yield']}")

with col2:
    st.markdown("### Molecular View")
    step_to_show = st.select_slider("Select step to view molecule", 
                                   options=list(range(len(cycle_data["steps"]))), 
                                   format_func=lambda x: f"Step {x+1}")
    
    step_smiles = cycle_data["steps"][step_to_show]["output"]
    mol = Chem.MolFromSmiles(step_smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))
        st.image(img, caption=f"Step {step_to_show+1} Structure")

# --- Navigation ---
st.markdown("---")
if st.button("← Back to Home"):
    st.switch_page("app.py")