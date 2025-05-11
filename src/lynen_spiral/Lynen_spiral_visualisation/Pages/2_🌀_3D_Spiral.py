import streamlit as st
import numpy as np
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import Draw
from enhanced_fatty_acid import EnhancedFattyAcidMetabolism

# --- Page Setup ---
st.set_page_config(layout="wide")
st.title("🌀 Lynen Spiral with Reaction Coloring")

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

# --- Color Scheme for Reaction Types ---
reaction_colors = {
    "Activation": "#FF6B6B",      # Red
    "Dehydrogenation": "#A74ECD", # Teal
    "Hydration": "#45B7D1",       # Blue
    "Oxidation": "#FFA07A",       # Light orange
    "Thiolysis": "#98D8C8"        # Mint green
}

# --- 3D Spiral Visualization ---
def create_lynens_spiral(cycles):
    fig = go.Figure()
    
    # Spiral parameters
    cycles_per_revolution = 1.25
    height_per_cycle = 1.0
    base_radius = 3.0
    
    # Create coordinates with proper step spacing
    x, y, z, colors, labels = [], [], [], [], []
    
    for cycle_idx, cycle in enumerate(cycles):
        steps = cycle["steps"]
        for step_idx, step in enumerate(steps):
            # Calculate position with even spacing
            angle = 2 * np.pi * (cycle_idx/cycles_per_revolution + step_idx/(len(steps)*cycles_per_revolution))
            radius = base_radius * (0.93 ** cycle_idx)  # Slight radius decrease
            
            x.append(radius * np.cos(angle))
            y.append(radius * np.sin(angle))
            z.append(cycle_idx * height_per_cycle)
            
            # Get reaction type (first word of step name)
            reaction_type = step["step"].split()[0]
            colors.append(reaction_colors.get(reaction_type, "#AAAAAA"))
            labels.append(f"Cycle {cycle_idx+1}: {step['step']}")

    # Create spiral path
    fig.add_trace(go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode='markers+lines',
        marker=dict(
            size=10,
            color=colors,
            opacity=0.9,
            line=dict(width=1, color='white')
        ),
        line=dict(
            color='rgba(1,1,1,0.3)',
            width=4
        ),
        text=labels,
        hoverinfo='text',
        name='Pathway'
    ))
    
    # Add cycle labels
    for cycle_idx in range(total_cycles):
        angle = 2 * np.pi * (cycle_idx/cycles_per_revolution)
        radius = base_radius * (0.93 ** cycle_idx) * 1.15
        fig.add_trace(go.Scatter3d(
            x=[radius * np.cos(angle)],
            y=[radius * np.sin(angle)],
            z=[cycle_idx * height_per_cycle],
            mode='text',
            text=[f"<b>Cycle {cycle_idx+1}</b>"],
            textfont=dict(size=14, color='black'),
            showlegend=False
        ))
    
    # Add legend for reaction types
    for reaction, color in reaction_colors.items():
        fig.add_trace(go.Scatter3d(
            x=[None],
            y=[None],
            z=[None],
            mode='markers',
            marker=dict(size=10, color=color),
            name=reaction,
            hoverinfo='none'
        ))
    
    # Clean layout
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False, showbackground=False),
            yaxis=dict(visible=False, showbackground=False),
            zaxis=dict(
                visible=False,
                showbackground=False,
                showgrid=False
            ),
            bgcolor='white',
            camera=dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0.2),
                eye=dict(x=1.5, y=1.5,z=0.8)
            ),
            aspectratio=dict(x=1, y=1, z=0.8)
        ),
        title=dict(
            text=f"<b>Lynen Spiral: {fa_data['name']}</b>",
            y=0.95,  # Position title higher
            x=0.5,   # Center the title
            xanchor='center',
            yanchor='top',
            font=dict(size=18)
        ),
        margin=dict(l=0, r=0, b=0, t=80),
        paper_bgcolor='white',
        height=800,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=0.8,
            xanchor="right",
            x=1
        )
    )
    return fig

# --- Main Display ---
col1, col2 = st.columns([3, 1])

with col1:
    st.markdown("<div style='height: 20px;'></div>", unsafe_allow_html=True)  # Add some space
    st.plotly_chart(create_lynens_spiral(cycles), use_container_width=True)  # This line was missing
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
        reaction_type = step["step"].split()[0]
        color = reaction_colors.get(reaction_type, "#AAAAAA")
        st.markdown(
            f"- <span style='color: {color}; font-weight: bold;'>{step['step']}</span>",
            unsafe_allow_html=True
        )
    
    st.markdown("---")
    st.subheader("Molecular View")
    step_to_show = st.select_slider(
        "Select step",
        options=list(range(len(cycle_data["steps"]))),
        format_func=lambda x: f"Step {x+1}: {cycle_data['steps'][x]['step']}"
    )
    
    step_smiles = cycle_data["steps"][step_to_show]["output"]
    # Adding the CoA group to the SMILES string for visualization 
    # (Workaround to visialize CoA group in the molecule without getting an error from RDKit))
    step_smiles = step_smiles + "*"
    step_mol = Chem.MolFromSmiles(step_smiles)
    for atom in step_mol.GetAtoms():
        if atom.GetSymbol() == "*":
            atom.SetProp("atomLabel", "CoA")
    mol_img = Draw.MolToImage(step_mol, size=(300, 300))
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