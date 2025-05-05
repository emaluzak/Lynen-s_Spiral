import streamlit as st
import numpy as np
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import sys

# Add path to your fatty_acid module
sys.path.append(r'C:\Users\ipeki\git\Lynen-s_Spiral\Lynen-s_Spiral\src\lynen_spiral')
from enhanced_fatty_acid import EnhancedFattyAcidMetabolism  # type: ignore

# --- Setup ---
st.set_page_config(layout="wide", page_title="Lynen Spiral β-Oxidation Visualizer")
st.title("🌀 Lynen Spiral: β-Oxidation of Fatty Acids")

# --- Sidebar: Input ---
fatty_acid_choice = st.sidebar.selectbox(
    "Choose a fatty acid or enter SMILES:",
    ["palmitic acid", "stearic acid", "oleic acid", "linoleic acid", "alpha-linolenic acid", "arachidonic acid", "Custom"]
)

if fatty_acid_choice == "Custom":
    input_smiles = st.sidebar.text_input("Enter SMILES string:", value="CCCCCCCCCCCCCCCC(=O)O")
else:
    input_smiles = fatty_acid_choice

# --- β-Oxidation Simulation ---
with st.spinner("Running β-oxidation..."):
    em = EnhancedFattyAcidMetabolism(input_smiles)
    em.run_complete_oxidation()
    viz_data = em.prepare_data_for_visualization()
    steps = viz_data["steps"]
    total_atp = viz_data["metadata"]["total_atp_yield"]
    cycles = viz_data["cycles"]

# --- Sidebar Navigation ---
st.sidebar.markdown(f"**Total ATP Yield:** {total_atp:.1f}")
st.sidebar.markdown(f"**Total Cycles:** {len(cycles)}")
st.sidebar.markdown("---")

# --- Select View ---
view_mode = st.sidebar.radio("View Mode", ["3D Spiral of Cycles", "2D View of Selected Cycle"])

# --- Render RDKit Molecule ---
def mol_to_image(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))
        buf = BytesIO()
        img.save(buf, format='PNG')
        return buf.getvalue()
    return None

# --- 3D Spiral View ---
def make_3d_spiral(cycles):
    fig = go.Figure()
    steps_per_cycle = len(cycles[0]["steps"])
    r = 1.0
    angle_step = 2 * np.pi / steps_per_cycle
    z_step_per_cycle = 1.0

    x_vals, y_vals, z_vals, names, cycle_nums = [], [], [], [], []

    for i, cycle in enumerate(cycles):
        for j, step in enumerate(cycle["steps"]):
            angle = j * angle_step
            x = r * np.cos(angle)
            y = r * np.sin(angle)
            z = i * z_step_per_cycle
            x_vals.append(x)
            y_vals.append(y)
            z_vals.append(z)
            names.append(step["step"])
            cycle_nums.append(f"Cycle {i + 1}")

    fig.add_trace(go.Scatter3d(
        x=x_vals,
        y=y_vals,
        z=z_vals,
        mode='lines+markers+text',
        text=names,
        textposition='top center',
        marker=dict(size=4, color=list(range(len(x_vals))), colorscale='Viridis', showscale=True),
        line=dict(width=3),
        name="Lynen Spiral"
    ))

    # Add cycle labels as annotations
    for i in range(len(cycles)):
        x = r * np.cos(0)
        y = r * np.sin(0)
        z = i * z_step_per_cycle
        fig.add_trace(go.Scatter3d(
            x=[x], y=[y], z=[z + 0.5],
            mode="text",
            text=[f"Cycle {i + 1}"],
            textposition="top center",
            showlegend=False
        ))

    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(title="Cycle Progression"),
        ),
        margin=dict(l=0, r=0, t=40, b=0),
        title="3D Spiral of β-Oxidation Cycles"
    )
    return fig

# --- 2D Cycle View ---
def make_2d_circle_for_cycle(cycle, cycle_id):
    steps_in_cycle = cycle["steps"]
    n = len(steps_in_cycle)
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x = np.cos(theta)
    y = np.sin(theta)
    labels = [step["step"] for step in steps_in_cycle]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=np.append(x, x[0]),
        y=np.append(y, y[0]),
        mode='lines+markers+text',
        text=labels,
        textposition='top center',
        marker=dict(size=10),
        line=dict(color='blue')
    ))
    fig.update_layout(
        xaxis=dict(visible=False), yaxis=dict(visible=False),
        width=600, height=600,
        margin=dict(l=0, r=0, t=30, b=0),
        title=f"Cycle {cycle_id} Steps"
    )
    return fig

# --- Render View ---
if view_mode == "3D Spiral of Cycles":
    st.plotly_chart(make_3d_spiral(cycles), use_container_width=True)
else:
    cycle_id = st.sidebar.selectbox("Select Cycle", range(1, len(cycles)+1))
    selected_cycle = cycles[cycle_id - 1]
    step_index = st.sidebar.slider("Step in Cycle", 1, len(selected_cycle["steps"]), 1) - 1
    step = selected_cycle["steps"][step_index]
    st.plotly_chart(make_2d_circle_for_cycle(selected_cycle, cycle_id), use_container_width=True)

    img_data = mol_to_image(step["output"])
    if img_data:
        st.image(img_data, caption=f"Step {step_index + 1}: {step['step']}")

    st.subheader("Step Details")
    st.markdown(f"**Step {step_index + 1}: {step['step']}**")
    st.markdown(f"- **Input SMILES:** `{step['input']}`")
    st.markdown(f"- **Output SMILES:** `{step['output']}`")
    st.markdown(f"- **SMARTS:** `{step['smarts']}`")
    st.markdown(f"- **ATP Yield:** `{step['atp_yield']}`")
