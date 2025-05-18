import streamlit as st
import numpy as np
import plotly.graph_objects as go
import io
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

smiles = fa_data.get('smiles')
if not smiles:
    st.error("No SMILES string found in session. Please restart.")
    st.stop()

mol = Chem.MolFromSmiles(smiles)
if mol is None:
    st.error("Invalid SMILES string stored. Please restart.")
    st.stop()

processor = EnhancedFattyAcidMetabolism(smiles)

# --- Run β-Oxidation ---
@st.cache_data
def run_beta_oxidation(smiles):

    """
    Simulate the complete β-oxidation process for a given fatty acid molecule.

    Parameters:
        smiles (str): The SMILES string representation of the fatty acid molecule.

    Returns:
        dict: A structured data dictionary prepared for visualization, containing
              detailed information about each reaction step, cycles, and overall metabolism.

    Notes:
        - This function uses the EnhancedFattyAcidMetabolism class to run the
          complete oxidation simulation.
        - The output is cached by Streamlit to avoid redundant calculations on repeated calls
          with the same SMILES input.
    """

    metabolism = EnhancedFattyAcidMetabolism(smiles)
    metabolism.run_complete_oxidation()
    return metabolism.prepare_data_for_visualization()

with st.spinner("Simulating β-oxidation cycles..."):
    viz_data = run_beta_oxidation(smiles)
    cycles = viz_data["cycles"]
    total_cycles = len(cycles)

# --- Color Scheme for Reaction Types ---
reaction_colors = {
    "Activation": "#FF6B6B",      # Red
    "Dehydrogenation": "#A74ECD", # Teal
    "Isomerase‑assisted dehydrogenation": "#FFA1E3", # Light pink
    "Hydration": "#45B7D1",       # Blue
    "Oxidation": "#FFA07A",       # Light orange
    "Thiolysis": "#98D8C8"        # Mint green
}

# --- 3D Spiral Visualization ---
def create_lynens_spiral(cycles):
    
    """
    Generate a 3D Lynen Spiral visualization of β-oxidation cycles using Plotly.

    This function maps each reaction step within each β-oxidation cycle onto a 3D
    spiral structure, where:
        - The spiral winds with a slight radius decrease per cycle.
        - Each step is positioned based on its cycle index and step index.
        - Steps are colored according to predefined reaction types.
        - Cycle labels are added at appropriate positions.
        - A legend clarifies the colors associated with different reaction types.
        - The layout is customized for clean visualization with hidden axes and a camera
          viewpoint optimized for viewing the spiral.

    Parameters:
        cycles (list): A list of cycle dictionaries, each containing reaction steps with
                       their names and other relevant data. (retrieved from the
                       EnhancedFattyAcidMetabolism class, prepare_data_for_visualization method)

    Returns:
        plotly.graph_objects.Figure: A 3D scatter plot figure representing the Lynen Spiral,
                                    with reaction steps colored by reaction type, cycle labels,
                                    and a legend for reaction types.

    Notes:
        - Requires global 'reaction_colors' dictionary for mapping reaction types to colors.
        - Uses global 'total_cycles' and 'fa_data' variables for labeling and title text.
    """

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
            reaction_type = step["step"]
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
col1, col2 = st.columns([2, 1])

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
        reaction_type = step["step"]
        color = reaction_colors.get(reaction_type, "#AAAAAA")
        st.markdown(
            f"- <span style='color: {color}; font-weight: bold;'>{step['step']}</span>",
            unsafe_allow_html=True
        )

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
    svg = mol_to_svg_image(step_mol)
    st.image(svg, caption=f"Step {step_to_show+1} Structure", use_container_width=True)


# --- Navigation ---
st.markdown("---")
col1, col2 = st.columns(2)
with col1:
    if st.button("← Back to Home", use_container_width=True):
        st.switch_page("app.py")
with col2:
    if st.button("View 2D Cycle Details →", use_container_width=True):
        st.switch_page("pages/3_📊_Detailed_Analysis.py")