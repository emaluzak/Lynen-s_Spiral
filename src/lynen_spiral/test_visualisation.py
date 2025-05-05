import streamlit as st
import plotly.graph_objects as go
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import pandas as pd
import io
import base64

# Set page config
st.set_page_config(layout="wide", page_title="β-Oxidation Visualizer")

# Title and introduction
st.title("Interactive β-Oxidation Visualizer")
st.markdown("""
This tool visualizes fatty acid breakdown through β-oxidation, showing both the Lynen Spiral 
and detailed reaction cycle. Navigate through each step to see molecular transformations.
""")

# Sidebar for user input
with st.sidebar:
    st.header("Fatty Acid Input")
    fatty_acid = st.selectbox("Select fatty acid", 
                            ["Palmitic acid (16:0)", "Stearic acid (18:0)", 
                             "Oleic acid (18:1 Δ⁹)", "Custom"])
    
    if fatty_acid == "Custom":
        custom_length = st.number_input("Carbon atoms", min_value=4, max_value=30, value=16)
        custom_unsat = st.number_input("Double bonds", min_value=0, max_value=10, value=0)
        custom_pos = st.text_input("Double bond positions (comma separated)", "9") if custom_unsat > 0 else ""
    
    st.markdown("---")
    st.markdown("**Visualization Options**")
    show_mechanism = st.checkbox("Show reaction mechanisms", True)
    show_energy = st.checkbox("Show energy calculations", True)

# Define the beta-oxidation steps
steps = [
    {
        "name": "Activation",
        "description": "Conversion to Acyl-CoA",
        "enzyme": "Acyl-CoA synthetase",
        "reactants": ["Fatty acid", "ATP", "CoA-SH"],
        "products": ["Acyl-CoA", "AMP", "PPi"],
        "energy_cost": 2,
        "mechanism": "activation_mechanism.png"
    },
    {
        "name": "Transport",
        "description": "Carnitine shuttle",
        "enzyme": "Carnitine palmitoyltransferase",
        "reactants": ["Acyl-CoA", "Carnitine"],
        "products": ["Acyl-carnitine", "CoA-SH"],
        "energy_cost": 0,
        "mechanism": "transport_mechanism.png"
    },
    {
        "name": "1st Dehydrogenation",
        "description": "Form trans-Δ²-enoyl-CoA",
        "enzyme": "Acyl-CoA dehydrogenase",
        "reactants": ["Acyl-CoA", "FAD"],
        "products": ["trans-Δ²-enoyl-CoA", "FADH₂"],
        "energy_cost": 0,
        "energy_gain": 1.5,
        "mechanism": "dehydrogenation1.png"
    },
    {
        "name": "Hydration",
        "description": "Form L-β-hydroxyacyl-CoA",
        "enzyme": "Enoyl-CoA hydratase",
        "reactants": ["trans-Δ²-enoyl-CoA", "H₂O"],
        "products": ["L-β-hydroxyacyl-CoA"],
        "energy_cost": 0,
        "mechanism": "hydration.png"
    },
    {
        "name": "2nd Dehydrogenation",
        "description": "Form β-ketoacyl-CoA",
        "enzyme": "β-hydroxyacyl-CoA dehydrogenase",
        "reactants": ["L-β-hydroxyacyl-CoA", "NAD⁺"],
        "products": ["β-ketoacyl-CoA", "NADH + H⁺"],
        "energy_cost": 0,
        "energy_gain": 2.5,
        "mechanism": "dehydrogenation2.png"
    },
    {
        "name": "Thiolysis",
        "description": "Cleavage to Acetyl-CoA",
        "enzyme": "β-ketothiolase",
        "reactants": ["β-ketoacyl-CoA", "CoA-SH"],
        "products": ["Acetyl-CoA", "Acyl-CoA (n-2)"],
        "energy_cost": 0,
        "energy_gain": 10,
        "mechanism": "thiolysis.png"
    }
]

# Create spiral visualization
def create_spiral(fatty_acid_length=16):
    theta = np.linspace(0, 2*np.pi*(fatty_acid_length//2 - 1), 100*(fatty_acid_length//2 - 1))
    z = np.linspace(0, fatty_acid_length, len(theta))
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter3d(
        x=np.cos(theta),
        y=np.sin(theta),
        z=z,
        mode='lines',
        line=dict(color='blue', width=4),
        name='Carbon chain'
    ))
    
    for i in range(fatty_acid_length//2 - 1):
        angle = 2*np.pi*i
        fig.add_trace(go.Scatter3d(
            x=[np.cos(angle)],
            y=[np.sin(angle)],
            z=[2*i + 2],
            mode='markers+text',
            marker=dict(size=8, color='red'),
            text=[f"Cycle {i+1}"],
            textposition="middle center",
            name=f"Cycle {i+1}"
        ))
    
    fig.update_layout(
        scene=dict(
            xaxis_title='',
            yaxis_title='',
            zaxis_title='Carbon Chain Length',
            camera=dict(eye=dict(x=1.5, y=1.5, z=0.8))
        ),
        title="Lynen's Spiral: β-Oxidation Pathway",
        height=700
    )
    return fig

# Create cycle diagram with better formatting
def draw_cycle_diagram(current_step):
    img = Image.new("RGB", (1000, 800), "white")
    draw = ImageDraw.Draw(img)
    
    # Try to load a nicer font
    try:
        font = ImageFont.truetype("arial.ttf", 24)
        small_font = ImageFont.truetype("arial.ttf", 18)
    except:
        font = ImageFont.load_default()
        small_font = ImageFont.load_default()
    
    center_x, center_y = 500, 400
    radius = 300
    
    # Draw the circular path
    draw.ellipse([(center_x-radius, center_y-radius), (center_x+radius, center_y+radius)], 
                outline="#555555", width=3)
    
    # Define molecules with positions
    molecules = [
        {"name": "Acyl-CoA", "formula": "R-CH₂-CH₂-C~S-CoA", "color": "#4E79A7", "angle": -np.pi/2},
        {"name": "trans-Δ²-Enoyl-CoA", "formula": "R-CH=CH-C~S-CoA", "color": "#F28E2B", "angle": np.pi/10},
        {"name": "L-β-Hydroxyacyl-CoA", "formula": "R-CH(OH)-CH₂-C~S-CoA", "color": "#E15759", "angle": 3*np.pi/5},
        {"name": "β-Ketoacyl-CoA", "formula": "R-CO-CH₂-C~S-CoA", "color": "#76B7B2", "angle": 7*np.pi/5},
        {"name": "Acetyl-CoA", "formula": "CH₃-CO~S-CoA", "color": "#59A14F", "angle": 11*np.pi/6}
    ]
    
    # Draw molecules and arrows
    for i, mol in enumerate(molecules):
        x = center_x + radius * np.cos(mol["angle"])
        y = center_y + radius * np.sin(mol["angle"])
        
        # Highlight current step
        outline = "red" if i == current_step else "black"
        text_color = "black"
        
        # Draw connecting lines (simplified arrows)
        next_angle = molecules[(i+1)%len(molecules)]["angle"]
        draw.line([(x, y), 
                  (center_x + radius*0.9*np.cos(next_angle), 
                   center_y + radius*0.9*np.sin(next_angle))],
                 fill="gray", width=2)
        
        # Draw molecule box
        box_pad = 20
        text_width = font.getlength(mol["formula"])
        draw.rectangle([(x-box_pad, y-box_pad), 
                       (x+text_width+box_pad, y+60+box_pad)], 
                      fill="white", outline=outline, width=2)
        
        # Draw molecule
        draw.text((x, y), mol["formula"], fill=text_color, font=font)
        draw.text((x, y+30), mol["name"], fill=mol["color"], font=small_font)
        
    return img

# Energy calculation
def calculate_energy(fatty_acid_length, unsaturated=False):
    cycles = (fatty_acid_length // 2) - 1
    acetyl_coa = fatty_acid_length // 2
    
    total_energy = {
        "FADH₂": cycles * 1,
        "NADH": cycles * 1,
        "Acetyl-CoA": acetyl_coa,
        "ATP Cost": 2,
        "Total ATP": cycles * (1.5 + 2.5) + acetyl_coa * 10 - 2
    }
    
    if unsaturated:
        total_energy["FADH₂"] -= 1
        total_energy["Total ATP"] -= 1.5
    
    return total_energy

# Main display with tabs
tab1, tab2 = st.tabs(["Lynen Spiral", "Reaction Cycle"])

with tab1:
    st.header("Lynen Spiral Visualization")
    spiral_fig = create_spiral(16 if fatty_acid != "Custom" else custom_length)
    st.plotly_chart(spiral_fig, use_container_width=True)
    
    st.subheader("Energy Balance")
    energy_data = calculate_energy(16 if fatty_acid != "Custom" else custom_length,
                                 fatty_acid.startswith("Oleic") or (fatty_acid == "Custom" and custom_unsat > 0))
    
    energy_df = pd.DataFrame({
        "Molecule": ["FADH₂", "NADH", "Acetyl-CoA", "ATP Cost", "Net ATP"],
        "Quantity": [energy_data["FADH₂"], energy_data["NADH"], 
                    energy_data["Acetyl-CoA"], energy_data["ATP Cost"], 
                    energy_data["Total ATP"]],
        "ATP Value": [1.5, 2.5, 10, -1, ""]
    })
    
    st.dataframe(energy_df, hide_index=True, use_container_width=True)

with tab2:
    st.header("β-Oxidation Reaction Cycle")
    
    # Initialize session state
    if 'current_step' not in st.session_state:
        st.session_state.current_step = 0
    
    # Create columns for diagram and details
    col1, col2 = st.columns([2, 1])
    
    with col1:
        cycle_img = draw_cycle_diagram(st.session_state.current_step)
        st.image(cycle_img, use_container_width=True)
    
    with col2:
        step = steps[st.session_state.current_step]
        st.markdown(f"### {step['name']}")
        st.markdown(f"**Enzyme:** {step['enzyme']}")
        
        st.markdown("#### Reactants:")
        for reactant in step["reactants"]:
            st.markdown(f"- {reactant}")
            
        st.markdown("#### Products:")
        for product in step["products"]:
            st.markdown(f"- {product}")
        
        if show_energy:
            st.markdown("#### Energy:")
            if 'energy_cost' in step and step['energy_cost'] > 0:
                st.markdown(f"- Cost: {step['energy_cost']} ATP")
            if 'energy_gain' in step and step['energy_gain'] > 0:
                st.markdown(f"- Produced: {step['energy_gain']} ATP")
        
        # Navigation controls
        st.markdown("---")
        cols = st.columns(3)
        with cols[0]:
            if st.button("⏮️ Previous"):
                if st.session_state.current_step > 0:
                    st.session_state.current_step -= 1
                    st.rerun()
        with cols[1]:
            if st.button("⏹️ Reset"):
                st.session_state.current_step = 0
                st.rerun()
        with cols[2]:
            if st.button("⏭️ Next"):
                if st.session_state.current_step < len(steps)-1:
                    st.session_state.current_step += 1
                    st.rerun()
        
        # Auto-play checkbox
        auto_play = st.checkbox("Auto-play animation")
        if auto_play:
            st.session_state.current_step = (st.session_state.current_step + 1) % len(steps)
            st.rerun()

# For demonstration - in production you would use actual mechanism images
def create_placeholder_image(text):
    img = Image.new("RGB", (400, 200), "#F0F0F0")
    draw = ImageDraw.Draw(img)
    draw.text((50, 80), text, fill="black")
    return img

# Reaction mechanism display (if enabled)
if show_mechanism:
    st.markdown("---")
    st.subheader("Reaction Mechanism")
    mechanism_img = create_placeholder_image(f"Mechanism for {steps[st.session_state.current_step]['name']}")
    st.image(mechanism_img)