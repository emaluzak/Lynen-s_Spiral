import streamlit as st
import py3Dmol
from stmol import showmol

st.title("3D Lynen Spiral Visualization")

# Get data from previous page
fa_data = st.session_state.get('fa_data', {})

# Create 3D visualization
view = py3Dmol.view(width=800, height=600)
# Add your 3D spiral visualization code here
showmol(view, height=600)

# Steps of beta-oxidation
steps = ["Activation", "First oxidation", "Hydration", 
         "Second oxidation", "Thiolysis"]

selected_step = st.selectbox("Select a step to see details:", steps)

if selected_step:
    # Show mechanism, molecules, ATP, etc.
    st.subheader(f"Step: {selected_step}")
    # Add detailed information for each step
    
if st.button("View 2D Cycle Representation"):
    st.switch_page("pages/3_📊_2D_Cycle.py")