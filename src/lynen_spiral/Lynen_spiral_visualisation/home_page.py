import streamlit as st

def show_home_page():
    """Displays the welcome/home page content"""
    st.title("**Welcome to the Fatty Acid Metabolism Visualizer**")
    st.markdown("---")

    col1, col2 = st.columns([1, 1])

    with col1:
        st.header(" What Is a Fatty Acid?")
        st.markdown("""
        A **fatty acid** is a type of fat molecule. It has:
        
        - A **long chain** of carbon (C) and hydrogen (H) atoms — this is the "tail."
        - A special end called a **carboxyl group (COOH)** — this is the "acid head."
        """)


    with col2:
        st.header(" Why Are Fatty Acids Important?")
        st.markdown("""
        Fatty acids serve **three key roles** in the body:

        1. **Fuel for Energy** – especially during fasting or exercise.
        2. **Cell Structure** – they form flexible cell membranes.
        3. **Messengers** – they act as chemical signals in processes like inflammation.
        """)

    st.markdown("---")

    col3, col4 = st.columns([1, 1])

    with col3:
        st.header(" Fuel for Energy")
        st.markdown("""
        The body burns fatty acids to produce energy, especially during:

        - Fasting
        - Prolonged exercise
        - Low-carb diets

       When fatty acids are broken down (in a process called β-oxidation, also known as the Lynen spiral), they are converted into acetyl-CoA, a high-energy molecule that then enters the Krebs cycle to be further oxidized for energy.
        **Each round of β-oxidation produces:**

        - 1 × Acetyl-CoA
        - 1 × FADH₂
        - 1 × NADH

        Breaking down a 16-carbon fatty acid like palmitate yields over **100 ATP**!

        """)
        st.image(r"src\lynen_spiral\Lynen_spiral_visualisation\coa.jpeg", caption="Acetyl-CoA Molecule", use_container_width=True)
        st.image(r"src\lynen_spiral\Lynen_spiral_visualisation\FAD.jpg", caption="FAD Molecule", use_container_width = True)

    with col4:
        st.header(" What Are NADH and FADH₂?")
        st.markdown("""
        **NADH** (*Nicotinamide Adenine Dinucleotide*) and **FADH₂** (*Flavin Adenine Dinucleotide*) are molecules that **carry high-energy electrons**.

        We can think of them as **rechargeable batteries** :

        - They collect energy from β-oxidation, glycolysis, and the Krebs cycle.
        - Then they **deliver electrons** to the **Electron Transport Chain (ETC)** in the mitochondria.
        - The ETC uses these electrons to generate **ATP**, the main energy currency of your body.

        **Energy yield:**

        - NADH → ~2.5 ATP
        - FADH₂ → ~1.5 ATP

        Without them:
        - Cells couldn't efficiently generate ATP.
        - Vital functions like movement, thinking, and healing would slow down.
        """)
        st.image(r"src\lynen_spiral\Lynen_spiral_visualisation\NADH.jpg", caption="NADH Molecule", use_container_width=True)
    

    st.markdown("---")
    st.subheader(" ATP – The Energy Currency")
    st.markdown("""
    **ATP (Adenosine Triphosphate)** is the final product of cellular respiration — it's the molecule the cells use to perform work.

    - Muscle contraction
    - Brain signaling
    - Cellular repair

    All rely on **ATP**, which is made using the energy stored in NADH and FADH₂.

    """)
    st.image(r"src\lynen_spiral\Lynen_spiral_visualisation\atp.jpg", caption="ATP Molecule", use_container_width=True)
