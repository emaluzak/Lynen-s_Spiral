import streamlit as st

def show_home_page():
    """Displays the welcome/home page content"""
    st.title("🧬 Welcome to the Fatty Acid Metabolism Visualizer")
    st.markdown("---")

    col1, col2 = st.columns([1, 1])

    with col1:
        st.header("🔍 What Is a Fatty Acid?")
        st.write("""
        A **fatty acid** is a type of fat molecule. It has:

        - A **long chain** of carbon (C) and hydrogen (H) atoms — this is the "tail."
        - A special end called a **carboxyl group (COOH)** — this is the "acid head."
        """)
        st.image("/Users/schultheis/project_practical_programming/Lynen_spiral_visualisation_8/fatty_acid_representation.jpg", caption="This is a representation of the general structure of a fatty acid ", use_container_width=True)


    with col2:
        st.header("🧠 Why Are Fatty Acids Important?")
        st.write("""
        Fatty acids do **three big jobs** in your body:

        1. 🏃‍♀️ **Fuel for Energy** – Your body burns them to produce energy, especially when fasting or exercising.
        2. 🧱 **Cell Structure** – They help form cell membranes and keep them flexible.
        3. 📬 **Messengers** – Some fatty acids act like signals to help control inflammation and healing.
        """)

    st.markdown("---")

    col3, col4 = st.columns([1, 1])

    with col3:
        st.header("🍔 Where Do You Find Fatty Acids?")
        st.write("""
        You eat fatty acids every day in foods like:

        - Butter, cheese, and animal fat
        - Olive oil, sunflower oil
        - Avocados, nuts, and seeds
        - Fish and meat

        Your body breaks down these fats into fatty acids for energy or storage.
        """)

    with col4:
        st.header("🧪 Are All Fatty Acids the Same?")
        st.write("""
        No! There are different types:

        | Type | Description | Example |
        |------|-------------|---------|
        | **Saturated** | No double bonds, solid at room temp | Butter, lard |
        | **Unsaturated** | Has double bonds, liquid at room temp | Olive oil, fish oil |

        Unsaturated fats can be **cis** (natural bend) or **trans** (artificial).
        """)

    st.markdown("---")

 