import streamlit as st
from pathlib import Path
import os

st.title("MD Data Folder Selector")

# Initialize session state
if 'data_path' not in st.session_state:
    st.session_state.data_path = ""
if 'variants_reps' not in st.session_state:
    st.session_state.variants_reps = {}

data_path = st.text_input(
    "Enter the full path to your /data folder:",
    value=st.session_state.data_path,
    help="Paste or type the absolute path to the data folder, e.g., /home/user/project/data"
)

if data_path and data_path != st.session_state.data_path:
    data_folder = Path(data_path)
    if data_folder.exists() and data_folder.is_dir():
        st.session_state.data_path = data_path
        # Find variants (immediate subdirs)
        variants = [v.name for v in data_folder.iterdir() if v.is_dir()]
        st.session_state.variants_reps = {}
        for var in variants:
            var_path = data_folder / var
            reps = [r.name for r in var_path.iterdir() if r.is_dir()]
            st.session_state.variants_reps[var] = reps
        st.success(f"Found data in {data_path}")
        st.rerun()
    else:
        st.error(f"Path does not exist or is not a directory: {data_path}")

if st.session_state.variants_reps:
    st.subheader("Found Variants and Replicates")
    for variant, reps in st.session_state.variants_reps.items():
        rep_str = ", ".join(reps)
        st.write(f"**{variant}:** {rep_str}")
    st.info(f"Total variants: {len(st.session_state.variants_reps)}, Total replicates: {sum(len(reps) for reps in st.session_state.variants_reps.values())}")