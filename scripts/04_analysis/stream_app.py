import streamlit as st
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode

st.set_page_config(layout="wide")

# Base directory
BASE_DIR = "../../data/1bf5"

# Discover available MUT directories
MUT_DIRS = sorted([d for d in os.listdir(BASE_DIR) 
                  if os.path.isdir(os.path.join(BASE_DIR, d)) and not d.startswith('.')])

# Let user select MUT value(s)
st.sidebar.title("Configuration")
selected_muts = st.sidebar.multiselect("Select MUT Values", MUT_DIRS, default=[MUT_DIRS[0]] if MUT_DIRS else [])

if not selected_muts:
    st.warning("Please select at least one MUT value from the sidebar.")
    st.stop()

### SECTIONS
SECTIONS = ["QC", "PCA", "Hydrogen Bonds", "Test zone"]

st.title(f"Analysis Dashboard for {', '.join(selected_muts)}")

section = st.selectbox("Select Section", SECTIONS, index=0)

### Section: QC
if section == "QC":
    # Plot types and their configurations
    PLOT_TYPES = {
        "RMSD": ("rmsd.csv", "Time (ps)", ["Backbone", 'C-alphas', 'Protein'], "RMSD vs Time"),
        "Radius of Gyration - Axis": ("radius_axis.csv", "Frame", ['all', 'x-axis', 'y-axis', 'z-axis'], "Radius of Gyration vs Time"),
        "RMSF - Protein Domains": ("rmsf_prot.csv", "Residue", ['RMSF_A', 'RMSF_B'], "RMSF Protein"),
        "RMSF - DNA Residues": ("rmsf_dna.csv", "Residue", ['RMSF_C', 'RMSF_D'], "RMSF DNA"),
        "Backbone RMSD vs Radius of Gyration": (("rmsd.csv", "radius_axis.csv"), "Backbone", "all", "Backbone RMSD vs Radius of Gyration")
    }

    st.title(f"QC Plots for {', '.join(selected_muts)}")

    # Dropdown to select plot type
    plot_type = st.selectbox("Select Plot Type", list(PLOT_TYPES.keys()), index=0)

    # RMSD-specific checkboxes
    if plot_type == "RMSD":
        st.write("Show:")
        col1, col2, col3 = st.columns(3)
        with col1:
            bb = st.checkbox("Backbone", value=True, key="bb_rmsd")
        with col2:
            ca = st.checkbox("C-alphas", value=True, key="ca_rmsd")
        with col3:
            prot = st.checkbox("Protein", value=True, key="prot_rmsd")
        col_y = [label for label, checked in [("Backbone", bb), ("C-alphas", ca), ("Protein", prot)] if checked]
    else:
        col_y = None

    # Display plots for all selected MUTs and their reps
    for mut in selected_muts:
        SAMPLE_DIRS = sorted([d for d in os.listdir(os.path.join(BASE_DIR, mut))
                            if os.path.isdir(os.path.join(BASE_DIR, mut, d)) and d.startswith("rep")])
        
        st.subheader(f"{plot_type} for {mut} - All Reps")
        row1_col1, row1_col2 = st.columns(2)
        row2_col1, row2_col2 = st.columns(2)

        # Get plot data
        file, x_col, y_col, plot_title = PLOT_TYPES[plot_type]

        # Apply RMSD filtering if applicable
        if plot_type == "RMSD":
            if not col_y:
                st.warning("No data selected. Select at least one property to plot.")
                continue
            y_col_to_plot = [col for col in y_col if col in col_y]
            yaxis_title = "RMSD (Å)"
        else:
            y_col_to_plot = y_col
            if plot_type == "Radius of Gyration - Axis":
                yaxis_title = "Radius of Gyration (Å)"
            elif plot_type == "RMSF - Protein Domains":
                yaxis_title = "RMSF (Å)"
            elif plot_type == "RMSF - DNA Residues":
                yaxis_title = "RMSF (Å)"
            elif plot_type == "Backbone RMSD vs Radius of Gyration":
                yaxis_title = "Radius of Gyration (Å)"
            else:
                yaxis_title = y_col if not isinstance(y_col, list) else "Value"

        # Compute global y-axis range
        y_min, y_max = float('inf'), float('-inf')
        if not isinstance(file, tuple):
            for rep in SAMPLE_DIRS[:4]:
                file_path = Path(BASE_DIR) / mut / rep / "analysis" / "001_qc" / file
                if file_path.exists():
                    df = pd.read_csv(file_path)
                    for col in y_col_to_plot:
                        if col in df.columns:
                            y_min = min(y_min, df[col].min())
                            y_max = max(y_max, df[col].max())
        else:
            for rep in SAMPLE_DIRS[:4]:
                rmsd_file, radius_file = file
                rmsd_path = Path(BASE_DIR) / mut / rep / "analysis" / "001_qc" / rmsd_file
                radius_path = Path(BASE_DIR) / mut / rep / "analysis" / "001_qc" / radius_file
                if rmsd_path.exists() and radius_path.exists():
                    df_radius = pd.read_csv(radius_path)
                    y_min = min(y_min, df_radius['all'].min())
                    y_max = max(y_max, df_radius['all'].max())
        padding = (y_max - y_min) * 0.05
        y_min -= padding
        y_max += padding

        # Domain annotations for RMSF - Protein Domains
        domain_shapes = [
            dict(type="rect", x0=137, x1=317, y0=0, y1=1, yref="paper", fillcolor="blue", opacity=0.2, layer="below", name="CCD"),
            dict(type="rect", x0=318, x1=488, y0=0, y1=1, yref="paper", fillcolor="red", opacity=0.2, layer="below", name="DBD"),
            dict(type="rect", x0=489, x1=576, y0=0, y1=1, yref="paper", fillcolor="green", opacity=0.2, layer="below", name="LD"),
            dict(type="rect", x0=577, x1=683, y0=0, y1=1, yref="paper", fillcolor="yellow", opacity=0.2, layer="below", name="SH2"),
            dict(type="rect", x0=684, x1=708, y0=0, y1=1, yref="paper", fillcolor="pink", opacity=0.2, layer="below", name="Tail")
        ]

        # Plot for each rep
        for i, rep in enumerate(SAMPLE_DIRS[:4]):
            data_dir = Path(BASE_DIR) / mut / rep / "analysis" / "001_qc"

            if isinstance(file, tuple):
                rmsd_file, radius_file = file
                rmsd_path = data_dir / rmsd_file
                radius_path = data_dir / radius_file
                if rmsd_path.exists() and radius_path.exists():
                    df_rmsd = pd.read_csv(rmsd_path)
                    df_radius = pd.read_csv(radius_path)
                    if len(df_rmsd) == len(df_radius):
                        df = pd.DataFrame({
                            "Backbone_RMSD": df_rmsd['Backbone'],
                            "Radius_All": df_radius['all']
                        })
                        fig = px.scatter(df, x="Backbone_RMSD", y="Radius_All", title=f"{plot_title} - {mut} {rep}", opacity=0.4)
                        fig.update_layout(
                            xaxis_title="Backbone RMSD (Å)",
                            yaxis_title=yaxis_title,
                            yaxis_range=[y_min, y_max]
                        )
                    else:
                        st.error(f"Data length mismatch for {mut} {rep}: {rmsd_file} ({len(df_rmsd)}) vs {radius_file} ({len(df_radius)})")
                        continue
                else:
                    st.warning(f"Files not found for {mut} {rep}: {rmsd_path} or {radius_path}")
                    continue
            else:
                file_path = data_dir / file
                if file_path.exists():
                    df = pd.read_csv(file_path)
                    if plot_type == "RMSF - Protein Domains":
                        df[x_col] = df[x_col] + 135
                    if isinstance(y_col_to_plot, list):
                        fig = px.line(df, x=x_col, y=y_col_to_plot, title=f"{plot_title} - {mut} {rep}")
                    elif x_col == "RMSD":
                        fig = px.scatter(df, x=x_col, y=y_col_to_plot, title=f"{plot_title} - {mut} {rep}")
                    else:
                        fig = px.line(df, x=x_col, y=y_col_to_plot, title=f"{plot_title} - {mut} {rep}")
                    
                    if plot_type == "RMSF - Protein Domains":
                        fig.update_layout(shapes=domain_shapes, showlegend=True)
                        for shape in domain_shapes:
                            fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines", line=dict(color=shape["fillcolor"], width=0),
                                                   fill="toself", opacity=0.2, name=shape["name"], legendgroup=shape["name"], showlegend=True))

                    fig.update_layout(xaxis_title=x_col, yaxis_title=yaxis_title, yaxis_range=[y_min, y_max])
                else:
                    st.warning(f"Data file not found for {mut} {rep}: {file_path}")
                    continue

            # Assign to 2x2 grid
            if i == 0:
                col = row1_col1
            elif i == 1:
                col = row1_col2
            elif i == 2:
                col = row2_col1
            elif i == 3:
                col = row2_col2
            col.plotly_chart(fig, use_container_width=True)

    st.write(f"Data sourced from: {BASE_DIR}/<mut>/<rep>/analysis/001_qc")

# Section: PCA
elif section == "PCA":
    st.subheader("Principal Component Analysis")
    for mut in selected_muts:
        st.write(f"PCA plots and tables for {mut} go here.")

# Section: Hydrogen Bonds
elif section == "Hydrogen Bonds":
    st.subheader("Hydrogen Bonds Analysis")
    
    VIEW_TYPES = {
        "Filtered Plot": ("Filtered Plot", None),
        "GMX External Table": ("hbonds_gmx_ext.csv", "GMX External Hydrogen Bonds"),
        "MDA External Table": ("hbonds_mda_ext.csv", "MDA External Hydrogen Bonds")
    }
    view_type = st.selectbox("Select View Type", list(VIEW_TYPES.keys()), index=0, key="hbond_view_select")
    file, title = VIEW_TYPES[view_type]
    
    for mut in selected_muts:
        SAMPLE_DIRS = sorted([d for d in os.listdir(os.path.join(BASE_DIR, mut))
                            if os.path.isdir(os.path.join(BASE_DIR, mut, d)) and d.startswith("rep")])
        
        col1, col2, col3 = st.columns(3)
        
        if view_type == "Filtered Plot":
            with st.sidebar:
                st.subheader(f"Filter Options for {mut}")
                sample_file = Path(BASE_DIR) / mut / SAMPLE_DIRS[0] / "analysis" / "002_hbonds" / "hbonds_mda_ext.csv"
                if sample_file.exists():
                    sample_df = pd.read_csv(sample_file)
                    donor_types = sorted(sample_df['donor_type'].unique())
                    acceptor_types = sorted(sample_df['acceptor_type'].unique())
                else:
                    donor_types = ["protein", "DNA"]
                    acceptor_types = ["protein", "DNA"]
                selected_donor_types = st.multiselect(f"Donor Type ({mut})", donor_types, default=donor_types, key=f"donor_{mut}")
                selected_acceptor_types = st.multiselect(f"Acceptor Type ({mut})", acceptor_types, default=acceptor_types, key=f"acceptor_{mut}")
                x_axis = st.selectbox(f"X-Axis ({mut})", ["Residue Number", "Atom Number"], index=0, key=f"xaxis_{mut}")
        
            st.subheader(f"Hydrogen Bonds for {mut}")
            for i, rep in enumerate(SAMPLE_DIRS[:3]):
                data_dir = Path(BASE_DIR) / mut / rep / "analysis" / "002_hbonds"
                file_path = data_dir / "hbonds_mda_ext.csv"
                if file_path.exists():
                    hbonds_df = pd.read_csv(file_path)
                    filtered_df = hbonds_df[
                        (hbonds_df['donor_type'].isin(selected_donor_types)) | 
                        (hbonds_df['acceptor_type'].isin(selected_acceptor_types))
                    ]
                    
                    if x_axis == "Residue Number":
                        df = pd.concat([
                            filtered_df[filtered_df['donor_type'].isin(selected_donor_types)]
                                .assign(Type="Donor")[["donor_res_num", "frame%", "Type"]]
                                .rename(columns={"donor_res_num": "X"}),
                            filtered_df[filtered_df['acceptor_type'].isin(selected_acceptor_types)]
                                .assign(Type="Acceptor")[["acceptor_res_num", "frame%", "Type"]]
                                .rename(columns={"acceptor_res_num": "X"})
                        ])
                        df["X"] = df["X"] + 135
                        plot_title = f"Hydrogen Bonds by Residue - {mut} {rep}"
                        x_label = "Residue Number"
                    else:
                        df = pd.concat([
                            filtered_df[filtered_df['donor_type'].isin(selected_donor_types)]
                                .assign(Type="Donor")[["donor", "frame%", "Type"]]
                                .rename(columns={"donor": "X"}),
                            filtered_df[filtered_df['acceptor_type'].isin(selected_acceptor_types)]
                                .assign(Type="Acceptor")[["acceptor", "frame%", "Type"]]
                                .rename(columns={"acceptor": "X"})
                        ])
                        plot_title = f"Hydrogen Bonds by Atom - {mut} {rep}"
                        x_label = "Atom Number"
                    
                    fig = px.scatter(df, x="X", y="frame%", color="Type", 
                                   title=plot_title,
                                   color_discrete_map={"Donor": "red", "Acceptor": "blue"},
                                   opacity=0.6)
                    fig.update_layout(xaxis_title=x_label, yaxis_title="Frame %")
                    
                    if x_axis == "Residue Number":
                        for domain, (start, end, color) in [
                            ("CCD", (136, 317, "blue")), ("DBD", (318, 488, "red")),
                            ("LD", (489, 576, "green")), ("SH2", (577, 683, "yellow")),
                            ("Tail", (684, 708, "pink"))
                        ]:
                            fig.add_vrect(x0=start, x1=end, fillcolor=color, opacity=0.2, line_width=0, annotation_text=domain)
                    
                    cols = [col1, col2, col3]
                    col = cols[i % 3]
                    col.plotly_chart(fig, use_container_width=True)
                else:
                    st.warning(f"Data file not found for {mut} {rep}: {file_path}")
        
        else:
            st.subheader(f"{title} for {mut}")
            for i, rep in enumerate(SAMPLE_DIRS[:3]):
                data_dir = Path(BASE_DIR) / mut / rep / "analysis" / "002_hbonds"
                file_path = data_dir / file
                if file_path.exists():
                    df = pd.read_csv(file_path)
                    gb = GridOptionsBuilder.from_dataframe(df)
                    gb.configure_default_column(sortable=True, filter=True, resizable=True)
                    gb.configure_grid_options(domLayout='normal')
                    grid_options = gb.build()
                    cols = [col1, col2, col3]
                    col = cols[i % 3]
                    col.subheader(f"{title} - {mut} {rep}")
                    AgGrid(df, gridOptions=grid_options, height=300, update_mode=GridUpdateMode.MANUAL,
                          fit_columns_on_grid_load=True)
                else:
                    st.warning(f"Data file not found for {mut} {rep}: {file_path}")
    
    st.write(f"Data sourced from: {BASE_DIR}/<mut>/<rep>/analysis/002_hbonds")

elif section == "Test zone":
    df = pd.read_csv(f'{BASE_DIR}/wt/rep0/analysis/002_hbonds/hbonds_mda_ext.csv')
    # Display the DataFrame in an AgGrid table
    table = AgGrid(df)
    
    # Retrieve the updated data from the table
    updated_data = table['data']
    
    # Create a chart based on the updated data
    fig = px.bar(updated_data, x='donor_res_num_adjusted', y='count')
    
    # Display the chart
    st.plotly_chart(fig)