# app.py
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px

st.set_page_config(
    page_title="Transcript Explorer",
    layout="wide",            # ← use the entire browser width
    #initial_sidebar_state="expanded"
)

# 1) Upload any CSV
uploaded_file = st.sidebar.file_uploader("Upload a CSV file", type=["csv"])

if uploaded_file is not None:
    @st.cache_data(show_spinner=False)
    def load_data(file_obj):
        df = pd.read_csv(file_obj)
        # clean headers: strip whitespace & surrounding single‐quotes
        df.columns = df.columns.str.strip().str.strip("'")
        return df

    df = load_data(uploaded_file)

    # 2) Sidebar controls
    st.sidebar.markdown("### Filters")
    genes = sorted(df["LLgeneSymbol"].unique())
    cells = sorted(df["celltype"].unique())

    selected_genes = st.sidebar.multiselect("Select gene(s)", genes, default=genes[:1])
    selected_cells = st.sidebar.multiselect("Select cell type(s)", cells, default=cells)

    # 2b) Transcript selector (populated after gene/cell filters)
    #    We’ll set default to _all_ transcripts present in the interim filter
    interim = df
    if selected_genes:
        interim = interim[interim["LLgeneSymbol"].isin(selected_genes)]
    if selected_cells:
        interim = interim[interim["celltype"].isin(selected_cells)]

    transcripts = sorted(interim["transcript_id"].unique())
    selected_transcripts = st.sidebar.multiselect(
        "Select transcript(s)",
        transcripts,
        default=transcripts
    )

    st.sidebar.markdown("### Plot settings")
    y_measure = st.sidebar.selectbox(
        "Y‐axis:",
        options=["IsoPct", "TPM", "FPKM"],
        index=0
    )
    log2_scale = st.sidebar.checkbox("Log₂ scale (log2(x+1))", value=False)

    # 3) Apply all filters
    filtered = df
    if selected_genes:
        filtered = filtered[filtered["LLgeneSymbol"].isin(selected_genes)]
    if selected_cells:
        filtered = filtered[filtered["celltype"].isin(selected_cells)]
    if selected_transcripts:
        filtered = filtered[filtered["transcript_id"].isin(selected_transcripts)]

    # 4) Prepare Y values
    if log2_scale:
        filtered = filtered.assign(
            __y = np.log2(filtered[y_measure] + 1)
        )
        y_field = "__y"
        y_label = f"log₂({y_measure} + 1)"
    else:
        y_field = y_measure
        y_label = y_measure

      # 5) Plot + Table
    if not filtered.empty:
        # BOX + POINTS by celltype, faceted by gene in a single column
        fig = px.box(
            filtered,
            x="transcript_id",
            y=y_field,
            color="celltype",
            points=False,               # show individual observations
            facet_col="LLgeneSymbol",
            facet_col_wrap=1,
            facet_row_spacing=0.15,
            height=700,
            #width=900,
            labels={y_field: y_label, "celltype": "Cell Type"}
        )
        # rotate x-labels if needed
        fig.update_xaxes(tickangle=-45, matches=None)
        st.plotly_chart(fig, use_container_width=True)

        st.subheader("Data for selected transcript(s)")
        cols = ["transcript_id", "LLgeneSymbol", "celltype", y_measure]
        if log2_scale:
            cols += ["__y"]
        st.dataframe(filtered[cols], use_container_width=True)

    else:
        st.info("No data to display for those selections.")

else:
    st.warning("Please upload a CSV file to get started.")
