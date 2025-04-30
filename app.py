# app.py
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

st.set_page_config(page_title="Transcript Dashboard", layout="wide")

# ──────────────────────────────────────────────────────────────────────────────
# Sidebar: upload & select
# ──────────────────────────────────────────────────────────────────────────────
# 1) Upload expression CSV
expr_file = st.sidebar.file_uploader("Upload expression CSV", type=["csv"], key="expr")
if expr_file:
    @st.cache_data(show_spinner=False)
    def load_expr(f):
        df = pd.read_csv(f)
        df.columns = df.columns.str.strip().str.strip("'")
        return df
    expr_df = load_expr(expr_file)

    # 2) Gene selector
    genes = sorted(expr_df["LLgeneSymbol"].unique())
    selected_genes = st.sidebar.multiselect("Select gene(s)", genes, default=genes[:1])

    # 3) Transcript selector (shared)
    transcripts_all = expr_df.query("LLgeneSymbol in @selected_genes")["transcript_id"].unique()
    transcripts = sorted(transcripts_all)
    selected_transcripts = st.sidebar.multiselect(
        "Select transcript(s)", transcripts, default=transcripts
    )
else:
    expr_df = None
    selected_genes = []
    selected_transcripts = []

# 4) Upload GTF
gtf_file = st.sidebar.file_uploader("Upload GTF file", type=["gtf","gtf.gz"], key="gtf")
if gtf_file:
    @st.cache_data(show_spinner=False)
    def load_gtf(buf):
        g = pd.read_csv(
            buf, sep="\t", comment="#", header=None,
            names=["seqname","src","feature","start","end","score","strand","frame","attributes"],
            dtype={"seqname": str}
        )
        g = g.query("feature=='exon'")
        g["transcript_id"] = g["attributes"].str.extract(r'transcript_id "([^"]+)"')[0]
        return g[["seqname","start","end","strand","transcript_id"]]
    exon_df = load_gtf(gtf_file)
else:
    exon_df = None

# ──────────────────────────────────────────────────────────────────────────────
# Tabs
# ──────────────────────────────────────────────────────────────────────────────
tab1, tab2 = st.tabs(["Expression Explorer", "Transcript Viewer"])

# ──────────────────────────────────────────────────────────────────────────────
with tab1:
    st.header("Expression Explorer")
    if expr_df is None:
        st.info("Please upload an expression CSV.")
    else:
        # filter by genes + transcripts
        df = expr_df
        if selected_genes:
            df = df[df["LLgeneSymbol"].isin(selected_genes)]
        if selected_transcripts:
            df = df[df["transcript_id"].isin(selected_transcripts)]

        # plot settings
        y_measure = st.sidebar.selectbox("Y-axis:", ["IsoPct","TPM","FPKM"], index=0)
        log2 = st.sidebar.checkbox("Log₂ scale (log2(x+1))", value=False)
        if log2:
            df = df.assign(__y=np.log2(df[y_measure]+1))
            y_field, y_label = "__y", f"log₂({y_measure}+1)"
        else:
            y_field, y_label = y_measure, y_measure

        if not df.empty:
            fig = px.box(
                df,
                x="transcript_id", y=y_field,
                color="celltype", points=False,
                facet_col="LLgeneSymbol", facet_col_wrap=1,
                facet_row_spacing=0.15,
                height=700, width=1200,
                labels={y_field: y_label, "celltype":"Cell Type"}
            )
            fig.update_xaxes(tickangle=-45, matches=None)
            st.plotly_chart(fig, use_container_width=True)

            st.subheader("Expression Data Table")
            cols = ["transcript_id","LLgeneSymbol","celltype",y_measure]
            if log2: cols.append("__y")
            st.dataframe(df[cols], use_container_width=True)
        else:
            st.info("No data for those selections.")

with tab2:
    st.header("Transcript Viewer")
    if exon_df is None:
        st.info("Please upload a GTF file.")
    elif expr_df is None:
        st.info("Please upload expression CSV first.")
    elif not selected_genes:
        st.info("Select at least one gene in the sidebar.")
    elif not selected_transcripts:
        st.info("Select at least one transcript in the sidebar.")
    else:
        # filter exons for chosen transcripts
        ed = exon_df[exon_df["transcript_id"].isin(selected_transcripts)]

        # group by gene to plot separate panels
        for gene in selected_genes:
            # transcripts of this gene that were selected
            gene_txs = [tx for tx in selected_transcripts
                        if tx in expr_df.query("LLgeneSymbol==@gene")["transcript_id"].values]
            if not gene_txs:
                st.write(f"**{gene}**: no transcripts selected.")
                continue

            st.subheader(f"Gene: {gene}")
            tx_exon_df = ed[ed["transcript_id"].isin(gene_txs)]

            # find shared exons
            exon_sets = {
                tx: set(map(tuple, grp[["start","end"]].values))
                for tx, grp in tx_exon_df.groupby("transcript_id")
            }
            shared = set.intersection(*exon_sets.values())

            # chromosome label
            chrs = set(tx_exon_df["seqname"])
            x_title = (f"{chrs.pop()})"
                       if len(chrs)==1 else "Genomic coordinate")

            # build figure
            fig = go.Figure()
            y_gap = 1.0
            for i, tx in enumerate(gene_txs):
                sub = tx_exon_df[tx_exon_df["transcript_id"]==tx]
                s, e = sub["start"].min(), sub["end"].max()
                strand = sub["strand"].iat[0]
                y = i * y_gap

                # intron line
                fig.add_trace(go.Scatter(
                    x=[s,e], y=[y,y], mode="lines",
                    line=dict(color="black", width=1),
                    showlegend=False, hoverinfo="skip"
                ))

                # strand arrow
                if strand=="+":
                    fig.add_annotation(
                        x=e, y=y,
                        ax=e-(e-s)*0.05, ay=y,
                        showarrow=True, arrowhead=3, arrowwidth=1
                    )
                else:
                    fig.add_annotation(
                        x=s, y=y,
                        ax=s+(e-s)*0.05, ay=y,
                        showarrow=True, arrowhead=3, arrowwidth=1
                    )

                # exon boxes
                for _, exon in sub.iterrows():
                    col = "steelblue" if (exon.start,exon.end) in shared else "red"
                    fig.add_shape(
                        type="rect",
                        x0=exon.start, x1=exon.end,
                        y0=y-0.2,     y1=y+0.2,
                        line=dict(color=col),
                        fillcolor=col
                    )

                # label each transcript
                fig.add_annotation(
                    x=s, y=y+0.4,
                    text=tx,
                    showarrow=False,
                    align="left",
                    font=dict(size=10)
                )

            fig.update_yaxes(
                showticklabels=False,
                range=[-0.5, len(gene_txs)*y_gap - 0.5]
            )
            fig.update_xaxes(title_text=x_title)
            fig.update_layout(
                height=max(200, len(gene_txs)*100),
                margin=dict(l=20, r=20, t=20, b=20)
            )
            st.plotly_chart(fig, use_container_width=True)
