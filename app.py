from shiny import App, ui, render, reactive, Inputs, Outputs, Session
import pandas as pd
import matplotlib.pyplot as plt
import os
import gseapy as gp
from scipy.stats import pearsonr
import plotly.express as px
from shinywidgets import render_widget, output_widget
import plotly.graph_objects as go
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import io

all_gseapy_libraries = gp.get_library_name()
# Pre-load MSigDB gene sets for both Human and Mouse
gene_sets_human = gp.get_library(name='GO_Biological_Process_2023', organism='Human')
gene_sets_mouse = gp.get_library(name='GO_Biological_Process_2023', organism='Mouse')


# Read the file
##mousediv10 = pd.read_excel("C:\Users\Lenovo\OneDrive\Documentos\Bioinformatics\mouse_div10.xlsx" , index_col=0)
#mousediv4 = pd.read_excel("C:\Users\Lenovo\OneDrive\Documentos\Bioinformatics\mouse_div4.xlsx" , index_col=0)
#humn = pd.read_excel("C:\Users\Lenovo\OneDrive\Documentos\Bioinformatics\human_esc.xlsx" , index_col=0)

app_ui = ui.page_fluid(
    ui.layout_sidebar(
        ui.sidebar(
            ui.input_file("file1", "Upload Dataset 1", accept=[".csv", ".txt", ".xls", ".xlsx"]),
            ui.input_file("file2", "Upload Dataset 2", accept=[".csv", ".txt", ".xls", ".xlsx"]),
            ui.input_numeric("padj_threshold", "padj threshold", 0.05, min=0, max=1, step=0.01),
            ui.input_checkbox("show_axes", "Show axes on V2 Plot", value=True),
            ui.input_select(
                "organism",
                "Select organism for gene sets",
                choices=["Human", "Mouse"],
                selected="Mouse"
            ),
            ui.input_select(
                "gene_set_library",
                "Select gene set library",
                choices=all_gseapy_libraries,
                selected="GO_Biological_Process_2023" if "GO_Biological_Process_2023" in all_gseapy_libraries else all_gseapy_libraries[0]
            ),
            ui.input_selectize(
                "gene_set_term",
                "Select a gene set",
                choices=[],
                multiple=True,
            ),
            ui.input_file("user_gene_set", "Upload gene set", accept=[".csv", ".txt", ".xls", ".xlsx"]),
        ),
        ui.navset_tab(
            ui.nav_panel("V2 Plot (Shape & Highlight)", ui.output_plot("v2_shape_highlight_plot"), ui.download_button("download_v2_shape_highlight_plot", "Download PNG"), ui.output_text("pearson_r_value"), ui.output_text("pearson_r_selected_gene_set"), ui.input_text("x_axis_label", "X Axis Label", value="Dataset 1 log2 fold change"), ui.input_text("y_axis_label", "Y Axis Label", value="Dataset 2 log2 fold change"), ui.input_text("plot_title", "Plot Title", value="V2 Plot (Shape & Highlight)"), ui.input_numeric("point_size", "Point Size", value=40, min=1, max=200, step=1),),
            ui.nav_panel("V2 Plot (Plotly Shape & Highlight)", output_widget("v2_shape_highlight_plotly")),
            ui.nav_panel("Dataset 1 & 2", ui.output_table("preview1"), ui.output_table("preview2")),
            ui.nav_panel(   
                "Gene Lists",
                ui.output_table("head_all_genes"),
                ui.markdown("**MSigDB Mouse Gene Sets (head):**"),
                ui.output_table("head_mouse_gene_sets"),
                ui.markdown("**MSigDB Human Gene Sets (head):**"),
                ui.output_table("head_human_gene_sets"),
                ui.markdown("**Filtered Gene Sets (head):**"),
                ui.output_table("filtered_gene_sets_preview"),
                ui.markdown("**User Uploaded Gene Names:**"),  
                ui.output_table("user_uploaded_gene_names"),   
        ),
        )
    )
)

def server(input: Inputs, output: Outputs, session: Session):

    def read_any_file(file_info):
        if not file_info:
            return None
        path = file_info[0]["datapath"]
        ext = os.path.splitext(path)[1].lower()
        try:
            if ext in [".csv", ".txt"]:
                df = pd.read_csv(path, sep=None, engine="python", index_col=0)
            elif ext in [".xls", ".xlsx"]:
                df = pd.read_excel(path, index_col=0)
            else:
                return None
            # Convert the index (gene names) to uppercase
            df.index = df.index.astype(str).str.upper()
            return df
        except Exception as e:
            print(f"Error reading file {path}: {e}")
            return None

    @reactive.Calc
    def df1():
        return read_any_file(input.file1())

    @reactive.Calc
    def df2():
        return read_any_file(input.file2())
    
    #Creating set of all genes from both datasets
    @reactive.Calc
    def all_genes():
        d1 = df1()
        d2 = df2()
        genes = set()
        if d1 is not None and d1.shape[0] > 0:
            genes.update(pd.Series(d1.index).dropna().astype(str).str.upper().unique())
        if d2 is not None and d2.shape[0] > 0:
            genes.update(pd.Series(d2.index).dropna().astype(str).str.upper().unique())
        return genes
    
    @reactive.Calc
    def gene_set_dict():
        lib = input.gene_set_library()
        org = input.organism()
        try:
            # Try to load the selected library for the selected organism
            return gp.get_library(name=lib, organism=org)
        except Exception as e:
            print(f"Error loading gene set library '{lib}' for organism '{org}': {e}")
            return {}
    
    @output
    @render.table
    def preview1():
        df = df1()
        return df.head(10) if df is not None else None

    @output
    @render.table
    def preview2():
        df = df2()
        return df.head(10) if df is not None else None
    
    #tab for all gene lists - all genes, gpseea mouse and human gene sets
    @output
    @render.table
    def head_all_genes():
        genes = list(all_genes())
        return pd.DataFrame(genes[:10], columns=["Gene"])

    @output
    @render.table
    def head_mouse_gene_sets():
        # Show the first 5 gene sets and their first 5 genes
        data = []
        for i, (k, v) in enumerate(gene_sets_mouse.items()):
            if i >= 5:
                break
            data.append({"Gene Set": k, "Genes": ", ".join(v[:5]) + ("..." if len(v) > 5 else "")})
        return pd.DataFrame(data)

    @output
    @render.table
    def head_human_gene_sets():
        # Show the first 5 gene sets and their first 5 genes
        data = []
        for i, (k, v) in enumerate(gene_sets_human.items()):
            if i >= 5:
                break
            data.append({"Gene Set": k, "Genes": ", ".join(v[:5]) + ("..." if len(v) > 5 else "")})
        return pd.DataFrame(data)
    
    @output
    @render.table
    def user_uploaded_gene_names():
        genes = list(user_uploaded_genes())
        if not genes:
            return pd.DataFrame({"Gene": ["(No genes uploaded)"]})
        return pd.DataFrame(genes, columns=["Gene"])
     
    #1. Filter gene sets to only those where all genes are present in the uploaded data    @reactive.Calc
    def filtered_gene_set_dict():
        genes_in_data = set(all_genes())
        # Include gene sets where AT LEAST ONE gene is present in the data
        return {
            k: v for k, v in gene_set_dict().items()
            if any(g.upper() in genes_in_data for g in v)
        }
    
    @output
    @render.table
    def filtered_gene_sets_preview():
        # Show the first 10 filtered gene sets
        data = []
        for i, (k, v) in enumerate(filtered_gene_set_dict().items()):
            if i >= 10:
                break
            data.append({"Gene Set": k, "Genes": ", ".join(v[:5]) + ("..." if len(v) > 5 else "")})
        return pd.DataFrame(data)
     
    # 2. Update the dropdown to show only these gene sets (limit to 1000)
    @reactive.effect
    def update_gene_set_dropdown():
        gene_set_data = [
            {"Gene Set": k, "Genes": ", ".join(v)}
            for k, v in filtered_gene_set_dict().items()
        ]
        gene_set_df = pd.DataFrame(gene_set_data, columns=["Gene Set", "Genes"])
        gene_set_names = gene_set_df["Gene Set"].tolist()[:1000] if not gene_set_df.empty else []
        ui.update_selectize(
            "gene_set_term",
            choices=gene_set_names,
            session=session,
        )

    #3. Get the selected gene set from the dropdown
    @reactive.Calc
    def selected_gene_set():
        terms = input.gene_set_term()
        if not terms:
            return set()
        # Ensure terms is a list (it will be if multiple=True)
        if isinstance(terms, str):
            terms = [terms]
        genes = set()
        for term in terms:
            genes.update([g.upper() for g in filtered_gene_set_dict().get(term, [])])
        return genes
    
    # Read user-uploaded gene set (first row as header, first column as gene names)
    @reactive.Calc
    def user_uploaded_genes():
        try:
            df = read_any_file(input.user_gene_set())
            if df is None or df.index.empty:
                return set()
            # Use the index for gene names, uppercase for matching
            genes = pd.Series(df.index).dropna().astype(str).str.upper().unique()
            return set(genes)
        except Exception as e:
            print(f"Error reading user gene set: {e}")
            return set()
        
    # Combine dropdown and uploaded gene set for highlighting
    @reactive.Calc
    def highlight_genes():
        dropdown_genes = selected_gene_set()
        uploaded_genes = user_uploaded_genes()
        return dropdown_genes.union(uploaded_genes)
    
    @output
    @render.text
    def pearson_r_selected_gene_set():
        d1, d2 = df1(), df2()
        if d1 is None or d2 is None:
            return "Upload both datasets to calculate Pearson R."
        merged = pd.merge(d1, d2, left_index=True, right_index=True, suffixes=("_1", "_2"))
        highlight = highlight_genes()
        # Only keep genes in the selected gene set
        filtered = merged.loc[merged.index.isin(highlight)].copy()
        if filtered.shape[0] < 2:
            return "Not enough significant genes in the selected gene set."
        if "log2_fold_change_1" not in filtered or "log2_fold_change_2" not in filtered:
            return "log2_fold_change columns not found in both datasets."
        r, p = pearsonr(filtered["log2_fold_change_1"], filtered["log2_fold_change_2"])
        return f"Pearson R (selected gene set): {r:.5f} (p={p:.5g})"
    
    @output
    @render.plot
    def v2_shape_highlight_plot():
        d1, d2 = df1(), df2()
        if d1 is None or d2 is None:
            return
        merged = pd.merge(d1, d2, left_index=True, right_index=True, suffixes=("_1", "_2"))
        merged["significant_1"] = merged["padj_1"] < input.padj_threshold()
        merged["significant_2"] = merged["padj_2"] < input.padj_threshold()
        shape_map = {
            (True, True): "o",      # Circle
            (True, False): "s",     # Square
            (False, True): "^",     # Triangle
            (False, False): "D",    # Diamond
        }
        merged["shape"] = merged[["significant_1", "significant_2"]].apply(lambda x: shape_map[tuple(x)], axis=1)
        highlight = highlight_genes()
        merged["highlight"] = merged.index.isin(highlight)
        fig, ax = plt.subplots()
        point_size = input.point_size() if input.point_size() else 40
        # Plot non-highlighted genes by shape
        for (sig1, sig2), marker in shape_map.items():
            mask = (merged["significant_1"] == sig1) & (merged["significant_2"] == sig2) & (~merged["highlight"])
            ax.scatter(
                merged.loc[mask, "log2_fold_change_1"],
                merged.loc[mask, "log2_fold_change_2"],
                marker=marker,
                c="gray",
                alpha=0.7,
                s=point_size,
                label=f"1:{'Yes' if sig1 else 'No'}, 2:{'Yes' if sig2 else 'No'}"
            )
        # Plot highlighted genes by shape
        for (sig1, sig2), marker in shape_map.items():
            mask = (merged["significant_1"] == sig1) & (merged["significant_2"] == sig2) & (merged["highlight"])
            ax.scatter(
                merged.loc[mask, "log2_fold_change_1"],
                merged.loc[mask, "log2_fold_change_2"],
                marker=marker,
                c="orange",
                edgecolor="black",
                alpha=0.9,
                s=point_size + 10,
                label=None if merged.loc[mask].empty else "Highlighted gene set"
            )
        x = merged["log2_fold_change_1"]
        y = merged["log2_fold_change_2"]
        min_val = min(x.min(), y.min())
        max_val = max(x.max(), y.max())
        ax.set_xlim(min_val, max_val)
        ax.set_ylim(min_val, max_val)
        ax.plot([min_val, max_val], [min_val, max_val], "k--", alpha=0.5)
        handles = [
            mlines.Line2D([], [], color='gray', marker='o', linestyle='None', markersize=8, label='1:Yes, 2:Yes'),
            mlines.Line2D([], [], color='gray', marker='s', linestyle='None', markersize=8, label='1:Yes, 2:No'),
            mlines.Line2D([], [], color='gray', marker='^', linestyle='None', markersize=8, label='1:No, 2:Yes'),
            mlines.Line2D([], [], color='gray', marker='D', linestyle='None', markersize=8, label='1:No, 2:No'),
            mlines.Line2D([], [], color='orange', marker='o', markeredgecolor='black', linestyle='None', markersize=8, label='Highlighted gene set'),
        ]
        ax.legend(handles=handles, title="Significance / Highlight", loc="best")       
         # Axis and legend
        if input.show_axes():
            ax.set_xlabel(input.x_axis_label() or "Dataset 1 log2 fold change")
            ax.set_ylabel(input.y_axis_label() or "Dataset 2 log2 fold change")
            ax.set_title(input.plot_title() or "V2 Plot (Shape & Highlight)")
            # Add crosshair axes at x=0 and y=0
            ax.axhline(0, color="black", linewidth=1, linestyle="--", alpha=0.7, zorder=0)
            ax.axvline(0, color="black", linewidth=1, linestyle="--", alpha=0.7, zorder=0)
        else:
            ax.set_xlabel("Dataset 1 log2 fold change")
            ax.set_ylabel("Dataset 2 log2 fold change")
            ax.set_title("V2 Plot (Shape & Highlight)")
        return fig

    @output
    @render.text
    def pearson_r_value():
        d1, d2 = df1(), df2()
        if d1 is None or d2 is None:
            return "Upload both datasets to calculate Pearson R."
        merged = pd.merge(d1, d2, left_index=True, right_index=True, suffixes=("_1", "_2"))
        # Use log2_fold_change columns for correlation
        if "log2_fold_change_1" not in merged or "log2_fold_change_2" not in merged:
            return "log2_fold_change columns not found in both datasets."
        r, p = pearsonr(merged["log2_fold_change_1"], merged["log2_fold_change_2"])
        return f"Pearson R: {r:.5f} (p={p:.5g})"

    @output
    @render.download(filename="v2_shape_highlight_plot.png")
    def download_v2_shape_highlight_plot():
        import io
        import matplotlib.lines as mlines
        d1, d2 = df1(), df2()
        if d1 is None or d2 is None:
            return
        merged = pd.merge(d1, d2, left_index=True, right_index=True, suffixes=("_1", "_2"))
        merged["significant_1"] = merged["padj_1"] < input.padj_threshold()
        merged["significant_2"] = merged["padj_2"] < input.padj_threshold()
        shape_map = {
            (True, True): "o",
            (True, False): "s",
            (False, True): "^",
            (False, False): "D",
        }
        merged["shape"] = merged[["significant_1", "significant_2"]].apply(lambda x: shape_map[tuple(x)], axis=1)
        highlight = highlight_genes()
        merged["highlight"] = merged.index.isin(highlight)
        fig, ax = plt.subplots()
        for (sig1, sig2), marker in shape_map.items():
            mask = (merged["significant_1"] == sig1) & (merged["significant_2"] == sig2) & (~merged["highlight"])
            ax.scatter(
                merged.loc[mask, "log2_fold_change_1"],
                merged.loc[mask, "log2_fold_change_2"],
                marker=marker,
                c="gray",
                alpha=0.7,
                label=f"1:{'Yes' if sig1 else 'No'}, 2:{'Yes' if sig2 else 'No'}"
            )
        for (sig1, sig2), marker in shape_map.items():
            mask = (merged["significant_1"] == sig1) & (merged["significant_2"] == sig2) & (merged["highlight"])
            ax.scatter(
                merged.loc[mask, "log2_fold_change_1"],
                merged.loc[mask, "log2_fold_change_2"],
                marker=marker,
                c="orange",
                edgecolor="black",
                alpha=0.9,
                label=None if merged.loc[mask].empty else "Highlighted gene set"
            )
        if input.show_axes():
            ax.set_xlabel("Dataset 1 log2 fold change")
            ax.set_ylabel("Dataset 2 log2 fold change")
            ax.set_title("V2 Plot (Shape & Highlight)")
        else:
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.set_title("")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)
        x = merged["log2_fold_change_1"]
        y = merged["log2_fold_change_2"]
        min_val = min(x.min(), y.min())
        max_val = max(x.max(), y.max())
        ax.set_xlim(min_val, max_val)
        ax.set_ylim(min_val, max_val)
        ax.plot([min_val, max_val], [min_val, max_val], "k--", alpha=0.5)
        handles = [
            mlines.Line2D([], [], color='gray', marker='o', linestyle='None', markersize=8, label='1:Yes, 2:Yes'),
            mlines.Line2D([], [], color='gray', marker='s', linestyle='None', markersize=8, label='1:Yes, 2:No'),
            mlines.Line2D([], [], color='gray', marker='^', linestyle='None', markersize=8, label='1:No, 2:Yes'),
            mlines.Line2D([], [], color='gray', marker='D', linestyle='None', markersize=8, label='1:No, 2:No'),
            mlines.Line2D([], [], color='orange', marker='o', markeredgecolor='black', linestyle='None', markersize=8, label='Highlighted gene set'),
        ]
        ax.legend(handles=handles, title="Significance / Highlight", loc="best")
        buf = io.BytesIO()
        fig.savefig(buf, format="png")
        buf.seek(0)
        yield buf.read()
        
        # Axis and legend
            
        if input.show_axes():
            ax.set_xlabel(input.x_axis_label() or "Dataset 1 log2 fold change")
            ax.set_ylabel(input.y_axis_label() or "Dataset 2 log2 fold change")
            ax.set_title(input.plot_title() or "V2 Plot (Shape & Highlight)")
            # Add crosshair axes at x=0 and y=0
            ax.axhline(0, color="black", linewidth=1, linestyle="-", alpha=0.7, zorder=0)
            ax.axvline(0, color="black", linewidth=1, linestyle="-", alpha=0.7, zorder=0)
        else:
            ax.set_xlabel("Dataset 1 log2 fold change")
            ax.set_ylabel("Dataset 2 log2 fold change")
            ax.set_title("V2 Plot (Shape & Highlight)")
        return fig
    @output
    @render_widget
    def v2_shape_highlight_plotly():
        d1, d2 = df1(), df2()
        if d1 is None or d2 is None:
            return

        merged = pd.merge(d1, d2, left_index=True, right_index=True, suffixes=("_1", "_2"))
        merged["Gene"] = merged.index
        merged["significant_1"] = merged["padj_1"] < input.padj_threshold()
        merged["significant_2"] = merged["padj_2"] < input.padj_threshold()
        shape_map = {
            (True, True): "circle",
            (True, False): "square",
            (False, True): "triangle-up",
            (False, False): "diamond",
        }
        merged["shape"] = merged[["significant_1", "significant_2"]].apply(lambda x: shape_map[tuple(x)], axis=1)
        highlight = highlight_genes()
        merged["highlight"] = merged.index.isin(highlight)
        merged["color"] = merged["highlight"].map({True: "orange", False: "gray"})

        # For legend grouping
        sig_labels = {
            (True, True): "1:Yes, 2:Yes",
            (True, False): "1:Yes, 2:No",
            (False, True): "1:No, 2:Yes",
            (False, False): "1:No, 2:No",
        }
        merged["sig_label"] = merged[["significant_1", "significant_2"]].apply(lambda x: sig_labels[tuple(x)], axis=1)

        fig = go.Figure()

        # Plot non-highlighted genes by shape
        for (sig1, sig2), shape in shape_map.items():
            mask = (merged["significant_1"] == sig1) & (merged["significant_2"] == sig2) & (~merged["highlight"])
            subset = merged[mask]
            if not subset.empty:
                fig.add_trace(go.Scattergl(
                    x=subset["log2_fold_change_1"],
                    y=subset["log2_fold_change_2"],
                    mode="markers",
                    marker=dict(
                        symbol=shape,
                        color="gray",
                        size=9,
                        line=dict(width=1, color="black")
                    ),
                    name=sig_labels[(sig1, sig2)],
                    customdata=subset[["Gene", "padj_1", "padj_2"]],
                    hovertemplate=(
                        "Gene: %{customdata[0]}<br>"
                        "log2FC 1: %{x:.3f}<br>"
                        "log2FC 2: %{y:.3f}<br>"
                        "padj 1: %{customdata[1]:.3g}<br>"
                        "padj 2: %{customdata[2]:.3g}<br>"
                        f"Significance: {sig_labels[(sig1, sig2)]}<extra></extra>"
                    ),
                    showlegend=True
                ))

        # Plot highlighted genes by shape (orange)
        for (sig1, sig2), shape in shape_map.items():
            mask = (merged["significant_1"] == sig1) & (merged["significant_2"] == sig2) & (merged["highlight"])
            subset = merged[mask]
            if not subset.empty:
                fig.add_trace(go.Scattergl(
                    x=subset["log2_fold_change_1"],
                    y=subset["log2_fold_change_2"],
                    mode="markers",
                    marker=dict(
                        symbol=shape,
                        color="orange",
                        size=11,
                        line=dict(width=2, color="black")
                    ),
                    name="Highlighted gene set",
                    customdata=subset[["Gene", "padj_1", "padj_2"]],
                    hovertemplate=(
                        "Gene: %{customdata[0]}<br>"
                        "log2FC 1: %{x:.3f}<br>"
                        "log2FC 2: %{y:.3f}<br>"
                        "padj 1: %{customdata[1]:.3g}<br>"
                        "padj 2: %{customdata[2]:.3g}<br>"
                        f"Significance: {sig_labels[(sig1, sig2)]}<br>"
                        "<b>Highlighted</b><extra></extra>"
                    ),
                    showlegend=True
                ))

        # Diagonal line
        min_val = min(merged["log2_fold_change_1"].min(), merged["log2_fold_change_2"].min())
        max_val = max(merged["log2_fold_change_1"].max(), merged["log2_fold_change_2"].max())
        fig.add_shape(
            type="line",
            x0=min_val, y0=min_val, x1=max_val, y1=max_val,
            line=dict(color="black", dash="dash"),
            layer="below"
        )

        # Axes and title toggle
        if input.show_axes():
            fig.update_layout(
                xaxis_title="Dataset 1 log2 fold change",
                yaxis_title="Dataset 2 log2 fold change",
                title="V2 Plot (Plotly Shape & Highlight)",
                legend_title="Significance / Highlight",
                xaxis=dict(range=[min_val, max_val]),
                yaxis=dict(range=[min_val, max_val])
            )
        else:
            fig.update_layout(
                xaxis=dict(showticklabels=False, showgrid=False, zeroline=False, title="", range=[min_val, max_val]),
                yaxis=dict(showticklabels=False, showgrid=False, zeroline=False, title="", range=[min_val, max_val]),
                showlegend=False,
                title=""
            )

        return fig
app = App(app_ui, server)