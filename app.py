from shiny import App, ui, render, reactive, Inputs, Outputs, Session
import pandas as pd
import matplotlib.pyplot as plt
import os
import gseapy as gp

# Read the file
##mousediv10 = pd.read_excel("C:\Users\Lenovo\OneDrive\Documentos\Bioinformatics\mouse_div10.xlsx" , index_col=0)
#mousediv4 = pd.read_excel("C:\Users\Lenovo\OneDrive\Documentos\Bioinformatics\mouse_div4.xlsx" , index_col=0)
#humn = pd.read_excel("C:\Users\Lenovo\OneDrive\Documentos\Bioinformatics\human_esc.xlsx" , index_col=0)

# Pre-load MSigDB gene sets for both Human and Mouse
gene_sets_human = gp.get_library(name='GO_Biological_Process_2023', organism='Human')
gene_sets_mouse = gp.get_library(name='GO_Biological_Process_2023', organism='Mouse')

app_ui = ui.page_fluid(
    
    ui.layout_sidebar(
        ui.sidebar(
            ui.input_file("file1", "Upload Dataset 1", accept=[".csv", ".txt", ".xls", ".xlsx"]),
            ui.input_file("file2", "Upload Dataset 2", accept=[".csv", ".txt", ".xls", ".xlsx"]),
            ui.input_numeric("padj_threshold", "padj threshold", 0.05, min=0, max=1, step=0.01),
            ui.input_select(
                "organism",
                "Select organism for gene sets",
                choices=["Human", "Mouse"],
                selected="Mouse"
            ),
            ui.input_selectize(
                "gene_set_term",
                "Select a gene set",
                choices=[],  # Will be updated reactively
            ),
            ui.input_file("user_gene_set", "Upload gene set (CSV or TXT)", accept=[".csv", ".txt"]),
        ),
        ui.navset_tab(
            ui.nav_panel("V2 Plot", ui.output_plot("v2_comparison_plot")),
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
        ),
            ui.nav_panel("Gene Set Highlight", ui.output_plot("gene_set_highlight_plot")),
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
        return gene_sets_mouse if input.organism() == "Mouse" else gene_sets_human
    
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

    #v2 comparison table
    @output
    @render.plot
    def v2_comparison_plot():
        d1, d2 = df1(), df2()
        if d1 is None or d2 is None:
            return
        merged = pd.merge(d1, d2, left_index=True, right_index=True, suffixes=("_1", "_2"))
        merged["significant"] = (merged["padj_1"] < input.padj_threshold()) & (merged["padj_2"] < input.padj_threshold())

        fig, ax = plt.subplots()
        ax.scatter(merged["log2_fold_change_1"], merged["log2_fold_change_2"],
                   c=merged["significant"].map({True: "red", False: "gray"}), alpha=0.7)
        ax.set_xlabel("Dataset 1 log2 fold change")
        ax.set_ylabel("Dataset 2 log2 fold change")
        ax.set_title("Expression Comparison")
        return fig
    
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
     
    #1. Filter gene sets to only those where all genes are present in the uploaded data    @reactive.Calc
    def filtered_gene_set_dict():
        genes_in_data = set(all_genes())
        # Only include gene sets where ALL genes in the set are present in the data
        return {
            k: v for k, v in gene_set_dict().items()
            if set([g.upper() for g in v]).issubset(genes_in_data)
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
        # Only show the "Gene Set" names as choices
        gene_set_names = [str(k) for i, (k, v) in enumerate(filtered_gene_set_dict().items()) if i < 1000]
        print("Dropdown gene set names:", gene_set_names[:5])  # Debug: print first 5
        session.send_input_message("gene_set_term", {"choices": gene_set_names})

    #3. Get the selected gene set from the dropdown
    @reactive.Calc
    def selected_gene_set():
        term = input.gene_set_term()
        if not term:
            return set()
        # Always uppercase for matching
        return set([g.upper() for g in filtered_gene_set_dict().get(term, [])])
    
    @output
    @render.plot
    def gene_set_highlight_plot():
        d1, d2 = df1(), df2()
        if d1 is None or d2 is None:
            return
        merged = pd.merge(d1, d2, left_index=True, right_index=True, suffixes=("_1", "_2"))
        highlight_genes = selected_gene_set()
        merged["highlight"] = merged.index.isin(highlight_genes)
        fig, ax = plt.subplots()
        ax.scatter(
            merged["log2_fold_change_1"],
            merged["log2_fold_change_2"],
            c=merged["highlight"].map({True: "orange", False: "gray"}),
            alpha=0.7
        )
        ax.set_xlabel("Dataset 1 log2 fold change")
        ax.set_ylabel("Dataset 2 log2 fold change")
        ax.set_title("Highlight: Selected Gene Set")
        return fig

app = App(app_ui, server)