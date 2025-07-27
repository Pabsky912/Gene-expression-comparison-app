from shiny import App, ui, render, reactive, Inputs, Outputs, Session
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import gseapy as gp
from gseapy import Msigdb
from scipy.stats import pearsonr
from shinywidgets import render_widget, output_widget
import plotly.graph_objects as go
import matplotlib.lines as mlines
import io

all_gseapy_libraries = gp.get_library_name()

# If different organisms, use ortholog mapping
ortho_path = "ensembl99_human_mouse_orthologs (1).csv"
ortho_df = pd.read_csv(ortho_path)
ortho_df["Gene name"] = ortho_df["Gene name"].str.upper()
ortho_df["Human gene name"] = ortho_df["Human gene name"].str.upper()

# Read the file
#mousediv10 = pd.read_excel("C:\Users\Lenovo\OneDrive\Documentos\Bioinformatics\mouse_div10.xlsx" , index_col=0)
#mousediv4 = pd.read_excel("C:\Users\Lenovo\OneDrive\Documentos\Bioinformatics\mouse_div4.xlsx" , index_col=0)
#humn = pd.read_excel("C:\Users\Lenovo\OneDrive\Documentos\Bioinformatics\human_esc.xlsx" , index_col=0)

app_ui = ui.page_fluid(
    ui.layout_sidebar(
        ui.sidebar(
            ui.h5("Graph Controls"),
            ui.input_numeric("point_size", "Point Size", value=9, min=1, max=100, step=1),
            ui.input_text("plot_title", "Plot Title", value="V2 Plot (Shape & Highlight)"),
            ui.input_text("x_axis_label", "X Axis Label", value="Dataset 1 log2 fold change"),
            ui.input_text("y_axis_label", "Y Axis Label", value="Dataset 2 log2 fold change"),
            ui.input_numeric("padj_threshold", "padj threshold", 0.05, min=0, max=1, step=0.01),
            ui.input_checkbox("show_axes", "Show axes on V2 Plot", value=True),
            ui.input_checkbox("show_diagonal", "Show diagonal line", value=True),
        ),
        ui.navset_tab(
            ui.nav_panel(
                "Data",
                output_widget("v2_shape_highlight_plotly"),
                ui.card(ui.output_text("pearson_r_value"),ui.output_text("pearson_r_selected_gene_set")), ui.output_text("pearson_r_enriched_genes"),
                    # Second: radio buttons row (full width)
                ui.hr(),    
                ui.input_radio_buttons(
                    "gene_set_source",
                    ui.h3("Choose gene set source for highlight:"),
                    choices=["Upload", "Database", "Enrichment", "Download"],
                    inline=True
                ),
                ui.panel_conditional(
                    "input.gene_set_source === 'Download'",
                    ui.markdown("**Preview of download (The size will autoadjust)**"),
                    ui.output_plot("v2_shape_highlight_plot"),
                    ui.card(ui.download_button("download_v2_shape_highlight_plot", "Download PNG"))
                ),
                # Panel conditionals for radio buttons
                ui.panel_conditional(
                    "input.gene_set_source === 'Upload'",
                    ui.input_file("user_gene_set", "Upload gene set", accept=[".csv", ".txt", ".xls", ".xlsx"]),
                    ui.input_selectize("search_gene", "Search and highlight gene on plot:", choices=[], multiple=True),
                ),
                ui.panel_conditional(
                    "input.gene_set_source === 'Database'",
                    ui.input_select(
                        "db_library",
                        "Choose Database Library",
                        choices=["Enrichr", "MsigDB"]
                    ),
                    ui.panel_conditional(
                        "input.db_library === 'Enrichr'",
                        ui.input_select(
                            "gene_set_library",
                            "Select gene set library",
                            choices=all_gseapy_libraries,
                            selected="GO_Biological_Process_2023" if "GO_Biological_Process_2023" in all_gseapy_libraries else all_gseapy_libraries[0]
                        ),
                        ui.input_select(
                            "organism",
                            "Select organism for gene set libraries",
                            choices=["Human", "Mouse", "Yeast", "Fly", "Fish", "Worm"],
                            selected="Mouse"
                        ),
                        ui.input_selectize(
                            "gene_set_term",
                            "Select a gene set",
                            choices=[],
                            multiple=True,
                        ),
                    ),
                    
                    ui.panel_conditional(
                        "input.db_library === 'MsigDB'",
                        ui.input_select(
                            "msigdb_version",
                            "Select MsigDB version",
                            choices=[], 
                            selected='2025.1.Mm'
                        ),
                        ui.input_select(
                            "msigdb_library",
                            "Select MsigDB library",
                            choices=[],  # Will be updated reactively
                            selected=None
                        ),
                        ui.input_selectize(
                            "msigdb_gene_set_term",
                            "Select a gene set",
                            choices=[],
                            multiple=True,
                        ),
                    )
                ),
                ui.panel_conditional(
                    "input.gene_set_source === 'Enrichment'",
                    ui.input_radio_buttons(
                        "enrich_significance",
                        "Which genes to use for enrichment?",
                        choices={
                            "d1": "Significant in Dataset 1",
                            "d2": "Significant in Dataset 2",
                            "both": "Significant in Both",
                        },
                        selected=None,
                        inline=True
                    ),
                    ui.input_select(
                        "enrich_gene_set_library",
                        "Gene set library for enrichment:",
                        choices=all_gseapy_libraries,
                        selected="GO_Biological_Process_2023" if "GO_Biological_Process_2023" in all_gseapy_libraries else all_gseapy_libraries[0]
                    ),
                    ui.input_action_button("run_enrichr", "Run Enrichment"),
                    ui.input_selectize("enriched_gene_set", "Select enriched gene set", choices=[], multiple=True),
                    ui.markdown("**Enrichment Results:**"),
                    ui.output_table("enrichr_results"),
                ),
            ),
            ui.nav_panel(
                "Upload",
                ui.layout_columns(
                    ui.column(
                        7,
                        ui.h5("Dataset 1"),
                        ui.input_file("file1", "Upload Dataset 1", accept=[".csv", ".txt", ".xls", ".xlsx"]),
                        ui.input_select(
                            "organism1",
                            "Select organism for gene sets",
                            choices=["Human", "Mouse"],
                            selected="Mouse"
                        )
                    ),
                    ui.column(
                        7,
                        ui.h5("Dataset 2"),
                        ui.input_file("file2", "Upload Dataset 2", accept=[".csv", ".txt", ".xls", ".xlsx"]),
                        ui.input_select(
                            "organism2",
                            "Select organism for gene sets",
                            choices=["Human", "Mouse"],
                            selected="Mouse"
                        ),
                    )
                ),
                ui.hr(),
                ui.input_checkbox_group(
                    "join_method",
                    "Join datasets by:",
                    choices={
                        "Gene Stable ID": "Gene Stable ID",
                        "Gene Name": "Gene Name"
                    },
                    selected= "Gene Name",
                    inline=True
                ),
                ui.markdown("**Merged Dataset Preview:**"),
                ui.output_table("merged_preview"),
            ),
            ui.nav_panel(
                "About & Gene Lists for potential errors",
                ui.markdown("*This is a panel for displaying information about the datasets and potential issues with gene lists.*"),
                ui.markdown("**Preview dataset 1:**"),
                ui.output_table("preview1"),
                ui.markdown("**Preview dataset 2:**"),
                ui.output_table("preview2"),
                ui.markdown("**All genes in both datasets list for gene set highlighting**"),
                ui.output_table("head_all_genes"),
                ui.markdown("**Enrichr Mouse Gene Sets (head):**"),
                ui.output_table("head_mouse_gene_sets"),
                ui.markdown("**Enrichr Human Gene Sets (head):**"),
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
                df = pd.read_csv(path, sep=None, engine="python")
            elif ext in [".xls", ".xlsx"]:
                df = pd.read_excel(path)
            else:
                return None
            if 'gene_name' in df.columns:
                df['gene_name'] = df['gene_name'].astype(str).str.upper()
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
            d1 = d1.copy()
            genes.update(pd.Series(d1['gene_name']).dropna().astype(str).str.upper().unique())
        if d2 is not None and d2.shape[0] > 0:
            d2 = d2.copy()
            genes.update(pd.Series(d2['gene_name']).dropna().astype(str).str.upper().unique())
        return genes
    
    @reactive.effect
    def update_search_gene_dropdown():
        merged_genes = merged_species_df().get('gene_name', pd.Series([])).dropna().astype(str).str.upper().unique()
        gene_choices = list(merged_genes)
        ui.update_selectize(
            "search_gene",
            choices=gene_choices[:1000],  # Show only the first 1000 initially
            session=session,
        )
    
    @reactive.Calc
    def search_gene_list():
        search_gene = input.search_gene()
        if isinstance(search_gene, (list, tuple)):
            search_gene_set = set(g.upper() for g in search_gene if g)
        else:
            search_gene_set = {search_gene.upper()} if search_gene else set()
        return search_gene_set
        
    @reactive.Calc
    def msigdb_versions():
        try:
            msig = Msigdb()
            df = msig.list_dbver()
            dates = df['Name'].tolist()
            return dates
        except Exception as e:
            print(f"Msigdb version error: {e}")
            return []

    @reactive.effect
    def update_msigdb_version_dropdown():
        ui.update_select(
            "msigdb_version",
            choices=msigdb_versions(),
            session=session,
        )

    @reactive.Calc
    def msigdb_libraries():
        try:
            msig = Msigdb()
            df = msig.list_category()
            return df
        except Exception as e:
            print(f"Msigdb library error: {e}")
            return []
    
    @reactive.effect
    def update_msigdb_library_dropdown():
        ui.update_select(
            "msigdb_library",
            choices=msigdb_libraries(),
            session=session,
        )

    @reactive.Calc
    def msigdb_gene_sets():
        # Get gene sets for selected category and version
        version = input.msigdb_version() 
        library = input.msigdb_library()
        try:
            if version and library:
                msig = Msigdb()
                return msig.get_gmt(category=library, dbver=version)
            else:
                return {}
        except Exception as e:
            print(f"Msigdb gene set error: {e}")
            return {}

    @reactive.Calc
    def filtered_msigdb_gene_sets():
        # Filter gene sets to those with at least one gene in uploaded data
        genes_in_data = set(all_genes())
        gene_sets = msigdb_gene_sets()
        if not isinstance(gene_sets, dict) or gene_sets is None:
            gene_sets = {}
        return {
            k: v for k, v in gene_sets.items()
            if any(g.upper() in genes_in_data for g in v)
        }
    
    @reactive.Effect
    def update_msigdb_gene_set_term_dropdown():
        terms = list(filtered_msigdb_gene_sets().keys())
        ui.update_selectize(
            "msigdb_gene_set_term",
            choices=terms,
            session=session,
        )

    @reactive.Calc
    def selected_msigdb_gene_set():
        terms = input.msigdb_gene_set_term()
        if not terms:
            return set()
        if isinstance(terms, str):
            terms = [terms]
        genes = set()
        for term in terms:
            genes.update([g.upper() for g in filtered_msigdb_gene_sets().get(term, [])])
        return genes

    @reactive.Calc
    def merged_species_df():
        def get_col(df, *names, default=None):
            """Try to get the first column that exists in df from names, else return default or raise KeyError."""
            for n in names:
                if n in df.columns:
                    return df[n]
            if default is not None:
                return default
            raise KeyError(f"None of columns {names} found in DataFrame columns: {df.columns}")

        d1, d2 = df1(), df2()
        if d1 is None or d2 is None:
            return pd.DataFrame({"Message": ["No merged data available."]})
        org1 = input.organism1() if "organism1" in input else "Mouse"
        org2 = input.organism2() if "organism2" in input else "Mouse"
        join_methods = input.join_method() if "join_method" in input else ["Gene Name"]

        if "Gene Name" in join_methods and not "Gene Stable ID" in join_methods:
            merged_name = pd.merge(d1, d2, left_on = 'gene_name', right_on = 'gene_name', suffixes=("_1", "_2"))
            return merged_name
        if "Gene Stable ID" in join_methods and not "Gene Name" in join_methods:
            if org1 == org2:
                if "gene" in d1.columns and "gene" in d2.columns:
                    merged_id = pd.merge(d1, d2, on="gene", suffixes=("_1", "_2"))
                    print(merged_id)
                return merged_id
            elif org1 == "Mouse" and org2 == "Human":
                # Map mouse gene to human gene using ortho_df
                if "gene" in d1.columns:
                    d1_id = d1.copy()
                    merged_ortho = pd.merge(d1_id, ortho_df, left_on="gene", right_on="Gene stable ID", how="inner")
                    merged_id = merged_ortho.merge(d2, left_on="Human gene stable ID", right_on="gene", suffixes=("_1", "_2"))
                else:
                    print("Gene Not Found In Columns")
                return merged_id
            elif org1 == "Human" and org2 == "Mouse":
                # Map human gene to mouse gene using ortho_df
                if "gene" in d2.columns:
                    d1_id = d1.copy()
                    merged_ortho = pd.merge(d1_id, ortho_df, left_on="gene", right_on="Human gene stable ID", how="inner")
                    merged_id = merged_ortho.merge(d2, left_on="Gene stable ID", right_on="gene", suffixes=("_1", "_2"))
                else:
                    print("Gene Not Found In Columns")
                return merged_id
        if "Gene Stable ID" in join_methods and "Gene Name" in join_methods:
            # Mouse-Mouse
            if org1 == "Mouse" and org2 == "Mouse":
                if "gene" in d1.columns and "gene" in d2.columns:
                    merged_name = pd.merge(d1, d2, left_on='gene_name', right_on='gene_name', suffixes=("_1", "_2"))
                    merged_id = pd.merge(d1, d2, on="gene", suffixes=("_1", "_2"))

                    merged_name_clean = pd.DataFrame({
                        'gene_name': get_col(merged_name, 'gene_name', 'gene_name_1', 'gene_name_2'),
                        'gene_mouse1': get_col(merged_name, 'gene_1', 'gene'),
                        'gene_mouse2': get_col(merged_name, 'gene_2', 'gene'),
                        'log2_fold_change_1': get_col(merged_name, 'log2_fold_change_1'),
                        'log2_fold_change_2': get_col(merged_name, 'log2_fold_change_2'),
                        'padj_1': get_col(merged_name, 'padj_1'),
                        'padj_2': get_col(merged_name, 'padj_2')
                    })

                    merged_id_clean = pd.DataFrame({
                        'gene_name': get_col(merged_id, 'gene_name', 'gene_name_1', 'gene_name_2'),
                        'gene_mouse1': get_col(merged_id, 'gene_1', 'gene'),
                        'gene_mouse2': get_col(merged_id, 'gene_2', 'gene'),
                        'log2_fold_change_1': get_col(merged_id, 'log2_fold_change_1'),
                        'log2_fold_change_2': get_col(merged_id, 'log2_fold_change_2'),
                        'padj_1': get_col(merged_id, 'padj_1'),
                        'padj_2': get_col(merged_id, 'padj_2')
                    })

                    merged_nameid = pd.concat([merged_id_clean, merged_name_clean], axis=0, ignore_index=True)
                    merged_nameid = merged_nameid.drop_duplicates(keep='first')
                    return merged_nameid
                else:
                    return pd.DataFrame({"Message": ["Gene column not found."]})

            # Human-Human
            if org1 == "Human" and org2 == "Human":
                if "gene" in d1.columns and "gene" in d2.columns:
                    merged_name = pd.merge(d1, d2, left_on='gene_name', right_on='gene_name', suffixes=("_1", "_2"))
                    merged_id = pd.merge(d1, d2, on="gene", suffixes=("_1", "_2"))

                    merged_name_clean = pd.DataFrame({
                        'gene_name': get_col(merged_name, 'gene_name', 'gene_name_1', 'gene_name_2'),
                        'gene_human1': get_col(merged_name, 'gene_1', 'gene'),
                        'gene_human2': get_col(merged_name, 'gene_2', 'gene'),
                        'log2_fold_change_1': get_col(merged_name, 'log2_fold_change_1'),
                        'log2_fold_change_2': get_col(merged_name, 'log2_fold_change_2'),
                        'padj_1': get_col(merged_name, 'padj_1'),
                        'padj_2': get_col(merged_name, 'padj_2')
                    })

                    merged_id_clean = pd.DataFrame({
                        'gene_name': get_col(merged_id, 'gene_name', 'gene_name_1', 'gene_name_2'),
                        'gene_human1': get_col(merged_id, 'gene_1', 'gene'),
                        'gene_human2': get_col(merged_id, 'gene_2', 'gene'),
                        'log2_fold_change_1': get_col(merged_id, 'log2_fold_change_1'),
                        'log2_fold_change_2': get_col(merged_id, 'log2_fold_change_2'),
                        'padj_1': get_col(merged_id, 'padj_1'),
                        'padj_2': get_col(merged_id, 'padj_2')
                    })

                    merged_nameid = pd.concat([merged_id_clean, merged_name_clean], axis=0, ignore_index=True)
                    merged_nameid = merged_nameid.drop_duplicates(keep='first')
                    return merged_nameid
                else:
                    return pd.DataFrame({"Message": ["Gene column not found."]})

            # Mouse-Human
            if org1 == "Mouse" and org2 == "Human":
                if "gene" in d1.columns:
                    d1_id = d1.copy()
                    merged_name = pd.merge(d1, d2, left_on='gene_name', right_on='gene_name', suffixes=("_1", "_2"))
                    merged_ortho = pd.merge(d1_id, ortho_df, left_on="gene", right_on="Gene stable ID", how="inner")
                    merged_id = pd.merge(merged_ortho, d2, left_on="Human gene stable ID", right_on="gene", suffixes=("_1", "_2"))

                    merged_name_clean = pd.DataFrame({
                        'gene_name': get_col(merged_name, 'gene_name', 'gene_name_1', 'gene_name_2'),
                        'gene_mouse': get_col(merged_name, 'gene_1', 'gene'),
                        'gene_human': get_col(merged_name, 'gene_2', 'gene'),
                        'log2_fold_change_1': get_col(merged_name, 'log2_fold_change_1'),
                        'log2_fold_change_2': get_col(merged_name, 'log2_fold_change_2'),
                        'padj_1': get_col(merged_name, 'padj_1'),
                        'padj_2': get_col(merged_name, 'padj_2')
                    })

                    merged_id_clean = pd.DataFrame({
                        'gene_name': get_col(merged_id, 'gene_name', 'gene_name_1', 'gene_name_2'),
                        'gene_mouse': get_col(merged_id, 'gene_1', 'gene'),
                        'gene_human': get_col(merged_id, 'gene_2', 'gene'),
                        'log2_fold_change_1': get_col(merged_id, 'log2_fold_change_1'),
                        'log2_fold_change_2': get_col(merged_id, 'log2_fold_change_2'),
                        'padj_1': get_col(merged_id, 'padj_1'),
                        'padj_2': get_col(merged_id, 'padj_2')
                    })

                    merged_nameid = pd.concat([merged_id_clean, merged_name_clean], axis=0, ignore_index=True)
                    merged_nameid = merged_nameid.drop_duplicates(keep='first')
                    return merged_nameid
                else:
                    return pd.DataFrame({"Message": ["Gene column not found."]})

            # Human-Mouse
            if org1 == "Human" and org2 == "Mouse":
                if "gene" in d1.columns:
                    d1_id = d1.copy()
                    merged_name = pd.merge(d1, d2, left_on='gene_name', right_on='gene_name', suffixes=("_1", "_2"))
                    merged_ortho = pd.merge(d1_id, ortho_df, left_on="gene", right_on="Human gene stable ID", how="inner")
                    merged_id = pd.merge(merged_ortho, d2, left_on="Gene stable ID", right_on="gene", suffixes=("_1", "_2"))

                    merged_name_clean = pd.DataFrame({
                        'gene_name': get_col(merged_name, 'gene_name', 'gene_name_1', 'gene_name_2'),
                        'gene_human': get_col(merged_name, 'gene_1', 'gene'),
                        'gene_mouse': get_col(merged_name, 'gene_2', 'gene'),
                        'log2_fold_change_1': get_col(merged_name, 'log2_fold_change_1'),
                        'log2_fold_change_2': get_col(merged_name, 'log2_fold_change_2'),
                        'padj_1': get_col(merged_name, 'padj_1'),
                        'padj_2': get_col(merged_name, 'padj_2')
                    })

                    merged_id_clean = pd.DataFrame({
                        'gene_name': get_col(merged_id, 'gene_name', 'gene_name_1', 'gene_name_2'),
                        'gene_human': get_col(merged_id, 'gene_1', 'gene'),
                        'gene_mouse': get_col(merged_id, 'gene_2', 'gene'),
                        'log2_fold_change_1': get_col(merged_id, 'log2_fold_change_1'),
                        'log2_fold_change_2': get_col(merged_id, 'log2_fold_change_2'),
                        'padj_1': get_col(merged_id, 'padj_1'),
                        'padj_2': get_col(merged_id, 'padj_2')
                    })

                    merged_nameid = pd.concat([merged_id_clean, merged_name_clean], axis=0, ignore_index=True)
                    merged_nameid = merged_nameid.drop_duplicates(keep='first')
                    return merged_nameid
                else:
                    return pd.DataFrame({"Message": ["Gene column not found."]})
    @output
    @render.table
    def merged_preview():
        df = merged_species_df()
        if df is None or df.empty:
            return pd.DataFrame({"Message": ["No merged data available."]})
        # Show first 20 rows for preview
        return df.tail(20)

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
        if df is not None:
            return df.head(10).reset_index()  # Show index as a column
        return None

    @output
    @render.table
    def preview2():
        df = df2()
        if df is not None:
            return df.head(10).reset_index()  # Show index as a column
        return None
    
    #tab for all gene lists - all genes, gpseea mouse and human gene sets
    @output
    @render.table
    def head_all_genes():
        genes = list(all_genes())
        return pd.DataFrame(genes[:10], columns=["Gene"])

    @output
    @render.table
    def head_mouse_gene_sets():
        # Show the first 5 gene sets and their first 5 genes from the current gene_set_dict (Mouse)
        data = []
        mouse_sets = {k: v for k, v in gene_set_dict().items() if "mouse" in input.organism().lower()}
        for i, (k, v) in enumerate(mouse_sets.items()):
            if i >= 5:
                break
            data.append({"Gene Set": k, "Genes": ", ".join(v[:5]) + ("..." if len(v) > 5 else "")})
        return pd.DataFrame(data)

    @output
    @render.table
    def head_human_gene_sets():
        # Show the first 5 gene sets and their first 5 genes from the current gene_set_dict (Human)
        data = []
        human_sets = {k: v for k, v in gene_set_dict().items() if "human" in input.organism().lower()}
        for i, (k, v) in enumerate(human_sets.items()):
            if i >= 5:
                break
            data.append({"Gene Set": k, "Genes": ", ".join(v[:5]) + ("..." if len(v) > 5 else "")})
        return pd.DataFrame(data)

     
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
    
    @output
    @render.table
    def user_uploaded_gene_names():
        genes = list(user_uploaded_genes())
        if not genes:
            return pd.DataFrame({"Gene": ["(No genes uploaded)"]})
        return pd.DataFrame(genes, columns=["Gene"])

    # Read user-uploaded gene set (first row as header, first column as gene names)
    @reactive.Calc
    def user_uploaded_genes():
        try:
            df = read_any_file(input.user_gene_set())
            if df is None or df.empty:
                return set()
            # Try to find a likely gene name column
            possible_cols = [col for col in df.columns if str(col).lower() in ["gene_name", "gene", "gene_id", "symbol"]]
            if possible_cols:
                genes = pd.Series(df[possible_cols[0]]).dropna().astype(str).str.upper().unique()
            else:
                # Use the first column if no known gene name column is found
                genes = pd.Series(df[df.columns[0]]).dropna().astype(str).str.upper().unique()
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
        merged = pd.merge(d1, d2, left_on = 'gene_name', right_on = 'gene_name', suffixes=("_1", "_2"))
        highlight = highlight_genes()
        # Only keep genes in the selected gene set
        filtered = merged.loc[merged["gene_name"].isin(highlight)].copy()
        if filtered.shape[0] < 2:
            return "Not enough significant genes in the selected gene set."
        if "log2_fold_change_1" not in filtered or "log2_fold_change_2" not in filtered:
            return "log2_fold_change columns not found in both datasets."
        r, p = pearsonr(filtered["log2_fold_change_1"], filtered["log2_fold_change_2"])
        return f"Pearson R (selected gene set): {r:.5f} (p={p:.5g})"
        
    @output
    @render.text
    def pearson_r_enriched_genes():
        d1, d2 = df1(), df2()
        if d1 is None or d2 is None:
            return "Upload both datasets to calculate Pearson R."
        merged = pd.merge(d1, d2, left_on='gene_name', right_on='gene_name', suffixes=("_1", "_2"))
        enriched_genes = genes_fromenrichedgeneset()
        filtered = merged.loc[merged["gene_name"].isin(enriched_genes)].copy()
        if filtered.shape[0] < 2:
            return "Not enough enriched genes for correlation."
        if "log2_fold_change_1" not in filtered or "log2_fold_change_2" not in filtered:
            return "log2_fold_change columns not found in both datasets."
        r, p = pearsonr(filtered["log2_fold_change_1"], filtered["log2_fold_change_2"])
        return f"Pearson R (enriched genes): {r:.5f} (p={p:.5g})"

    @reactive.Calc
    def significant_genes():
        d1, d2 = df1(), df2()
        if d1 is None or d2 is None:
            return pd.DataFrame({"Gene": [], "Source": []})

        sig1 = set(d1[d1["padj"] < input.padj_threshold()]["gene_name"]) if "padj" in d1 else set()
        sig2 = set(d2[d2["padj"] < input.padj_threshold()]["gene_name"]) if "padj" in d2 else set()

        which = input.enrich_significance() if "enrich_significance" in input else "both"

        if which == "d1":
            return sig1
        elif which == "d2":
            return sig2
        else:  # both
            return sig1 & sig2

    @reactive.Calc
    def enrichr_results_df():
        if input.run_enrichr() == 0:
            return pd.DataFrame({"Message": ["Click 'Run Enrichr Analysis' to perform enrichment."]})
        genes = list(significant_genes())
        if not genes:
            return pd.DataFrame({"Message": ["No significant genes found."]})
        lib = input.enrich_gene_set_library()
        org = input.organism().lower()
        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=lib,
                organism=org,
                outdir=None
            )
            enriched_results = enr.results
            return enriched_results
        except Exception as e:
            return pd.DataFrame({"Error": [str(e)]})

    @output
    @render.table
    def enrichr_results():
        return enrichr_results_df().head(20)
    
    @reactive.Effect
    def enriched_gene_sets():
        df = enrichr_results_df()
        if df is not None and "Term" in df.columns:
            choices = list(df["Term"].dropna().str.strip().str.upper().unique())[:1000]
        else:
            choices = []
        ui.update_selectize(
            "enriched_gene_set",
            choices=choices,
            session=session,
        )

    def genes_fromenrichedgeneset():
        selected_terms = input.enriched_gene_set()
        if not selected_terms:
            return set()
        df = enrichr_results_df()
        if df is None or df.empty or "Genes" not in df.columns:
            return set()
        # Get genes for the selected terms
        genes = set()
        for term in selected_terms:
            term_genes = df[df["Term"].str.strip().str.upper() == term]["Genes"]
            if not term_genes.empty:
                genes.update(term_genes.iloc[0].split(";"))
        return set(g.strip().upper() for g in genes)
    
    @reactive.Calc
    def enriched_term_genes():
        df = enrichr_results_df()
        if df is not None and not df.empty and "Genes" in df.columns:
            all_genes = set()
            for genes_str in df["Genes"]:
                genes = genes_str.split(";")
                all_genes.update(g.strip().upper() for g in genes)
            return all_genes
        return set()
    
    # Combine MsigDB selection with other highlight sources if needed:
    @reactive.Calc
    def highlight_genes():
        # Combine Enrichr, MsigDB, and uploaded genes
        enrichr_genes = selected_gene_set()  # from Enrichr
        msigdb_genes = selected_msigdb_gene_set()
        uploaded_genes = user_uploaded_genes()
        return enrichr_genes.union(msigdb_genes).union(uploaded_genes)

    @output
    @render.plot
    def v2_shape_highlight_plot():
        d1, d2 = df1(), df2()
        if d1 is None or d2 is None:
            return
        merged = pd.merge(d1, d2, left_on = 'gene_name', right_on = 'gene_name', suffixes=("_1", "_2"))
        merged["significant_1"] = merged["padj_1"] < input.padj_threshold()
        merged["significant_2"] = merged["padj_2"] < input.padj_threshold()
        merged["gene_name"] = merged["gene_name"].astype(str).str.upper()
        shape_map = {
            (True, True): "o",      # Circle
            (True, False): "s",     # Square
            (False, True): "^",     # Triangle
            (False, False): "D",    # Diamond
        }
        merged["shape"] = merged[["significant_1", "significant_2"]].apply(
            lambda x: shape_map.get(tuple(x), "o"), axis=1
        )
        highlight = highlight_genes()
        red_highlight = genes_fromenrichedgeneset()
        merged["highlight"] = merged["gene_name"].isin(highlight)
        merged["red_highlight"] = merged["gene_name"].isin(red_highlight)
        fig, ax = plt.subplots()
        point_size = input.point_size() if input.point_size() else 9
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
        red_mask = merged["red_highlight"]
        if red_mask.any():
            ax.scatter(
                merged.loc[red_mask, "log2_fold_change_1"],
                merged.loc[red_mask, "log2_fold_change_2"],
                marker=marker,
                c="red",
                edgecolor="black",
                alpha=0.9,
                s=point_size + 10,
                label="Enriched gene set (red)"
            )
        search_gene = input.search_gene()
        if isinstance(search_gene, (list, tuple)):
            search_gene_set = set(g.upper() for g in search_gene if g)
        else:
            search_gene_set = {search_gene.upper()} if search_gene else set()
        merged["search_highlight"] = merged["gene_name"].str.upper().isin(search_gene_set)
        if search_gene:
            for (sig1, sig2), shape in shape_map.items():
                mask = (merged["significant_1"] == sig1) & (merged["significant_2"] == sig2) & (merged["search_highlight"])
                subset = merged[mask]
                if not subset.empty:
                    fig.add_trace(go.Scattergl(
                        x=subset["log2_fold_change_1"],
                        y=subset["log2_fold_change_2"],
                        mode="markers+text",
                        marker=dict(
                            symbol=shape,
                            color="blue",
                            size=point_size,
                            line=dict(width=3, color="black")
                        ),
                        name="Searched gene",
                        customdata=subset[["gene_name", "padj_1", "padj_2"]],
                        hovertemplate=...,
                        text=subset["gene_name"],
                        textposition="top center",
                        showlegend=True
                    ))
        x = merged["log2_fold_change_1"]
        y = merged["log2_fold_change_2"]
        min_val = min(x.min(), y.min())
        max_val = max(x.max(), y.max())
        ax.set_xlim(min_val, max_val)
        ax.set_ylim(min_val, max_val)
        if input.show_diagonal():
            ax.plot([min_val, max_val], [min_val, max_val], "k--", alpha=0.5)
        handles = [
            mlines.Line2D([], [], color='gray', marker='o', linestyle='None', markersize=8, label='1:Yes, 2:Yes'),
            mlines.Line2D([], [], color='gray', marker='s', linestyle='None', markersize=8, label='1:Yes, 2:No'),
            mlines.Line2D([], [], color='gray', marker='^', linestyle='None', markersize=8, label='1:No, 2:Yes'),
            mlines.Line2D([], [], color='gray', marker='D', linestyle='None', markersize=8, label='1:No, 2:No'),
        ]
        if merged["highlight"].any():
            handles.append(
                mlines.Line2D([], [], color='orange', marker='o', markeredgecolor='black', linestyle='None', markersize=8, label='Highlighted gene set')
            )
        if merged["red_highlight"].any():
            handles.append(
                mlines.Line2D([], [], color='red', marker='o', markeredgecolor='black', linestyle='None', markersize=8, label='Enriched gene set (red)')
            )
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
        merged = pd.merge(d1, d2, left_on = 'gene_name', right_on = 'gene_name', suffixes=("_1", "_2"))
        # Use log2_fold_change columns for correlation
        if "log2_fold_change_1" not in merged or "log2_fold_change_2" not in merged:
            return "log2_fold_change columns not found in both datasets."
        r, p = pearsonr(merged["log2_fold_change_1"], merged["log2_fold_change_2"])
        return f"Pearson R: {r:.5f} (p={p:.5g})"

    @output
    @render.download(filename="v2_shape_highlight_plot.png")
    def download_v2_shape_highlight_plot():

        d1, d2 = df1(), df2()
        if d1 is None or d2 is None:
            return
        merged = pd.merge(d1, d2, left_on = 'gene_name', right_on = 'gene_name', suffixes=("_1", "_2"))
        merged["gene_name"] = merged["gene_name"].astype(str).str.upper()
        merged["significant_1"] = merged["padj_1"] < input.padj_threshold()
        merged["significant_2"] = merged["padj_2"] < input.padj_threshold()
        shape_map = {
            (True, True): "o",      # Circle
            (True, False): "s",     # Square
            (False, True): "^",     # Triangle
            (False, False): "D",    # Diamond
        }
        merged["shape"] = merged[["significant_1", "significant_2"]].apply(
            lambda x: shape_map.get(tuple(x), "o"), axis=1
        )
        highlight = highlight_genes()
        red_highlight = genes_fromenrichedgeneset()
        merged["highlight"] = merged['gene_name'].isin(highlight)
        merged["red_highlight"] = merged['gene_name'].isin(red_highlight)
        fig, ax = plt.subplots()
        point_size = input.point_size() if input.point_size() else 9
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
        red_mask = merged["red_highlight"]
        if red_mask.any():
            ax.scatter(
                merged.loc[red_mask, "log2_fold_change_1"],
                merged.loc[red_mask, "log2_fold_change_2"],
                marker=marker,
                c="red",
                edgecolor="black",
                alpha=1.0,
                s=point_size + 20,
                label="Enriched gene set (red)"
            )
        search_gene = input.search_gene()
        if isinstance(search_gene, (list, tuple)):
            search_gene_set = set(g.upper() for g in search_gene if g)
        else:
            search_gene_set = {search_gene.upper()} if search_gene else set()
        merged["search_highlight"] = merged["gene_name"].str.upper().isin(search_gene_set)
        if search_gene_set:
            for (sig1, sig2), marker in shape_map.items():
                shape_mask = (
                    (merged["significant_1"] == sig1)
                    & (merged["significant_2"] == sig2)
                    & (merged["gene_name"].str.upper().isin(search_gene_set))
                )
                if shape_mask.any():
                    ax.scatter(
                        merged.loc[shape_mask, "log2_fold_change_1"],
                        merged.loc[shape_mask, "log2_fold_change_2"],
                        marker=marker,
                        c="blue",
                        edgecolor="black",
                        alpha=1.0,
                        s=point_size + 20,
                        label="Searched gene"
                    )
        x = merged["log2_fold_change_1"]
        y = merged["log2_fold_change_2"]
        min_val = min(x.min(), y.min())
        max_val = max(x.max(), y.max())
        ax.set_xlim(min_val, max_val)
        ax.set_ylim(min_val, max_val)
        if input.show_diagonal():
            ax.plot([min_val, max_val], [min_val, max_val], "k--", alpha=0.5)
        handles = [
            mlines.Line2D([], [], color='gray', marker='o', linestyle='None', markersize=8, label='1:Yes, 2:Yes'),
            mlines.Line2D([], [], color='gray', marker='s', linestyle='None', markersize=8, label='1:Yes, 2:No'),
            mlines.Line2D([], [], color='gray', marker='^', linestyle='None', markersize=8, label='1:No, 2:Yes'),
            mlines.Line2D([], [], color='gray', marker='D', linestyle='None', markersize=8, label='1:No, 2:No'),
        ]
        if merged["highlight"].any():
            handles.append(
                mlines.Line2D([], [], color='orange', marker='o', markeredgecolor='black', linestyle='None', markersize=8, label='Highlighted gene set')
            )
        if merged["red_highlight"].any():
            handles.append(
                mlines.Line2D([], [], color='red', marker='o', markeredgecolor='black', linestyle='None', markersize=8, label='Enriched gene set (red)')
            )
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
        buf = io.BytesIO()
        fig.savefig(buf, format="png")
        buf.seek(0)
        yield buf.read()

    @output
    @render_widget
    def v2_shape_highlight_plotly():
        d1, d2 = df1(), df2()
        if d1 is None or d2 is None:
            return

        merged = merged_species_df()
        merged["significant_1"] = merged["padj_1"] < input.padj_threshold()
        merged["significant_2"] = merged["padj_2"] < input.padj_threshold()
        merged["gene_name"] = merged["gene_name"].astype(str).str.upper()
        shape_map = {
            (True, True): "circle",
            (True, False): "square",
            (False, True): "triangle-up",
            (False, False): "diamond",
        }
        merged["shape"] = merged[["significant_1", "significant_2"]].apply(
            lambda x: shape_map.get(tuple(x), "o"), axis=1
        )
        highlight = highlight_genes()
        red_highlight = genes_fromenrichedgeneset()
        blue_highlight = search_gene_list()
        merged["highlight"] = merged["gene_name"].isin(highlight)
        merged["red_highlight"] = merged["gene_name"].isin(red_highlight)  
        merged["color"] = merged["highlight"].map({True: "orange", False: "gray"})
        merged['blue_highlight'] = merged['gene_name'].isin(blue_highlight)
        
        # For legend grouping
        sig_labels = {
            (True, True): "1:Yes, 2:Yes",
            (True, False): "1:Yes, 2:No",
            (False, True): "1:No, 2:Yes",
            (False, False): "1:No, 2:No",
        }
        merged["sig_label"] = merged[["significant_1", "significant_2"]].apply(lambda x: sig_labels[tuple(x)], axis=1)

        fig = go.Figure()
        point_size = input.point_size() if input.point_size() else 40

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
                        size=point_size,
                        line=dict(width=1, color="black")
                    ),
                    name=sig_labels[(sig1, sig2)],
                    customdata=subset[["gene_name", "padj_1", "padj_2"]],
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
                        size=point_size,
                        line=dict(width=2, color="black")
                    ),
                    name="Highlighted gene set",
                    customdata=subset[["gene_name", "padj_1", "padj_2"]],
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
        enrich_choice = input.enrich_significance() if "enrich_significance" in input else "both"

        if enrich_choice == "d1":
            red_mask = merged["significant_1"] & merged["red_highlight"]
        elif enrich_choice == "d2":
            red_mask = merged["significant_2"] & merged["red_highlight"]
        else:  # both
            red_mask = merged["significant_1"] & merged["significant_2"] & merged["red_highlight"]
    # Plot enriched genes (red)
        for (sig1, sig2), shape in shape_map.items():
            mask = (merged["significant_1"] == sig1) & (merged["significant_2"] == sig2) & red_mask
            subset = merged[mask]
            if not subset.empty:
                fig.add_trace(go.Scattergl(
                    x=subset["log2_fold_change_1"],
                    y=subset["log2_fold_change_2"],
                    mode="markers",
                    marker=dict(
                        symbol=shape,
                        color="red",
                        size=point_size,
                        line=dict(width=2, color="black")
                    ),
                    name="Enriched gene set (red)",
                    customdata=subset[["gene_name", "padj_1", "padj_2"]],
                    hovertemplate=(
                        "<b>Enriched Gene</b><br>"
                        "Gene: %{customdata[0]}<br>"
                        "log2FC 1: %{x:.3f}<br>"
                        "log2FC 2: %{y:.3f}<br>"
                        "padj 1: %{customdata[1]:.3g}<br>"
                        "padj 2: %{customdata[2]:.3g}<br>"
                        f"Significance: {sig_labels[(sig1, sig2)]}<br>"
                        "<extra></extra>"
                    ),
                    showlegend=True
                ))
        # Plot searched gene (blue)
        for (sig1, sig2), shape in shape_map.items():
            mask = (merged["significant_1"] == sig1) & (merged["significant_2"] == sig2) & (merged["blue_highlight"])
            subset = merged[mask]
            if not subset.empty:
                fig.add_trace(go.Scattergl(
                    x=subset["log2_fold_change_1"],
                    y=subset["log2_fold_change_2"],
                    mode="markers+text",
                    marker=dict(
                        symbol=shape,
                        color="blue",
                        size=point_size,
                        line=dict(width=3, color="black")
                    ),
                    name="Searched gene",
                    customdata=subset[["gene_name", "padj_1", "padj_2"]],
                    hovertemplate=(
                        "<b>Searched Gene</b><br>"
                        "Gene: %{customdata[0]}<br>"
                        "log2FC 1: %{x:.3f}<br>"
                        "log2FC 2: %{y:.3f}<br>"
                        "padj 1: %{customdata[1]:.3g}<br>"
                        "padj 2: %{customdata[2]:.3g}<br>"
                        "<extra></extra>"
                    ),
                    text=subset["gene_name"],
                    textposition="top center",
                    showlegend=True
                ))
        # Diagonal line
        min_val = min(merged["log2_fold_change_1"].min(), merged["log2_fold_change_2"].min())
        max_val = max(merged["log2_fold_change_1"].max(), merged["log2_fold_change_2"].max())
        if input.show_diagonal():
            fig.add_shape(
                type="line",
                x0=min_val, y0=min_val, x1=max_val, y1=max_val,
                line=dict(color="black", dash="dash"),
                layer="above"
            )

        # Crosshair axes at x=0 and y=0 if axes are shown
        if input.show_axes():
            fig.add_shape(
                type="line",
                x0=min_val, y0=0, x1=max_val, y1=0,
                line=dict(color="black", dash="dash"),
                layer="above"
            )
            fig.add_shape(
                type="line",
                x0=0, y0=min_val, x1=0, y1=max_val,
                line=dict(color="black", dash="dash"),
                layer="above"
            )
        if input.show_axes():
            fig.update_layout(
                xaxis_title=input.x_axis_label() or "Dataset 1 log2 fold change",
                yaxis_title=input.y_axis_label() or "Dataset 2 log2 fold change",
                title=input.plot_title() or "V2 Plot (Shape & Highlight)",
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