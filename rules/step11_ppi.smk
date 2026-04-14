# ===========================================================================
# STEP 11: PPI network
# ===========================================================================

rule step11_ppi:
    """Protein-protein interaction network visualization."""
    input:
        edge_file = config["step11"]["ppi_edge_file"],
    output:
        pdf = f"{OUT_DIR}/11_ppi/PPI_network.pdf",
    params:
        layout      = config["step11"]["igraph_layout"],
        top_modules = config["step11"].get("top_modules", 15),
    log:
        f"{LOG_DIR}/11_ppi.log",
    shell:
        """
        Rscript R/11_ppi.R \
            --edge_file {input.edge_file} \
            --outdir $(dirname {output.pdf}) \
            --layout {params.layout} \
            --top_modules {params.top_modules} \
            2>&1 | tee {log}
        """
