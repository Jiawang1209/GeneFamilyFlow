# ===========================================================================
# STEP 14 (optional): qRT-PCR relative expression
# ===========================================================================

if STEP14_ENABLED:
    rule step14_qrt_pcr:
        """qRT-PCR relative expression bar plot with optional significance stars."""
        input:
            expression = config["step14_qrtpcr"]["expression_data"],
        output:
            pdf = f"{OUT_DIR}/14_qrt_pcr/qRT_PCR.pdf",
            summary = f"{OUT_DIR}/14_qrt_pcr/qRT_PCR_summary.csv",
        params:
            sheet       = config["step14_qrtpcr"].get("sheet", "Sheet5"),
            id_col      = config["step14_qrtpcr"].get("id_col", "ID"),
            group_cols  = ",".join(config["step14_qrtpcr"].get("group_cols", []) or []),
            comparisons = ",".join(config["step14_qrtpcr"].get("comparisons", []) or []),
            ncol        = config["step14_qrtpcr"].get("facet_ncol", 4),
            width       = config["step14_qrtpcr"].get("figure_width", 12),
            height      = config["step14_qrtpcr"].get("figure_height", 8),
        log:
            f"{LOG_DIR}/14_qrt_pcr.log",
        shell:
            """
            Rscript R/14_qrt_pcr.R \
                --expression {input.expression} \
                --sheet "{params.sheet}" \
                --id_col "{params.id_col}" \
                --group_cols "{params.group_cols}" \
                --comparisons "{params.comparisons}" \
                --outdir $(dirname {output.pdf}) \
                --ncol {params.ncol} \
                --width {params.width} \
                --height {params.height} \
                2>&1 | tee {log}
            """
