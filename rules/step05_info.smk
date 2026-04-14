# ===========================================================================
# STEP 5: Gene family information statistics
# ===========================================================================

rule step05_genefamily_info:
    """Compute physicochemical properties and gene info."""
    input:
        fasta = f"{OUT_DIR}/04_identification/identify.ID.clean.fa",
        bed   = config["step5"]["gene_bed_file"],
    output:
        xlsx = f"{OUT_DIR}/05_genefamily_info/Gene_Information.xlsx",
        csv  = f"{OUT_DIR}/05_genefamily_info/Gene_Information.csv",
    params:
        species_map = build_species_map(),
    log:
        f"{LOG_DIR}/05_genefamily_info.log",
    shell:
        """
        Rscript R/05_genefamily_info.R \
            --input_fasta {input.fasta} \
            --input_bed {input.bed} \
            --species_map "{params.species_map}" \
            --outdir $(dirname {output.xlsx}) \
            2>&1 | tee {log}
        """
