# ===========================================================================
# STEP 10: Promoter cis-element analysis
# ===========================================================================
# Promoter extraction is automated; PlantCARE submission is manual (online).
# Set step10.plantcare_dir to the downloaded PlantCARE results directory.

if STEP10_COMPUTE_PROMOTER:

    rule step10_extract_promoter:
        """Extract upstream promoter sequences from genome."""
        input:
            genome = f"{DB_DIR}/{TARGET}.genome.fasta",
            bed    = config["step5"]["gene_bed_file"],
        output:
            fasta = f"{WORK_DIR}/10_promoter/{TARGET}_upstream.fasta",
        params:
            upstream = config["step10"].get("upstream_distance", 2000),
            gene_ids = config["step9"]["species_config"][TARGET]["npf_id_file"],
        log:
            f"{LOG_DIR}/10_extract_promoter.log",
        shell:
            """
            # Filter BED to target gene family members
            grep -F -f {params.gene_ids} {input.bed} > {output.fasta}.members.bed

            # Generate upstream region BED
            awk -v up={params.upstream} 'BEGIN{{OFS="\\t"}} {{
                if ($6 == "+") {{
                    s = $2 - up; if (s < 0) s = 0;
                    print $1, s, $2, $4, 0, $6
                }} else {{
                    print $1, $3, $3 + up, $4, 0, $6
                }}
            }}' {output.fasta}.members.bed > {output.fasta}.upstream.bed

            # Extract sequences
            bedtools getfasta -fi {input.genome} -bed {output.fasta}.upstream.bed \
                -fo {output.fasta} -name -s 2>&1 | tee {log}

            rm -f {output.fasta}.members.bed {output.fasta}.upstream.bed
            """


rule step10_promoter:
    """Analyze promoter cis-elements from PlantCARE output."""
    input:
        element_desc = config["step10"]["element_annotation_file"],
        gene_ids     = config["step9"]["species_config"][TARGET]["npf_id_file"],
    output:
        pdf = f"{OUT_DIR}/10_promoter/promoter_elements.pdf",
    params:
        plantcare_dir = config["step10"].get("plantcare_dir",
                            "example/10.promoter/PlantCARE_410_SB_plantCARE"),
    log:
        f"{LOG_DIR}/10_promoter.log",
    shell:
        """
        Rscript R/10_promoter.R \
            --plantcare_dir {params.plantcare_dir} \
            --element_desc {input.element_desc} \
            --gene_ids {input.gene_ids} \
            --outdir $(dirname {output.pdf}) \
            2>&1 | tee {log}
        """
