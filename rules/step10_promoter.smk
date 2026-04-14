# ===========================================================================
# STEP 10: Promoter cis-element analysis
# ===========================================================================
# Two scan methods:
#   "plantcare" (default, backward compat): user manually submits the promoter
#     FASTA to https://bioinformatics.psb.ugent.be/webtools/plantcare/html/,
#     downloads the results, and points `step10.plantcare_dir` at the folder.
#   "local": fully offline literal motif scan against a shipped motif library
#     (example/10.promoter/plantcare_motifs.tsv). Produces a synthetic
#     PlantCARE-format .tab that the existing R parser consumes unchanged.
#
# Promoter extraction is automated when `step10.compute_promoter` is true.

if STEP10_COMPUTE_PROMOTER:

    rule step10_extract_promoter:
        """Extract upstream promoter sequences from genome."""
        input:
            genome   = f"{DB_DIR}/{TARGET}.genome.fasta",
            bed      = config["step5"]["gene_bed_file"],
            gene_ids = family_ids_for(TARGET),
        output:
            fasta = f"{WORK_DIR}/10_promoter/{TARGET}_upstream.fasta",
        params:
            upstream = config["step10"].get("upstream_distance", 2000),
        log:
            f"{LOG_DIR}/10_extract_promoter.log",
        shell:
            """
            # Filter BED to target gene family members
            grep -F -f {input.gene_ids} {input.bed} > {output.fasta}.members.bed

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


if STEP10_SCAN_METHOD == "local":

    rule step10_local_motif_scan:
        """Scan promoter FASTA with a local motif library (offline PlantCARE)."""
        input:
            fasta  = promoter_fasta_for_local(),
            motifs = config["step10"].get(
                "motif_library", "example/10.promoter/plantcare_motifs.tsv"
            ),
        output:
            tab = f"{WORK_DIR}/10_promoter/local_plantcare/plantCARE_output_local.tab",
        log:
            f"{LOG_DIR}/10_local_motif_scan.log",
        shell:
            """
            python3 scripts/scan_promoter_motifs.py \
                --fasta {input.fasta} \
                --motifs {input.motifs} \
                -o {output.tab} \
                2>&1 | tee {log}
            """


def _step10_plantcare_dir():
    if STEP10_SCAN_METHOD == "local":
        return f"{WORK_DIR}/10_promoter/local_plantcare"
    return config["step10"].get(
        "plantcare_dir", "example/10.promoter/PlantCARE_410_SB_plantCARE"
    )


def _step10_plantcare_inputs():
    """Extra file-level inputs the DAG should depend on for each scan method."""
    if STEP10_SCAN_METHOD == "local":
        return [
            f"{WORK_DIR}/10_promoter/local_plantcare/plantCARE_output_local.tab",
        ]
    return []


rule step10_promoter:
    """Analyze promoter cis-elements from PlantCARE (or local scan) output."""
    input:
        element_desc   = config["step10"]["element_annotation_file"],
        gene_ids       = family_ids_for(TARGET),
        plantcare_tabs = _step10_plantcare_inputs(),
    output:
        pdf = f"{OUT_DIR}/10_promoter/promoter_elements.pdf",
    params:
        plantcare_dir = _step10_plantcare_dir(),
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
