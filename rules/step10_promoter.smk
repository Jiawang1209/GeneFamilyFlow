# ===========================================================================
# STEP 10: Promoter cis-element analysis
# ===========================================================================
# Three scan methods (config: step10.scan_method):
#   "jaspar" (recommended): fully offline PWM scan with FIMO against the
#     JASPAR 2024 CORE plants non-redundant bundle (~805 plant TF motifs).
#     fetch_jaspar_plants.py materializes the MEME bundle, FIMO does the
#     scan, scan_promoter_fimo.py converts the result into a PlantCARE
#     8-column .tab so R/10_promoter.R consumes it unchanged. Eliminates
#     the network round-trip and is more comprehensive than the literal
#     fixture-based "local" mode.
#   "local": offline literal-substring scan against a shipped motif library
#     (example/10.promoter/plantcare_motifs.tsv, ~142 motifs). Lighter and
#     dependency-free but recall is much lower than jaspar/plantcare.
#   "plantcare" (backward compat): user manually submits the promoter
#     FASTA to https://bioinformatics.psb.ugent.be/webtools/plantcare/html/,
#     downloads the results, and points `step10.plantcare_dir` at the folder.
#     Most authoritative but blocks the pipeline on a manual web step.
#
# Promoter extraction is automated when `step10.compute_promoter` is true.

_VALID_STEP10_SCAN_METHODS = {"jaspar", "local", "plantcare"}
if STEP10_SCAN_METHOD not in _VALID_STEP10_SCAN_METHODS:
    raise ValueError(
        f"step10.scan_method must be one of {sorted(_VALID_STEP10_SCAN_METHODS)}, "
        f"got {STEP10_SCAN_METHOD!r}"
    )

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


if STEP10_SCAN_METHOD == "jaspar":

    rule step10_fetch_jaspar:
        """Materialize the JASPAR 2024 CORE plants MEME bundle (cached)."""
        output:
            meme = config["step10"].get(
                "jaspar_meme", "example/10.promoter/JASPAR2024_plants.meme"
            ),
        params:
            url = config["step10"].get(
                "jaspar_url",
                "https://jaspar.elixir.no/download/data/2024/CORE/"
                "JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt",
            ),
        log:
            f"{LOG_DIR}/10_fetch_jaspar.log",
        shell:
            """
            python3 scripts/fetch_jaspar_plants.py \
                --output {output.meme} \
                --url '{params.url}' \
                2>&1 | tee {log}
            """

    rule step10_fimo_scan:
        """Scan promoter FASTA with FIMO against the JASPAR plants bundle."""
        input:
            fasta = promoter_fasta_for_local(),
            meme  = config["step10"].get(
                "jaspar_meme", "example/10.promoter/JASPAR2024_plants.meme"
            ),
        output:
            tab = f"{WORK_DIR}/10_promoter/jaspar_fimo/plantCARE_output_jaspar.tab",
        params:
            threshold = config["step10"].get("fimo_threshold", 1e-4),
        log:
            f"{LOG_DIR}/10_fimo_scan.log",
        shell:
            """
            python3 scripts/scan_promoter_fimo.py \
                --fasta {input.fasta} \
                --meme {input.meme} \
                --threshold {params.threshold} \
                -o {output.tab} \
                2>&1 | tee {log}
            """


def _step10_plantcare_dir():
    if STEP10_SCAN_METHOD == "local":
        return f"{WORK_DIR}/10_promoter/local_plantcare"
    if STEP10_SCAN_METHOD == "jaspar":
        return f"{WORK_DIR}/10_promoter/jaspar_fimo"
    return config["step10"].get(
        "plantcare_dir", "example/10.promoter/PlantCARE_410_SB_plantCARE"
    )


def _step10_plantcare_inputs():
    """Extra file-level inputs the DAG should depend on for each scan method."""
    if STEP10_SCAN_METHOD == "local":
        return [
            f"{WORK_DIR}/10_promoter/local_plantcare/plantCARE_output_local.tab",
        ]
    if STEP10_SCAN_METHOD == "jaspar":
        return [
            f"{WORK_DIR}/10_promoter/jaspar_fimo/plantCARE_output_jaspar.tab",
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
