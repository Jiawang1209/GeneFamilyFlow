# ===========================================================================
# STEP 7: Motif and gene structure analysis
# ===========================================================================

rule step07_meme:
    """Discover motifs with MEME."""
    input:
        fasta = f"{OUT_DIR}/04_identification/identify.ID.clean.fa",
    output:
        meme_txt = f"{WORK_DIR}/07_motif/meme_out/meme.txt",
    params:
        outdir  = f"{WORK_DIR}/07_motif/meme_out",
        mod     = config["step7"]["meme_mod"],
        nmotifs = config["step7"]["meme_nmotifs"],
        minw    = config["step7"]["meme_minw"],
        maxw    = config["step7"]["meme_maxw"],
    log:
        f"{LOG_DIR}/07_meme.log",
    shell:
        """
        meme {input.fasta} \
            -oc {params.outdir} \
            -mod {params.mod} -protein \
            -nmotifs {params.nmotifs} \
            -minw {params.minw} -maxw {params.maxw} \
            2>&1 | tee {log}
        """


rule step07_parse_meme:
    """Parse MEME output to extract motif info and locations."""
    input:
        meme_txt = f"{WORK_DIR}/07_motif/meme_out/meme.txt",
        fasta    = f"{OUT_DIR}/04_identification/identify.ID.clean.fa",
    output:
        info     = f"{OUT_DIR}/07_motif/meme_info.txt",
        location = f"{OUT_DIR}/07_motif/meme_location.txt",
    params:
        pattern = config["step7"].get("species_prefix", ""),
    log:
        f"{LOG_DIR}/07_parse_meme.log",
    run:
        cmd = (
            f"python3 scripts/parse_meme_output.py {input.meme_txt} "
            f"--fasta {input.fasta} "
            f"-o $(dirname {output.info})"
        )
        if params.pattern:
            cmd += f" --gene-id-pattern '{params.pattern}'"
        shell(cmd + f" 2>&1 | tee {log}")


def _step07_plot_inputs():
    inputs = {
        "location": f"{OUT_DIR}/07_motif/meme_location.txt",
        "tree": f"{WORK_DIR}/06_tree/tree.nwk",
    }
    if USE_PFAM_SCAN:
        inputs["pfam_scan"] = f"{WORK_DIR}/04_identification/Pfam_scan.out"
    return inputs


rule step07_plot:
    """Composite: Tree + Domain + Motif + Gene Structure."""
    input:
        unpack(lambda wc: _step07_plot_inputs()),
    output:
        pdf = f"{OUT_DIR}/07_motif/Tree_Domain_Motif_GeneStructure.pdf",
    params:
        group_file = config["step6"].get("group_file", ""),
        gff3       = config.get("step7", {}).get("gff3_file", ""),
    log:
        f"{LOG_DIR}/07_motif_plot.log",
    run:
        cmd = (
            f"Rscript R/07_domain_motif_structure.R "
            f"--treefile {input.tree} "
            f"--motif_file {input.location} "
            f"--outdir $(dirname {output.pdf})"
        )
        if USE_PFAM_SCAN:
            cmd += f" --pfam_scan {input.pfam_scan}"
        if params.group_file:
            cmd += f" --group_file {params.group_file}"
        if params.gff3:
            cmd += f" --gff3 {params.gff3}"
        shell(cmd + f" 2>&1 | tee {log}")
