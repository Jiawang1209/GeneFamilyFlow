# ===========================================================================
# STEP 1: Data preparation (filter longest transcript)
# ===========================================================================

rule step01_filter_longest_transcript:
    """Filter longest transcript from protein FASTA per species."""
    input:
        fasta = f"{DB_DIR}/{{species}}.pep.fasta",
    output:
        fasta = f"{WORK_DIR}/01_database/{{species}}.longest.pep.fasta",
    log:
        f"{LOG_DIR}/01_filter_{{species}}.log",
    shell:
        """
        python3 scripts/filter_longest_transcript.py \
            {input.fasta} {output.fasta} \
            2>&1 | tee {log}
        """


rule step01_merge_proteins:
    """Merge all species protein files into one."""
    input:
        expand(f"{WORK_DIR}/01_database/{{sp}}.longest.pep.fasta", sp=SPECIES),
    output:
        f"{WORK_DIR}/01_database/all_species.pep.fasta",
    shell:
        "cat {input} > {output}"
