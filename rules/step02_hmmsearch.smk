# ===========================================================================
# STEP 2: HMM search (per domain)
# ===========================================================================

rule step02_hmmsearch_pfam:
    """Search species proteins against Pfam HMM model for each domain."""
    input:
        fasta = f"{WORK_DIR}/01_database/{{species}}.longest.pep.fasta",
        hmm   = lambda wc: get_domain_config(wc.domain)["pfam_hmm_file"],
    output:
        domtblout = f"{OUT_DIR}/02_hmmsearch/{{domain}}/{{species}}.pfam.domtblout",
    log:
        f"{LOG_DIR}/02_hmmsearch_pfam_{{domain}}_{{species}}.log",
    threads: 4
    shell:
        """
        hmmsearch --cut_tc \
            --domtblout {output.domtblout} \
            -o /dev/null \
            {input.hmm} {input.fasta} \
            2>&1 | tee {log}
        """


rule step02_extract_hmm_ids:
    """Extract gene IDs passing E-value from hmmsearch results (per domain)."""
    input:
        expand(f"{OUT_DIR}/02_hmmsearch/{{{{domain}}}}/{{sp}}.pfam.domtblout", sp=SPECIES),
    output:
        ids = f"{OUT_DIR}/02_hmmsearch/{{domain}}/hmm_candidate_ids.txt",
    params:
        evalue = config["step2"]["hmmsearch_evalue"],
    shell:
        """
        python3 scripts/parse_hmmsearch.py {input} \
            -e {params.evalue} -o {output.ids}
        """


rule step02_clustalw_hmmbuild:
    """Build custom HMM from first-round HMM hits (per domain)."""
    input:
        ids   = f"{OUT_DIR}/02_hmmsearch/{{domain}}/hmm_candidate_ids.txt",
        fasta = f"{WORK_DIR}/01_database/all_species.pep.fasta",
    output:
        hmm = f"{OUT_DIR}/02_hmmsearch/{{domain}}/custom.hmm",
        aln = f"{OUT_DIR}/02_hmmsearch/{{domain}}/1st_round.aln",
    log:
        f"{LOG_DIR}/02_clustalw_hmmbuild_{{domain}}.log",
    threads: MAX_THREADS
    shell:
        """
        seqkit grep -f {input.ids} {input.fasta} -o {output.aln}.fa 2>> {log}
        clustalw -infile={output.aln}.fa -output=clustal \
            -type=PROTEIN -outfile={output.aln} >> {log} 2>&1
        hmmbuild {output.hmm} {output.aln} >> {log} 2>&1
        """


rule step02_hmmsearch_custom:
    """Search with custom HMM model — 2nd round (per domain)."""
    input:
        fasta = f"{WORK_DIR}/01_database/{{species}}.longest.pep.fasta",
        hmm   = f"{OUT_DIR}/02_hmmsearch/{{domain}}/custom.hmm",
    output:
        domtblout = f"{OUT_DIR}/02_hmmsearch/{{domain}}/{{species}}.custom.domtblout",
    log:
        f"{LOG_DIR}/02_hmmsearch_custom_{{domain}}_{{species}}.log",
    shell:
        """
        hmmsearch --domtblout {output.domtblout} \
            -o /dev/null \
            {input.hmm} {input.fasta} \
            2>&1 | tee {log}
        """


rule step02_final_hmm_ids:
    """Collect 2nd-round HMM IDs (per domain)."""
    input:
        expand(f"{OUT_DIR}/02_hmmsearch/{{{{domain}}}}/{{sp}}.custom.domtblout", sp=SPECIES),
    output:
        f"{OUT_DIR}/02_hmmsearch/{{domain}}/hmm_2nd_ids.txt",
    params:
        evalue = config["step2"]["hmmsearch_evalue"],
    shell:
        """
        python3 scripts/parse_hmmsearch.py {input} \
            -e {params.evalue} -o {output}
        """
