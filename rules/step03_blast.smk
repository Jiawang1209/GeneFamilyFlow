# ===========================================================================
# STEP 3: BLAST search (per domain)
# ===========================================================================

rule step03_extract_seed:
    """Extract BLAST seed sequences for a Pfam domain from reference species.

    Uses scripts/extract_seed_by_pfam.py to filter each reference species'
    clean 2-col domains TSV by the target Pfam ID, pulls matching proteins
    from the gene-level FASTA, and concatenates the per-reference outputs.
    """
    input:
        domains_tables = lambda wc: [
            config["step3"]["reference_species"][sp]["domains_table"]
            for sp in get_seed_references(wc.domain)
        ],
        proteomes = lambda wc: [
            config["step3"]["reference_species"][sp]["proteome"]
            for sp in get_seed_references(wc.domain)
        ],
    output:
        seed = f"{OUT_DIR}/03_blast/{{domain}}/seed.fa",
    log:
        f"{LOG_DIR}/03_extract_seed_{{domain}}.log",
    run:
        import subprocess
        from pathlib import Path
        out_path = Path(output.seed)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        refs = get_seed_references(wildcards.domain)
        ref_cfg = config["step3"]["reference_species"]
        parts = []
        with open(log[0], "w") as log_fh:
            for sp in refs:
                cfg = ref_cfg[sp]
                part = out_path.parent / f"seed.{sp}.fa"
                log_fh.write(f"[extract] {sp} -> {part}\n")
                log_fh.flush()
                subprocess.run(
                    [
                        "python3", "scripts/extract_seed_by_pfam.py",
                        "--domains-table", cfg["domains_table"],
                        "--proteome",      cfg["proteome"],
                        "--pfam-id",       wildcards.domain,
                        "--output",        str(part),
                    ],
                    check=True,
                    stdout=log_fh, stderr=subprocess.STDOUT,
                )
                parts.append(part)
        with open(out_path, "w") as fh_out:
            for part in parts:
                fh_out.write(part.read_text())


rule step03_makeblastdb:
    """Build BLAST database from extracted seed sequences (per domain)."""
    input:
        seed = f"{OUT_DIR}/03_blast/{{domain}}/seed.fa",
    output:
        pin = f"{OUT_DIR}/03_blast/{{domain}}/seed_db.pin",
    params:
        db = f"{OUT_DIR}/03_blast/{{domain}}/seed_db",
    shell:
        """
        makeblastdb -in {input.seed} -dbtype prot \
            -out {params.db} -parse_seqids
        """


rule step03_blastp:
    """BLAST each species against seed database (per domain)."""
    input:
        fasta = f"{WORK_DIR}/01_database/{{species}}.longest.pep.fasta",
        pin   = f"{OUT_DIR}/03_blast/{{domain}}/seed_db.pin",
    output:
        blast = f"{OUT_DIR}/03_blast/{{domain}}/{{species}}.blast",
    params:
        db     = f"{OUT_DIR}/03_blast/{{domain}}/seed_db",
        evalue = config["step3"]["blast_evalue"],
        fmt    = config["step3"]["blast_outfmt"],
    threads: config["step3"].get("blast_num_threads", 10)
    log:
        f"{LOG_DIR}/03_blastp_{{domain}}_{{species}}.log",
    shell:
        """
        blastp -query {input.fasta} -db {params.db} \
            -outfmt '{params.fmt}' \
            -out {output.blast} \
            -evalue {params.evalue} \
            -num_threads {threads} \
            -num_alignments {config[step3][blast_num_alignments]} \
            2>&1 | tee {log}
        """


rule step03_collect_blast_ids:
    """Collect BLAST hit IDs (per domain)."""
    input:
        expand(f"{OUT_DIR}/03_blast/{{{{domain}}}}/{{sp}}.blast", sp=SPECIES),
    output:
        f"{OUT_DIR}/03_blast/{{domain}}/blast_ids.txt",
    shell:
        "cut -f1 {input} | sort -u > {output}"
