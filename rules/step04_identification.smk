# ===========================================================================
# STEP 4: Gene family identification (merge across domains + Pfam verify)
# ===========================================================================

rule step04_merge_per_domain:
    """Merge HMM and BLAST candidates per domain (intersection or union)."""
    input:
        hmm   = f"{OUT_DIR}/02_hmmsearch/{{domain}}/hmm_2nd_ids.txt",
        blast = f"{OUT_DIR}/03_blast/{{domain}}/blast_ids.txt",
    output:
        ids = f"{WORK_DIR}/04_identification/{{domain}}.merged.ID",
    params:
        method = config["step4"]["merge_method"],
    shell:
        """
        python3 scripts/merge_gene_ids.py \
            {input.hmm} {input.blast} \
            -m {params.method} -o {output.ids}
        """


rule step04_combine_domains:
    """Union all per-domain gene IDs and extract FASTA."""
    input:
        ids   = expand(f"{WORK_DIR}/04_identification/{{dom}}.merged.ID", dom=DOMAIN_IDS),
        fasta = f"{WORK_DIR}/01_database/all_species.pep.fasta",
    output:
        ids   = f"{WORK_DIR}/04_identification/all_domains.ID",
        fasta = f"{WORK_DIR}/04_identification/all_domains.ID.fa",
    shell:
        """
        cat {input.ids} | sort -u > {output.ids}
        seqkit grep -f {output.ids} {input.fasta} -o {output.fasta}
        """


if USE_PFAM_SCAN:

    rule step04_pfam_scan:
        """Run pfam_scan.pl on combined candidates."""
        input:
            fasta = f"{WORK_DIR}/04_identification/all_domains.ID.fa",
        output:
            out = f"{WORK_DIR}/04_identification/Pfam_scan.out",
        params:
            pfam_dir = config["step4"]["pfam_database_dir"],
            evalue   = config["step4"]["pfam_scan_evalue"],
            cpu      = config["step4"]["pfam_scan_cpu"],
        log:
            f"{LOG_DIR}/04_pfam_scan.log",
        threads: config["step4"].get("pfam_scan_cpu", 30)
        shell:
            """
            pfam_scan.pl -fasta {input.fasta} \
                -dir {params.pfam_dir} \
                -cpu {params.cpu} \
                -e_seq {params.evalue} \
                -out {output.out} 2>&1 | tee {log}
            """


    rule step04_pfam_filter:
        """Verify candidates by pfam_scan and combine target domains.

        The `domain_combine` mode decides how multi-domain results are merged:
          - "any" → gene kept if it has at least one listed domain (union)
          - "all" → gene kept only if it has every listed domain (intersection)
        """
        input:
            pfam_out = f"{WORK_DIR}/04_identification/Pfam_scan.out",
            fasta    = f"{WORK_DIR}/04_identification/all_domains.ID.fa",
        output:
            ids   = f"{WORK_DIR}/04_identification/pfam_scan.id",
            fasta = f"{OUT_DIR}/04_identification/identify.ID.clean.fa",
        params:
            pfam_ids = ",".join(DOMAIN_IDS),
            mode     = config["step4"].get("domain_combine", "any"),
        shell:
            """
            python3 scripts/parse_pfam_scan.py {input.pfam_out} \
                --pfam-id {params.pfam_ids} \
                --mode {params.mode} \
                -o {output.ids}
            seqkit grep -f {output.ids} {input.fasta} -o {output.fasta}
            """

else:

    rule step04_hmm_filter:
        """Skip pfam_scan.pl — trust per-domain hmmsearch hits directly.

        Applies the `domain_combine` mode across per-domain ID files from step 2:
          - "any" → union of per-domain hits
          - "all" → intersection of per-domain hits (all domains must be present)
        """
        input:
            per_domain_ids = expand(
                f"{WORK_DIR}/04_identification/{{dom}}.merged.ID", dom=DOMAIN_IDS
            ),
            fasta = f"{WORK_DIR}/04_identification/all_domains.ID.fa",
        output:
            ids   = f"{WORK_DIR}/04_identification/pfam_scan.id",
            fasta = f"{OUT_DIR}/04_identification/identify.ID.clean.fa",
        params:
            mode = config["step4"].get("domain_combine", "any"),
        log:
            f"{LOG_DIR}/04_hmm_filter.log",
        run:
            from pathlib import Path

            per_domain_sets = [
                set(Path(f).read_text().split()) for f in input.per_domain_ids
            ]
            if not per_domain_sets:
                keep = set()
            elif params.mode == "all":
                keep = set.intersection(*per_domain_sets)
            else:
                keep = set.union(*per_domain_sets)

            Path(output.ids).write_text("\n".join(sorted(keep)) + "\n")
            with open(log[0], "w") as fh:
                fh.write(
                    f"domain_combine mode: {params.mode}\n"
                    f"per-domain counts: "
                    f"{[len(s) for s in per_domain_sets]}\n"
                    f"kept ({params.mode}): {len(keep)}\n"
                )
            shell(
                f"seqkit grep -f {output.ids} {input.fasta} -o {output.fasta} "
                f"2>> {log[0]}"
            )
