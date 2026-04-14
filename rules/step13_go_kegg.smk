# ===========================================================================
# STEP 13 (optional): GO/KEGG over-representation enrichment
# ===========================================================================

if STEP13_ENABLED:

    rule step13_extract_family_ids:
        """Per-species family member IDs = intersection of step4 final set with
        that species' longest-transcript proteome."""
        input:
            family_fa  = f"{OUT_DIR}/04_identification/identify.ID.clean.fa",
            species_fa = f"{WORK_DIR}/01_database/{{species}}.longest.pep.fasta",
        output:
            ids = f"{WORK_DIR}/13_go_kegg/{{species}}.family.ID",
        log:
            f"{LOG_DIR}/13_extract_family_ids_{{species}}.log",
        shell:
            """
            set -eo pipefail
            seqkit seq -n -i {input.species_fa} > {output.ids}.species.tmp 2>> {log}
            seqkit seq -n -i {input.family_fa}  > {output.ids}.family.tmp  2>> {log}
            grep -Fx -f {output.ids}.species.tmp {output.ids}.family.tmp \
                > {output.ids} 2>> {log} || true
            rm -f {output.ids}.species.tmp {output.ids}.family.tmp
            echo "[step13] {wildcards.species} family members: $(wc -l < {output.ids})" \
                | tee -a {log}
            """

    if STEP13_RUN_EGGNOG:
        rule step13_eggnog_run:
            """Run eggNOG-mapper on the full longest-transcript proteome so the
            enrichment universe covers every annotated gene in the species."""
            input:
                fasta = f"{WORK_DIR}/01_database/{{species}}.longest.pep.fasta",
            output:
                annotations = f"{WORK_DIR}/13_go_kegg/eggnog/{{species}}.emapper.annotations",
            params:
                data_dir = _STEP13_CFG["eggnog_data_dir"],
                out_dir  = f"{WORK_DIR}/13_go_kegg/eggnog",
                prefix   = "{species}",
            log:
                f"{LOG_DIR}/13_eggnog_{{species}}.log",
            threads: _STEP13_CFG.get("eggnog_cpu", 20)
            conda:
                "../envs/eggnog.yaml"
            shell:
                """
                mkdir -p {params.out_dir}
                emapper.py \
                    -i {input.fasta} \
                    --output {params.prefix} \
                    --output_dir {params.out_dir} \
                    --data_dir {params.data_dir} \
                    --cpu {threads} \
                    --itype proteins \
                    --override \
                    2>&1 | tee {log}
                """

    rule step13_parse_eggnog:
        """Convert eggNOG-mapper annotations into TERM2GENE tables + universe."""
        input:
            annotations = _step13_annotations,
        output:
            go       = f"{WORK_DIR}/13_go_kegg/{{species}}.go.tsv",
            kegg     = f"{WORK_DIR}/13_go_kegg/{{species}}.kegg.tsv",
            universe = f"{WORK_DIR}/13_go_kegg/{{species}}.universe.txt",
        log:
            f"{LOG_DIR}/13_parse_eggnog_{{species}}.log",
        shell:
            """
            python3 scripts/parse_eggnog.py {input.annotations} \
                --go-out {output.go} \
                --kegg-out {output.kegg} \
                --universe-out {output.universe} \
                2>&1 | tee {log}
            """

    rule step13_go_kegg:
        """Per-species GO/KEGG over-representation enrichment."""
        input:
            go        = f"{WORK_DIR}/13_go_kegg/{{species}}.go.tsv",
            kegg      = f"{WORK_DIR}/13_go_kegg/{{species}}.kegg.tsv",
            universe  = f"{WORK_DIR}/13_go_kegg/{{species}}.universe.txt",
            gene_list = f"{WORK_DIR}/13_go_kegg/{{species}}.family.ID",
        output:
            pdf = f"{OUT_DIR}/13_go_kegg/{{species}}_enrichment.pdf",
        params:
            species    = "{species}",
            outdir     = f"{OUT_DIR}/13_go_kegg",
            kegg_names = _STEP13_CFG.get("kegg_names", "") or "",
            pvalue     = _STEP13_CFG.get("pvalue", 0.05),
            qvalue     = _STEP13_CFG.get("qvalue", 0.2),
            min_size   = _STEP13_CFG.get("min_size", 5),
            top_n      = _STEP13_CFG.get("top_n", 20),
        log:
            f"{LOG_DIR}/13_go_kegg_{{species}}.log",
        shell:
            """
            mkdir -p {params.outdir}
            KEGG_ARG=""
            if [ -n "{params.kegg_names}" ] && [ -f "{params.kegg_names}" ]; then
                KEGG_ARG="--kegg_names {params.kegg_names}"
            fi
            Rscript R/13_go_kegg.R \
                --term2gene_go   {input.go} \
                --term2gene_kegg {input.kegg} \
                --gene_list      {input.gene_list} \
                --universe       {input.universe} \
                --species        {params.species} \
                --outdir         {params.outdir} \
                --pvalue         {params.pvalue} \
                --qvalue         {params.qvalue} \
                --min_size       {params.min_size} \
                --top_n          {params.top_n} \
                $KEGG_ARG \
                2>&1 | tee {log}
            """
