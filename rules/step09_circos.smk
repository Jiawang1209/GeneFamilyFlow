# ===========================================================================
# STEP 9: MCScanX synteny and Circos
# ===========================================================================
# Set step9.precomputed = false to run MCScanX from scratch.
# Computation requires: {DB_DIR}/{species}.gff3

if not STEP9_PRECOMPUTED:

    rule step09_prep_gff:
        """Convert GFF3 to MCScanX simplified GFF format."""
        input:
            gff3 = f"{DB_DIR}/{{species}}.gff3",
        output:
            gff = f"{WORK_DIR}/09_mcscanx/{{species}}.gff",
        shell:
            r"""
            awk -F'\t' '$3=="gene" {{
                id=""; split($9, attrs, ";");
                for(i in attrs) {{
                    if(attrs[i] ~ /^ID=/) {{
                        sub(/^ID=/, "", attrs[i]);
                        id=attrs[i]
                    }}
                }};
                if(id != "") print $1"\t"id"\t"$4"\t"$5
            }}' {input.gff3} > {output.gff}
            """

    rule step09_self_blast:
        """All-vs-all BLAST within species for MCScanX."""
        input:
            fasta = f"{WORK_DIR}/01_database/{{species}}.longest.pep.fasta",
        output:
            blast = f"{WORK_DIR}/09_mcscanx/{{species}}.blast",
        params:
            db = f"{WORK_DIR}/09_mcscanx/{{species}}_db",
        threads: config["step3"].get("blast_num_threads", 10)
        log:
            f"{LOG_DIR}/09_self_blast_{{species}}.log",
        shell:
            """
            makeblastdb -in {input.fasta} -dbtype prot \
                -out {params.db} 2>> {log}
            blastp -query {input.fasta} -db {params.db} \
                -outfmt 6 -evalue 1e-10 \
                -num_threads {threads} \
                -num_alignments 5 \
                -out {output.blast} \
                2>&1 | tee -a {log}
            """

    rule step09_mcscanx:
        """Run MCScanX for synteny detection."""
        input:
            gff   = f"{WORK_DIR}/09_mcscanx/{{species}}.gff",
            blast = f"{WORK_DIR}/09_mcscanx/{{species}}.blast",
        output:
            collinearity = f"{WORK_DIR}/09_mcscanx/{{species}}.collinearity",
            tandem       = f"{WORK_DIR}/09_mcscanx/{{species}}.tandem",
        params:
            prefix  = f"{WORK_DIR}/09_mcscanx/{{species}}",
            mcscanx = config["step9"]["mcscanx_path"],
        log:
            f"{LOG_DIR}/09_mcscanx_{{species}}.log",
        shell:
            "{params.mcscanx} {params.prefix} 2>&1 | tee {log}"

    rule step09_dup_classifier:
        """Classify duplicate genes with MCScanX."""
        input:
            gff          = f"{WORK_DIR}/09_mcscanx/{{species}}.gff",
            blast        = f"{WORK_DIR}/09_mcscanx/{{species}}.blast",
            collinearity = f"{WORK_DIR}/09_mcscanx/{{species}}.collinearity",
        output:
            gene_type = f"{WORK_DIR}/09_mcscanx/{{species}}.gene_type",
        params:
            prefix     = f"{WORK_DIR}/09_mcscanx/{{species}}",
            classifier = config["step9"]["duplicate_classifier_path"],
        log:
            f"{LOG_DIR}/09_dup_classifier_{{species}}.log",
        shell:
            "{params.classifier} {params.prefix} 2>&1 | tee {log}"


rule step09_circos:
    """Per-species Circos plot from MCScanX results."""
    input:
        chr_fai      = lambda wc: config["step9"]["species_config"][wc.species]["genome_length_file"],
        bed          = config["step5"]["gene_bed_file"],
        gene_ids     = lambda wc: family_ids_for(wc.species),
        gene_type    = lambda wc: _step09_file(wc, "gene_type"),
        tandem       = lambda wc: _step09_file(wc, "tandem"),
        collinearity = lambda wc: _step09_file(wc, "collinearity"),
    output:
        pdf = f"{OUT_DIR}/09_circos/circos_{{species}}.pdf",
    params:
        n_chr = lambda wc: config["step9"]["species_config"][wc.species]["n_chromosomes"],
    log:
        f"{LOG_DIR}/09_circos_{{species}}.log",
    shell:
        """
        Rscript R/09_circos.R \
            --chr_fai {input.chr_fai} \
            --bed {input.bed} \
            --gene_ids {input.gene_ids} \
            --gene_type {input.gene_type} \
            --tandem {input.tandem} \
            --collinearity {input.collinearity} \
            --outdir $(dirname {output.pdf}) \
            --n_chr {params.n_chr} \
            --output_name $(basename {output.pdf}) \
            2>&1 | tee {log}
        """
