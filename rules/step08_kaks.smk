# ===========================================================================
# STEP 8: JCVI synteny and Ka/Ks
# ===========================================================================
# Set step8.precomputed = false to run JCVI from scratch.
# Computation requires: {DB_DIR}/{species}.gff3, {DB_DIR}/{species}.cds.fasta

if not STEP8_PRECOMPUTED:

    rule step08_gff_to_bed:
        """Convert GFF3 to JCVI BED format."""
        input:
            gff3 = f"{DB_DIR}/{{species}}.gff3",
        output:
            bed = f"{WORK_DIR}/08_jcvi/{{species}}.bed",
        log:
            f"{LOG_DIR}/08_gff2bed_{{species}}.log",
        shell:
            """
            python -m jcvi.formats.gff bed --type=mRNA --key=Name \
                {input.gff3} -o {output.bed} 2>&1 | tee {log}
            """

    rule step08_jcvi_ortholog:
        """JCVI pairwise ortholog + synteny calling."""
        input:
            bed_a = f"{WORK_DIR}/08_jcvi/{{sp_a}}.bed",
            bed_b = f"{WORK_DIR}/08_jcvi/{{sp_b}}.bed",
            cds_a = f"{DB_DIR}/{{sp_a}}.cds.fasta",
            cds_b = f"{DB_DIR}/{{sp_b}}.cds.fasta",
        output:
            anchors = f"{WORK_DIR}/08_jcvi/{{sp_a}}.{{sp_b}}.anchors",
            simple  = f"{WORK_DIR}/08_jcvi/{{sp_a}}.{{sp_b}}.anchors.simple",
        params:
            workdir = f"{WORK_DIR}/08_jcvi",
            dbtype  = config["step8"].get("jcvi_database_type", "prot"),
            evalue  = config["step8"].get("jcvi_evalue", "1e-10"),
            minspan = config["step8"].get("synteny_minspan", 30),
            no_strip = "--no_strip_names" if config["step8"].get("jcvi_no_strip_names", True) else "",
        log:
            f"{LOG_DIR}/08_jcvi_{{sp_a}}_{{sp_b}}.log",
        threads: MAX_THREADS
        shell:
            """
            mkdir -p {params.workdir}
            ln -sf $(realpath {input.cds_a}) {params.workdir}/{wildcards.sp_a}.cds
            ln -sf $(realpath {input.cds_b}) {params.workdir}/{wildcards.sp_b}.cds
            cd {params.workdir} && \
            python -m jcvi.compara.catalog ortholog \
                {wildcards.sp_a} {wildcards.sp_b} \
                {params.no_strip} --dbtype={params.dbtype} \
                2>&1 | tee $(realpath {log})
            cd {params.workdir} && \
            python -m jcvi.compara.synteny screen \
                --minspan={params.minspan} --simple \
                {wildcards.sp_a}.{wildcards.sp_b}.anchors \
                {wildcards.sp_a}.{wildcards.sp_b}.anchors.new \
                2>> $(realpath {log})
            """

    rule step08_kaks_calc:
        """Calculate Ka/Ks for syntenic gene pairs using JCVI."""
        input:
            simple = [
                f"{WORK_DIR}/08_jcvi/{a}.{b}.anchors.simple"
                for a, b in SPECIES_PAIRS
            ],
            cds = expand(f"{DB_DIR}/{{sp}}.cds.fasta", sp=SPECIES),
        output:
            kaks = f"{WORK_DIR}/08_jcvi/kaks.tab.xls",
        params:
            workdir = f"{WORK_DIR}/08_jcvi",
        log:
            f"{LOG_DIR}/08_kaks_calc.log",
        shell:
            """
            cd {params.workdir}
            echo -e "Seq_1\\tSeq_2\\tKa\\tKs\\tKa_Ks" > {output.kaks}
            for f in *.anchors.simple; do
                python -m jcvi.apps.ks kaks "$f" 2>> $(realpath {log}) || true
            done
            for f in *.kaks; do
                [ -f "$f" ] && tail -n +2 "$f" >> {output.kaks} 2>/dev/null || true
            done
            """


rule step08_kaks_plot:
    """Ka/Ks distribution boxplot from JCVI collinearity results."""
    input:
        kaks = STEP8_KAKS_FILE,
    output:
        pdf = f"{OUT_DIR}/08_jcvi_kaks/08_jcvi_kaks.pdf",
    params:
        species_map = build_species_map(),
    log:
        f"{LOG_DIR}/08_jcvi_kaks.log",
    shell:
        """
        Rscript R/08_jcvi_kaks.R \
            --kaks_file {input.kaks} \
            --species_map "{params.species_map}" \
            --outdir $(dirname {output.pdf}) \
            2>&1 | tee {log}
        """
