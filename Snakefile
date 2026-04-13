"""
GeneFamilyFlow — Multi-species gene family analysis pipeline.

Usage:
    snakemake --configfile config/default_config.yaml -n -p           # dry-run
    snakemake --configfile config/default_config.yaml -j 8 --use-conda # full run
    snakemake --configfile config/default_config.yaml --until step05   # partial run
"""

configfile: "config/default_config.yaml"

# ---------------------------------------------------------------------------
# Globals from config
# ---------------------------------------------------------------------------
SPECIES       = config["species"]
TARGET        = config["target_species"]
DB_DIR        = config["database_dir"]
OUT_DIR       = config["output_dir"]
WORK_DIR      = config["work_dir"]
LOG_DIR       = config.get("log_dir", "logs")
MAX_THREADS   = config.get("max_threads", 8)

# ---------------------------------------------------------------------------
# Helper: species-pair combinations for JCVI synteny
# ---------------------------------------------------------------------------
import itertools
SPECIES_PAIRS = list(itertools.combinations(SPECIES, 2))

# ---------------------------------------------------------------------------
# Multi-domain support
# ---------------------------------------------------------------------------
DOMAINS     = config["step2"]["domains"]
DOMAIN_IDS  = [d["pfam_id"] for d in DOMAINS]

def get_domain_config(pfam_id):
    """Look up per-domain config by Pfam ID."""
    for d in DOMAINS:
        if d["pfam_id"] == pfam_id:
            return d
    raise ValueError(f"Domain {pfam_id} not found in config step2.domains")


def get_seed_references(pfam_id):
    """Return the list of reference species keys for a domain's BLAST seed.

    If the domain entry omits `seed_references`, fall back to every species
    declared in step3.reference_species.
    """
    domain_cfg = get_domain_config(pfam_id)
    all_refs = config["step3"]["reference_species"]
    refs = domain_cfg.get("seed_references")
    if refs is None:
        return list(all_refs.keys())
    missing = [r for r in refs if r not in all_refs]
    if missing:
        raise ValueError(
            f"Domain {pfam_id} references unknown species {missing}; "
            f"add them to step3.reference_species."
        )
    return refs

# Tree tool
TREE_TOOL = config["step6"].get("tree_tool", "iqtree")

# ---------------------------------------------------------------------------
# rule all — define final targets
# ---------------------------------------------------------------------------
rule all:
    input:
        # Step 4 (final gene family set)
        f"{OUT_DIR}/04_identification/identify.ID.clean.fa",
        # Step 5
        f"{OUT_DIR}/05_genefamily_info/Gene_Information.xlsx",
        # Step 6
        f"{OUT_DIR}/06_tree/phylogenetic_tree.pdf",
        # Step 7
        f"{OUT_DIR}/07_motif/meme_location.txt",
        f"{OUT_DIR}/07_motif/Tree_Domain_Motif_GeneStructure.pdf",
        # Step 8
        f"{OUT_DIR}/08_jcvi_kaks/08_jcvi_kaks.pdf",
        # Step 9
        expand(f"{OUT_DIR}/09_circos/circos_{{sp}}.pdf", sp=SPECIES),
        # Step 10
        f"{OUT_DIR}/10_promoter/promoter_elements.pdf",
        # Step 11
        f"{OUT_DIR}/11_ppi/PPI_network.pdf",


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


# ===========================================================================
# STEP 5: Gene family information statistics
# ===========================================================================

rule step05:
    """Compute physicochemical properties and gene info."""
    input:
        fasta = config["step5"]["gene_fasta_file"],
        bed   = config["step5"]["gene_bed_file"],
    output:
        xlsx = f"{OUT_DIR}/05_genefamily_info/Gene_Information.xlsx",
        csv  = f"{OUT_DIR}/05_genefamily_info/Gene_Information.csv",
    params:
        species_map = ",".join(
            f"{sp[:2] if sp != 'Osativa' else 'LOC_Os'}={sp}"
            for sp in SPECIES
        ),
    log:
        f"{LOG_DIR}/05_genefamily_info.log",
    shell:
        """
        Rscript R/05_genefamily_info.R \
            --input_fasta {input.fasta} \
            --input_bed {input.bed} \
            --species_map "{params.species_map}" \
            --outdir $(dirname {output.xlsx}) \
            2>&1 | tee {log}
        """


# ===========================================================================
# STEP 6: Phylogenetic tree construction
# ===========================================================================

rule step06_alignment:
    """Multiple sequence alignment with MUSCLE."""
    input:
        fasta = config["step5"]["gene_fasta_file"],
    output:
        aln = f"{WORK_DIR}/06_tree/alignment.muscle",
    log:
        f"{LOG_DIR}/06_muscle.log",
    threads: config["step6"].get("alignment_threads", 8)
    shell:
        """
        muscle -in {input.fasta} -out {output.aln} \
            2>&1 | tee {log}
        """


if TREE_TOOL == "iqtree":
    rule step06_tree_build:
        """Build phylogenetic tree with IQ-TREE."""
        input:
            aln = f"{WORK_DIR}/06_tree/alignment.muscle",
        output:
            tree = f"{WORK_DIR}/06_tree/tree.nwk",
        params:
            prefix = f"{WORK_DIR}/06_tree/iqtree_out",
            model     = config["step6"]["iqtree_model"],
            bootstrap = config["step6"]["iqtree_bootstrap"],
            seed      = config["step6"].get("iqtree_seed", 12345),
        log:
            f"{LOG_DIR}/06_tree_build.log",
        threads: config["step6"].get("iqtree_threads", 8)
        shell:
            """
            iqtree -s {input.aln} \
                -m {params.model} \
                -bb {params.bootstrap} -bnni \
                -nt {threads} \
                -seed {params.seed} \
                --prefix {params.prefix} \
                -redo \
                2>&1 | tee {log}
            cp {params.prefix}.treefile {output.tree}
            """
else:
    rule step06_tree_build:
        """Build phylogenetic tree with FastTree."""
        input:
            aln = f"{WORK_DIR}/06_tree/alignment.muscle",
        output:
            tree = f"{WORK_DIR}/06_tree/tree.nwk",
        params:
            options = config["step6"].get("fasttree_options", "-gamma"),
        log:
            f"{LOG_DIR}/06_tree_build.log",
        shell:
            """
            FastTree {params.options} {input.aln} > {output.tree} \
                2> {log}
            """


rule step06_tree_plot:
    """Visualize phylogenetic tree with ggtree."""
    input:
        tree = f"{WORK_DIR}/06_tree/tree.nwk",
    output:
        pdf = f"{OUT_DIR}/06_tree/phylogenetic_tree.pdf",
    params:
        species_map = ",".join(
            f"{sp[:2] if sp != 'Osativa' else 'LOC_Os'}={sp}"
            for sp in SPECIES
        ),
        group_file = config["step6"].get("group_file", ""),
        layout     = config["step6"].get("tree_layout", "circular"),
    log:
        f"{LOG_DIR}/06_tree_plot.log",
    run:
        cmd = (
            f"Rscript R/06_tree.R "
            f"--treefile {input.tree} "
            f"--species_map '{params.species_map}' "
            f"--layout {params.layout} "
            f"--outdir $(dirname {output.pdf})"
        )
        if params.group_file:
            cmd += f" --group_file {params.group_file}"
        shell(cmd + f" 2>&1 | tee {log}")


# ===========================================================================
# STEP 7: Motif and gene structure analysis
# ===========================================================================

rule step07_meme:
    """Discover motifs with MEME."""
    input:
        fasta = config["step5"]["gene_fasta_file"],
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
            f"-o $(dirname {output.info})"
        )
        if params.pattern:
            cmd += f" --gene-id-pattern '{params.pattern}'"
        shell(cmd + f" 2>&1 | tee {log}")


rule step07_plot:
    """Composite: Tree + Domain + Motif + Gene Structure."""
    input:
        location  = f"{OUT_DIR}/07_motif/meme_location.txt",
        tree      = f"{WORK_DIR}/06_tree/tree.nwk",
        pfam_scan = f"{WORK_DIR}/04_identification/Pfam_scan.out",
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
            f"--pfam_scan {input.pfam_scan} "
            f"--outdir $(dirname {output.pdf})"
        )
        if params.group_file:
            cmd += f" --group_file {params.group_file}"
        if params.gff3:
            cmd += f" --gff3 {params.gff3}"
        shell(cmd + f" 2>&1 | tee {log}")


# ===========================================================================
# STEP 8: JCVI synteny and Ka/Ks
# ===========================================================================
# Set step8.precomputed = false to run JCVI from scratch.
# Computation requires: {DB_DIR}/{species}.gff3, {DB_DIR}/{species}.cds.fasta

STEP8_PRECOMPUTED = config["step8"].get("precomputed", True)
STEP8_KAKS_FILE = (
    config["step8"]["kaks_file"] if STEP8_PRECOMPUTED
    else f"{WORK_DIR}/08_jcvi/kaks.tab.xls"
)

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
        species_map = ",".join(
            f"{sp[:2] if sp != 'Osativa' else 'LOC_Os'}={sp}"
            for sp in SPECIES
        ),
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


# ===========================================================================
# STEP 9: MCScanX synteny and Circos
# ===========================================================================
# Set step9.precomputed = false to run MCScanX from scratch.
# Computation requires: {DB_DIR}/{species}.gff3

STEP9_PRECOMPUTED = config["step9"].get("precomputed", True)

def _step09_file(wc, ext):
    """Resolve MCScanX output path (precomputed or computed)."""
    if STEP9_PRECOMPUTED:
        return f"example/9.mcscanx/{wc.species}.{ext}"
    return f"{WORK_DIR}/09_mcscanx/{wc.species}.{ext}"

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
        gene_ids     = lambda wc: config["step9"]["species_config"][wc.species]["npf_id_file"],
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


# ===========================================================================
# STEP 10: Promoter cis-element analysis
# ===========================================================================
# Promoter extraction is automated; PlantCARE submission is manual (online).
# Set step10.plantcare_dir to the downloaded PlantCARE results directory.

STEP10_COMPUTE_PROMOTER = config["step10"].get("compute_promoter", False)

if STEP10_COMPUTE_PROMOTER:

    rule step10_extract_promoter:
        """Extract upstream promoter sequences from genome."""
        input:
            genome = f"{DB_DIR}/{TARGET}.genome.fasta",
            bed    = config["step5"]["gene_bed_file"],
        output:
            fasta = f"{WORK_DIR}/10_promoter/{TARGET}_upstream.fasta",
        params:
            upstream = config["step10"].get("upstream_distance", 2000),
            gene_ids = config["step9"]["species_config"][TARGET]["npf_id_file"],
        log:
            f"{LOG_DIR}/10_extract_promoter.log",
        shell:
            """
            # Filter BED to target gene family members
            grep -F -f {params.gene_ids} {input.bed} > {output.fasta}.members.bed

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


rule step10_promoter:
    """Analyze promoter cis-elements from PlantCARE output."""
    input:
        element_desc = config["step10"]["element_annotation_file"],
        gene_ids     = config["step9"]["species_config"][TARGET]["npf_id_file"],
    output:
        pdf = f"{OUT_DIR}/10_promoter/promoter_elements.pdf",
    params:
        plantcare_dir = config["step10"].get("plantcare_dir",
                            "example/10.promoter/PlantCARE_410_SB_plantCARE"),
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


# ===========================================================================
# STEP 11: PPI network
# ===========================================================================

rule step11_ppi:
    """Protein-protein interaction network visualization."""
    input:
        edge_file = config["step11"]["ppi_edge_file"],
    output:
        pdf = f"{OUT_DIR}/11_ppi/PPI_network.pdf",
    params:
        layout      = config["step11"]["igraph_layout"],
        top_modules = config["step11"].get("top_modules", 15),
    log:
        f"{LOG_DIR}/11_ppi.log",
    shell:
        """
        Rscript R/11_ppi.R \
            --edge_file {input.edge_file} \
            --outdir $(dirname {output.pdf}) \
            --layout {params.layout} \
            --top_modules {params.top_modules} \
            2>&1 | tee {log}
        """


# ===========================================================================
# OPTIONAL: Step 13-14 (GO/KEGG, qRT-PCR)
# ===========================================================================

if config.get("step13_go_kegg", {}).get("enabled", False):
    rule step13_go_kegg:
        input:
            config["step13_go_kegg"]["annotation_file"],
        output:
            f"{OUT_DIR}/13_go_kegg/enrichment.pdf",
        shell:
            """
            Rscript R/13_go_kegg.R \
                --input {input} \
                --outdir $(dirname {output})
            """

if config.get("step14_qrtpcr", {}).get("enabled", False):
    rule step14_qrtpcr:
        input:
            config["step14_qrtpcr"]["expression_data"],
        output:
            f"{OUT_DIR}/14_qrtpcr/qrtpcr_plot.pdf",
        shell:
            """
            Rscript R/14_qrt_pcr.R \
                --input {input} \
                --outdir $(dirname {output})
            """
