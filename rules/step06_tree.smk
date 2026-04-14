# ===========================================================================
# STEP 6: Phylogenetic tree construction
# ===========================================================================

rule step06_alignment:
    """Multiple sequence alignment with MUSCLE (v5 CLI)."""
    input:
        fasta = f"{OUT_DIR}/04_identification/identify.ID.clean.fa",
    output:
        aln = f"{WORK_DIR}/06_tree/alignment.muscle",
    log:
        f"{LOG_DIR}/06_muscle.log",
    threads: config["step6"].get("alignment_threads", 8)
    shell:
        """
        muscle -align {input.fasta} -output {output.aln} \
            -threads {threads} \
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
    """Visualize phylogenetic tree with ggtree — circular + rectangular."""
    input:
        tree = f"{WORK_DIR}/06_tree/tree.nwk",
    output:
        circular    = f"{OUT_DIR}/06_tree/phylogenetic_tree_circular.pdf",
        rectangular = f"{OUT_DIR}/06_tree/phylogenetic_tree_rectangular.pdf",
    params:
        species_map = build_species_map(),
        group_file = config["step6"].get("group_file", ""),
        layout     = config["step6"].get("tree_layout", "both"),
        n_subfamilies   = config["step6"].get("auto_n_subfamilies", 8),
        subfamily_prefix = config["step6"].get("auto_subfamily_prefix", "Sub"),
        outdir = f"{OUT_DIR}/06_tree",
    log:
        f"{LOG_DIR}/06_tree_plot.log",
    run:
        cmd = (
            f"Rscript R/06_tree.R "
            f"--treefile {input.tree} "
            f"--species_map '{params.species_map}' "
            f"--layout {params.layout} "
            f"--n_subfamilies {params.n_subfamilies} "
            f"--subfamily_prefix '{params.subfamily_prefix}' "
            f"--outdir {params.outdir}"
        )
        if params.group_file:
            cmd += f" --group_file {params.group_file}"
        shell(cmd + f" 2>&1 | tee {log}")
