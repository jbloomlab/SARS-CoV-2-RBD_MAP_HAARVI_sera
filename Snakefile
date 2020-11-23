"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import glob
import itertools
import os.path
import textwrap
import urllib.request

import Bio.SeqIO

import dms_variants.codonvarianttable
import dms_variants.illuminabarcodeparser

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: make_rulegraph,
            make_summary

# Functions -------------------------------------------------------------------
def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# Information on samples and barcode runs -------------------------------------
barcode_runs = pd.read_csv(config['barcode_runs'])

# combination of the *library* and *sample* columns should be unique.
assert len(barcode_runs.groupby(['library', 'sample'])) == len(barcode_runs)

# *sample* should be the hyphen separated concatenation of
# *experiment*, *antibody*, *concentration*, and *selection*.
sample_vs_expect = (
    barcode_runs
    .assign(expect=lambda x: x[['experiment', 'antibody', 'concentration',
                                'selection']]
                             .apply(lambda r: '-'.join(r.values.astype(str)),
                                    axis=1),
            equal=lambda x: x['sample'] == x['expect'],
            )
    )
assert sample_vs_expect['equal'].all(), sample_vs_expect.query('equal != True')

# barcode runs with R1 files expanded by glob
barcode_runs_expandR1 = (
    barcode_runs
    .assign(R1=lambda x: x['R1'].str.split('; ').map(
                    lambda y: list(itertools.chain(*map(glob.glob, y)))),
            n_R1=lambda x: x['R1'].map(len),
            sample_lib=lambda x: x['sample'] + '_' + x['library'],
            )
    )
assert barcode_runs_expandR1['sample_lib'].nunique() == len(barcode_runs_expandR1)
if any(barcode_runs_expandR1['n_R1'] < 1):
    raise ValueError(f"no R1 for {barcode_runs_expandR1.query('n_R1 < 1')}")

# Rules -----------------------------------------------------------------------

# this is the target rule (in place of `all`) since it first rule listed
rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        rulegraph=os.path.join(config['summary_dir'], 'rulegraph.svg'),
        codon_variant_table=config['codon_variant_table'],
        aggregate_variant_counts=nb_markdown('aggregate_variant_counts.ipynb'),
        variant_counts=config['variant_counts'],
        counts_to_cells_ratio=nb_markdown('counts_to_cells_ratio.ipynb'),
        counts_to_cells_csv=config['counts_to_cells_csv'],
        counts_to_scores=nb_markdown('counts_to_scores.ipynb'),
        scores_to_frac_escape=nb_markdown('scores_to_frac_escape.ipynb'),
        escape_fracs=config['escape_fracs'],
        escape_fracs_homologs=config['escape_fracs_homologs'],
        call_strong_escape_sites=nb_markdown('call_strong_escape_sites.ipynb'),
        strong_escape_sites=config['strong_escape_sites'],
        escape_profiles=nb_markdown('escape_profiles.ipynb'),
        mds_escape_profiles=nb_markdown('mds_escape_profiles.ipynb'),
        output_pdbs='results/summary/output_pdbs.md',
        gisaid_rbd_mutations=nb_markdown('gisaid_rbd_mutations.ipynb'),
        gisaid_mutation_counts=config['gisaid_mutation_counts'],
        natural_mutations=nb_markdown('natural_mutations.ipynb'),
        homolog_escape=nb_markdown('homolog_escape.ipynb'),
        escape_selections=nb_markdown('escape_selections.ipynb'),
        make_supp_data=nb_markdown('make_supp_data.ipynb'),
        structural_contacts='results/summary/annotate_structural_contacts.md',
    output:
        summary = os.path.join(config['summary_dir'], 'summary.md')
    run:
        def path(f):
            """Get path relative to `summary_dir`."""
            return os.path.relpath(f, config['summary_dir'])
        with open(output.summary, 'w') as f:
            f.write(textwrap.dedent(f"""
            # Summary

            Analysis run by [Snakefile]({path(workflow.snakefile)})
            using [this config file]({path(workflow.configfiles[0])}).
            See the [README in the top directory]({path('README.md')})
            for details.

            Here is the rule graph of the computational workflow:
            ![{path(input.rulegraph)}]({path(input.rulegraph)})

            Here is the Markdown output of each notebook in the workflow:

            1. Get codon-variant-table from [here]({config['codon_variant_table_url']}).

            2. Count variants and then
               [aggregate counts]({path(input.aggregate_variant_counts)}) to create
               to create [variant counts file]({path(input.variant_counts)}).

            3. [Analyze sequencing counts to cells ratio]({path(input.counts_to_cells_ratio)});
               this prints a list of any samples where this ratio too low. Also
               creates [a CSV]({path(input.counts_to_cells_csv)}) with the
               sequencing counts, number of sorted cells, and ratios for
               all samples.

            4. [Escape scores from variant counts]({path(input.counts_to_scores)}).

            5. [Escape fractions for mutations and homologs]({path(input.scores_to_frac_escape)});
               creating [mutation escape fraction file]({path(input.escape_fracs)})
               and [homolog escape fraction file]({path(input.escape_fracs_homologs)}).

            6. [Call sites of strong escape]({path(input.call_strong_escape_sites)}),
               and write to [a CSV file]({path(input.strong_escape_sites)}).

            7. Plot [escape profiles]({path(input.escape_profiles)}).

            9. [Multidimensional scaling]({path(input.mds_escape_profiles)}) on escape profiles.

            10. Map escape profiles to ``*.pdb`` files using [this notebook]({path(input.output_pdbs)})
            
            11. Output a list of RBD structural contacts for each antibody with a high-resolution structure using [this notebook]({path(input.structural_contacts)})

            12. [Count mutations in GISAID RBD sequences]({path(input.gisaid_rbd_mutations)}).
               to create [this counts file]({path(input.gisaid_mutation_counts)}).

            13. [Analyze GISAID mutations at sites of escape]({path(input.natural_mutations)}).

            14. [Escape from binding by homologs]({path(input.homolog_escape)}).

            15. [Analyze viral escape selections]({path(input.escape_selections)}).

            16. [Make supplementary data files]({path(input.make_supp_data)}),
                which are [here]({path(config['supp_data_dir'])}). These include
                `dms-view` input files.

            """
            ).strip())

rule make_rulegraph:
    # error message, but works: https://github.com/sequana/sequana/issues/115
    input:
        workflow.snakefile
    output:
        os.path.join(config['summary_dir'], 'rulegraph.svg')
    shell:
        "snakemake --forceall --rulegraph | dot -Tsvg > {output}"

rule structural_contacts:
    input:
        config['structural_contacts_config']
    output:
        config['structural_contacts'],
        md='results/summary/annotate_structural_contacts.md'
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='annotate_structural_contacts.Rmd',
        md='annotate_structural_contacts.md'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md}
        """

rule make_supp_data:
    input:
        config['escape_profiles_config'],
        config['output_pdbs_config'],
        config['escape_fracs'],
        config['escape_profiles_dms_colors']
    output:
        nb_markdown=nb_markdown('make_supp_data.ipynb'),
        outdir=directory(config['supp_data_dir']),
    params:
        nb='make_supp_data.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule output_pdbs:
    input:
        config['escape_fracs'],
        config['output_pdbs_config'],
    output:
        nb_markdown=nb_markdown('output_pdbs.ipynb'),
        outdir=directory(config['pdb_outputs_dir']),
    params:
        nb='output_pdbs.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule escape_selections:
    input:
        config['escape_selection_results'],
        config['escape_fracs'],
    output:
        nb_markdown=nb_markdown('escape_selections.ipynb')
    params:
        nb='escape_selections.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule homolog_escape:
    input:
        config['escape_fracs_homologs'],
        config['escape_profiles_config'],
    output:
        nb_markdown=nb_markdown('homolog_escape.ipynb')
    params:
        nb='homolog_escape.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule natural_mutations:
    input:
        config['gisaid_mutation_counts'],
        config['strong_escape_sites'],
        config['escape_fracs'],
        config['escape_profiles_config'],
    output:
        nb_markdown=nb_markdown('natural_mutations.ipynb')
    params:
        nb='natural_mutations.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule gisaid_rbd_mutations:
    input:
        config['gisaid_spikes'],
        config['wildtype_sequence'],
    output:
        config['gisaid_mutation_counts'],
        nb_markdown=nb_markdown('gisaid_rbd_mutations.ipynb'),
    params:
        nb='gisaid_rbd_mutations.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule mds_escape_profiles:
    """Multi-dimensional scaling on antibody escape profiles."""
    input:
        escape_fracs=config['escape_fracs'],
        mds_config=config['mds_config'],
        site_color_schemes=config['site_color_schemes'],
    output:
        nb_markdown=nb_markdown('mds_escape_profiles.ipynb'),
    params:
        nb='mds_escape_profiles.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule escape_profiles:
    """Make stacked logo plots of antibody escape profiles."""
    input:
        escape_fracs=config['escape_fracs'],
        escape_profiles_config=config['escape_profiles_config'],
        site_color_schemes=config['site_color_schemes'],
        wildtype_sequence=config['wildtype_sequence'],
        mut_bind_expr=config['mut_bind_expr'],
        strong_escape_sites=config['strong_escape_sites'],
    output:
        nb_markdown=nb_markdown('escape_profiles.ipynb'),
        escape_profiles_dms_colors=config['escape_profiles_dms_colors'],
    params:
        nb='escape_profiles.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule call_strong_escape_sites:
    """Call sites of strong escape."""
    input:
        escape_fracs=config['escape_fracs'],
    output:
        nb_markdown=nb_markdown('call_strong_escape_sites.ipynb'),
        strong_escape_sites=config['strong_escape_sites'],
    params:
        nb='call_strong_escape_sites.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule scores_to_frac_escape:
    """Estimate mutation- and homolog-level escape scores."""
    input:
        escape_score_samples=config['escape_score_samples'],
        escape_scores=config['escape_scores'],
        escape_scores_homologs=config['escape_scores_homologs'],
    output:
        nb_markdown=nb_markdown('scores_to_frac_escape.ipynb'),
        escape_fracs=config['escape_fracs'],
        escape_fracs_homologs=config['escape_fracs_homologs'],
    params:
        nb='scores_to_frac_escape.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule counts_to_scores:
    """Analyze variant counts to compute escape scores."""
    input:
        config['variant_counts'],
        config['wildtype_sequence'],
        config['mut_bind_expr'],
        config['variant_expr'],
        config['variant_bind'],
    output:
        nb_markdown=nb_markdown('counts_to_scores.ipynb'),
        escape_scores=config['escape_scores'],
        escape_scores_homologs=config['escape_scores_homologs'],
        escape_score_samples=config['escape_score_samples'],
    params:
        nb='counts_to_scores.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule counts_to_cells_ratio:
    input:
        config['variant_counts'],
        config['barcode_runs'],
        config['wildtype_sequence'],
    output:
        nb_markdown=nb_markdown('counts_to_cells_ratio.ipynb'),
        counts_to_cells_csv=config['counts_to_cells_csv'],
    params:
        nb='counts_to_cells_ratio.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule aggregate_variant_counts:
    input:
        counts=expand(os.path.join(config['counts_dir'],
                                   "{sample_lib}_counts.csv"),
                      sample_lib=barcode_runs_expandR1['sample_lib']),
        fates=expand(os.path.join(config['counts_dir'],
                                  "{sample_lib}_fates.csv"),
                     sample_lib=barcode_runs_expandR1['sample_lib']),
        variant_table=config['codon_variant_table'],
        wt_seq=config['wildtype_sequence'],
        barcode_runs=config['barcode_runs'],
    output:
        config['variant_counts'],
        nb_markdown=nb_markdown('aggregate_variant_counts.ipynb')
    params:
        nb='aggregate_variant_counts.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule count_variants:
    """Count variants for a specific sample."""
    input:
        variant_table=config['codon_variant_table'],
        wt_seq=config['wildtype_sequence'],
        r1s=lambda wildcards: (barcode_runs_expandR1
                               .set_index('sample_lib')
                               .at[wildcards.sample_lib, 'R1']
                               ),
    output:
        counts=os.path.join(config['counts_dir'], "{sample_lib}_counts.csv"),
        fates=os.path.join(config['counts_dir'], "{sample_lib}_fates.csv"),
    params:
        sample_lib="{sample_lib}"
    run:
        # parse sample and library from `sample_lib` wildcard
        lib = params.sample_lib.split('_')[-1]
        sample = params.sample_lib[: -len(lib) - 1]
        assert sample == (barcode_runs_expandR1
                          .set_index('sample_lib')
                          .at[params.sample_lib, 'sample']
                          )
        assert lib == (barcode_runs_expandR1
                       .set_index('sample_lib')
                       .at[params.sample_lib, 'library']
                       )
        # initialize `CodonVariantTable` (used to get valid barcodes)
        wt_seqrecord = Bio.SeqIO.read(input.wt_seq, 'fasta')
        geneseq = str(wt_seqrecord.seq)
        primary_target = wt_seqrecord.name
        variants=dms_variants.codonvarianttable.CodonVariantTable(
                    geneseq=geneseq,
                    barcode_variant_file=input.variant_table,
                    substitutions_are_codon=True,
                    substitutions_col='codon_substitutions',
                    primary_target=primary_target)
        # initialize `IlluminaBarcodeParser`
        parser = dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    valid_barcodes=variants.valid_barcodes(lib),
                    **config['illumina_barcode_parser_params'])
        # parse barcodes
        counts, fates = parser.parse(input.r1s,
                                     add_cols={'library': lib,
                                               'sample': sample})
        # write files
        counts.to_csv(output.counts, index=False)
        fates.to_csv(output.fates, index=False)

rule get_codon_variant_table:
    """Download codon variant table from URL."""
    output:
        codon_variant_table=config['codon_variant_table']
    run:
        urllib.request.urlretrieve(config['codon_variant_table_url'],
                                   output.codon_variant_table)

rule get_mut_bind_expr:
    """Download SARS-CoV-2 mutation ACE2-binding and expression from URL."""
    output:
        file=config['mut_bind_expr']
    run:
        urllib.request.urlretrieve(config['mut_bind_expr_url'], output.file)

rule get_variant_expr:
    """Download SARS-CoV-2 variant expression from URL."""
    output:
        file=config['variant_expr']
    run:
        urllib.request.urlretrieve(config['variant_expr_url'], output.file)

rule get_variant_bind:
    """Download SARS-CoV-2 variant ACE2-binding from URL."""
    output:
        file=config['variant_bind']
    run:
        urllib.request.urlretrieve(config['variant_bind_url'], output.file)
