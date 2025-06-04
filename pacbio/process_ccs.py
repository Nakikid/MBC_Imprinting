# Code from J. Bloom et al.
# https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS/
# Modified by Jian F. 2021/10/19

import alignparse
import alignparse.ccs
import alignparse.minimap2
import alignparse.targets
import alignparse.consensus
from alignparse.constants import CBPALETTE

import matplotlib
matplotlib.use('Agg')

import dms_variants
import dms_variants.utils
import dms_variants.codonvarianttable
import dms_variants.plotnine_themes

import numpy
import pandas as pd
import argparse

import plotnine
from plotnine import theme, element_blank, ggplot, facet_wrap, geom_point, scale_y_log10, scale_x_log10, geom_text, position_dodge, facet_grid, scale_x_continuous
from plotnine import aes, scale_fill_manual, element_text, geom_bar, geom_histogram, geom_vline, ylab, xlab, scale_color_manual,scale_y_continuous
from plotnine.ggplot import save_as_pdf_pages
import math
import os,sys

parser = argparse.ArgumentParser()
parser.add_argument('--seqsfile', type=str, default="pacbio_amplicons.gb")
parser.add_argument('--configfile', type=str, default="parse_specs.yaml")
parser.add_argument('--threads', '-j', type=int)
parser.add_argument('--pbruns', type=str, default="pacbio_runs.csv")
parser.add_argument('--output', '-o', type=str)
parser.add_argument('--target', type=str)

args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)
pdf = os.path.join(args.output, 'process_ccs_plots.pdf') 
plots = []
processed_ccs = os.path.join(args.output, 'processed_ccs.csv')

targets = alignparse.targets.Targets(seqsfile=args.seqsfile,
                                     feature_parse_specs=args.configfile)

pacbio_runs = pd.read_csv(args.pbruns)
ccs_summaries = alignparse.ccs.Summaries(pacbio_runs, 
                                        #  report_col = 'report',
                                        report_col = None,
                                        ncpus = args.threads)
if ccs_summaries.has_zmw_stats():
    p = ccs_summaries.plot_zmw_stats()
    p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
    plots.append(p)
else:
    print('No ZMW stats available.')

for variable in ['length', 'passes', 'accuracy']:
    if ccs_summaries.has_stat(variable):
        p = ccs_summaries.plot_ccs_stats(variable, maxcol=7, bins=25)
        p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
        plots.append(p)
    else:
        print(f"No {variable} statistics available.")

options = alignparse.minimap2.OPTIONS_CODON_DMS
options = list(options)
options[4] = '--end-bonus=23'
options = tuple(options)
mapper = alignparse.minimap2.Mapper(options)

readstats, aligned, filtered = targets.align_and_parse(
    df = pacbio_runs,
    mapper = mapper,
    outdir = args.output,
    name_col = 'run',
    group_cols = ['library', 'name'],
    queryfile_col = 'fastq',
    overwrite = True,
    ncpus = args.threads)

######################

aligned_df = (
    pd.concat(df.assign(target=target) for target, df in aligned.items())
    .drop(columns=['query_clip5', 'query_clip3', 'run','name'])
    .rename(columns={'barcode_sequence': 'barcode'})
    )

aligned_df.to_csv(processed_ccs, index=False)

##########################
##### build variants #####
##########################

pdf = os.path.join(args.output, 'build_variants_plots.pdf') 
plots = []

processed_ccs = aligned_df
nlibs = processed_ccs['library'].nunique()
ntargets = processed_ccs['target'].nunique()

print(f"Read {len(processed_ccs)} CCSs from {nlibs} libraries and {ntargets} targets.")

error_rate_floor = 1e-7
max_error_rate = 1e-4

processed_ccs = (
    processed_ccs
    .assign(barcode_error=lambda x: numpy.clip(1 - x['barcode_accuracy'],
                                               error_rate_floor, None),
            gene_error=lambda x: numpy.clip(1 - x['gene_accuracy'],
                                            error_rate_floor, None)
            ).assign(retained=lambda x: ((x['gene_error'] < max_error_rate) &
                                (x['barcode_error'] < max_error_rate)))
    )



processed_ccs = alignparse.consensus.add_mut_info_cols(processed_ccs,
                                                       mutation_col='gene_mutations',
                                                       n_indel_col='n_indels')

processed_ccs = processed_ccs.assign(has_indel=lambda x: x['n_indels'] > 0)
high_acc = max_error_rate / 10
empirical_acc = []

for desc, query_str in [
        ('retained', 'retained'),
        ('retained, no indel', 'retained and not has_indel'),
        ('10X accuracy',
         f"(gene_error < {high_acc}) and (barcode_error < {high_acc})"),
        ('10X accuracy, no indel',
         f"(gene_error < {high_acc}) and (barcode_error < {high_acc}) and not has_indel")
        ]:
    # get just CCSs in that category
    df = processed_ccs.query(query_str)
    
    # compute empirical accuracy
    empirical_acc.append(
        alignparse.consensus.empirical_accuracy(df,
                                                mutation_col='gene_mutations')
        .assign(description=desc)
        .merge(df
               .groupby('library')
               .size()
               .rename('number_CCSs')
               .reset_index()
               )
        )

# make description categorical to preserve order, and annotate as "actual"
# the category ("retained, no indel") that we will use for building variants.
empirical_acc = (
    pd.concat(empirical_acc, ignore_index=True, sort=False)
    .assign(description=lambda x: pd.Categorical(x['description'],
                                                 x['description'].unique(),
                                                 ordered=True),
            actual=lambda x: numpy.where(x['description'] == 'retained, no indel',
                                         True, False),
            )
    )

consensus, dropped = alignparse.consensus.simple_mutconsensus(
                        processed_ccs.query('retained'),
                        group_cols=('library', 'barcode', 'target'),
                        mutation_col='gene_mutations',
                        )
consensus = alignparse.consensus.add_mut_info_cols(
                    consensus,
                    mutation_col='gene_mutations',
                    sub_str_col='substitutions',
                    n_indel_col='number_of_indels',
                    overwrite_cols=True)



consensus = consensus.query('number_of_indels == 0')

lib_target_counts = (
    consensus
    .groupby(['library', 'target'])
    .size()
    .rename('consensus sequences')
    .reset_index()
    )
p = (ggplot(lib_target_counts.assign(xlabel=lambda x: x['target'] + ', ' + x['library']),
            aes('xlabel', 'consensus sequences')) +
     geom_point(size=3) +
     theme(figure_size=(0.5 * nlibs * ntargets, 1.75),
           axis_text_x=element_text(angle=90)) +
     xlab('') +
     scale_y_log10()
     )
plots.append(p)

primary_target = args.target
print(f"Dropping variants with mutations for all targets except {primary_target}")

consensus = (
    consensus
    .assign(has_substitutions=lambda x: x['substitutions'].str.len().astype(bool))
    )

has_subs_by_target = (
        consensus
        .groupby(['target', 'library', 'has_substitutions'])
        .aggregate(n_barcodes=pd.NamedAgg('barcode', 'count'))
        .reset_index()
        )


nt_variant_table_file = os.path.join(args.output, 'nt_variant_table.csv')

print(f"Culling the {len(consensus)} barcodes to remove mutated non-primary targets")

consensus = consensus.query('(target == @primary_target) or (has_substitutions == False)')

print(f"Retained {len(consensus)} barcodes after culling")
dup_barcodes = (
    consensus
    .groupby(['library', 'barcode'])
    .size()
    .rename('duplicate_count')
    .reset_index()
    .query('duplicate_count > 1')
    )

print(dup_barcodes.head())
print(f"\nRemoving the {len(dup_barcodes)} duplicated barcodes."
      f"Started with {len(consensus)} barcodes:")

consensus = (
    consensus
    .merge(dup_barcodes, on=['library', 'barcode'], how='outer')
    .query('duplicate_count.isnull()', engine='python')
    )
print(f"After removing duplicates, there are {len(consensus)} barcodes.")
print(f"Writing nucleotide variants to {nt_variant_table_file}")

consensus[['target', 'library', 'barcode', 'substitutions', 'variant_call_support']].to_csv(nt_variant_table_file, index=False)

print(dropped.head())

max_nseqs = 8  # plot together all barcodes with >= this many sequences

p = (
 ggplot(
    dropped.assign(nseqs=lambda x: numpy.clip(x['nseqs'], None, max_nseqs)),
    aes('nseqs')) + 
 geom_bar() + 
 scale_x_continuous(limits=(1, None)) +
 xlab('number of sequences for barcode') +
 ylab('number of barcodes') +
 facet_grid('library ~ drop_reason') +
 theme(figure_size=(10, 1.5 * nlibs),
       panel_grid_major_x=element_blank(),
       )
 )
plots.append(p)

geneseq = targets.get_target(primary_target).get_feature('gene').seq

print(f"Read gene of {len(geneseq)} nts for {primary_target} from {args.seqsfile}")

variants = dms_variants.codonvarianttable.CodonVariantTable(
                barcode_variant_file=nt_variant_table_file,
                geneseq=geneseq,
                primary_target=primary_target,
                )

print(variants.n_variants_df(samples=None).pivot_table(index=['target'], columns='library', values='count').head())

max_support = 10  # group variants with >= this much support
p = variants.plotVariantSupportHistogram(max_support=max_support,
                                         widthscale=1.1,
                                         heightscale=0.9)
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
plots.append(p)
max_muts = 7  # group all variants with >= this many mutations
for mut_type in ['aa', 'codon']:
    p = variants.plotNumMutsHistogram(mut_type, samples=None, max_muts=max_muts,
                                      widthscale=1.1,
                                      heightscale=0.9)
    p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
    plots.append(p)

p = variants.plotNumCodonMutsByType(variant_type='all', samples=None,
                                    ylabel='mutations per variant',
                                    heightscale=0.8)
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
plots.append(p)
p = variants.plotNumCodonMutsByType(variant_type='all', samples=None,
                                    ylabel='mutations per variant', 
                                    min_support=2, heightscale=0.8)
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
plots.append(p)

for mut_type in ['aa', 'codon']:
    p = variants.plotMutHeatmap('all', mut_type, samples=None, #libraries='all_only',
                                widthscale=2)
    plots.append(p)


#save_as_pdf_pages(plots, filename=pdf)
######################
codon_variant_table_file = os.path.join(args.output, 'codon_variant_table.csv')
print(f"Writing codon-variant table to {codon_variant_table_file}")

variants.barcode_variant_df.to_csv(codon_variant_table_file, index=False)