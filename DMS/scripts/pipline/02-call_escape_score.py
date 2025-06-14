import pandas as pd
from scipy.stats import binom_test, poisson
import numpy as np
import argparse
import sys, os
import Bio.SeqIO
import time
import dms_variants.globalepistasis
import dms_variants.codonvarianttable
import dms_variants.binarymap
import collections, itertools
from lets_plot import GGBunch, ylim, ggplot, ggsave, aes, geom_boxplot, geom_point, geom_text, scale_x_discrete, ggtitle

parser = argparse.ArgumentParser(description="Step 2: calculate escape score for every barcode.")

parser.add_argument('--reference_variant_counts', '-r', type=str)
parser.add_argument('--escape_variant_counts', '-e', type=str)
parser.add_argument('--experiment', '-exp', type=str, default='MACS')
parser.add_argument('--wildtype', '-w', type=str)
parser.add_argument('--table', '-t', type=str)
parser.add_argument('--frac_escape', type=float, default = None)
parser.add_argument('--cells_sorted', type=int, default = None)
parser.add_argument('--output_dir', '-o', type=str)
parser.add_argument('--primary_target', '-p', type=str, default='SARS-CoV-2')
parser.add_argument('--wtseqstart', '-S', type=int, default=330)
parser.add_argument('--mutbindexpr', type=str)
parser.add_argument('--variant_bind', type=str)
parser.add_argument('--variant_expr', type=str)
parser.add_argument('--usevarfilt', action="store_true")
parser.add_argument('--afterfitfilter', action="store_true")
parser.add_argument('--exprmin', type=float, default=-1.0)
parser.add_argument('--bindmin', type=float, default=-2.35)
parser.add_argument("--mincount", type=int, default=5)
parser.add_argument("--min_var_support_call", type=int, default=1)
parser.add_argument("--protchain", type=str, default='E')

def build_variant_table(table_path, wt_seq_path):
    wt_seqrecord = Bio.SeqIO.read(wt_seq_path, 'fasta')
    geneseq = str(wt_seqrecord.seq)
    primary_target = wt_seqrecord.name
    variants = dms_variants.codonvarianttable.CodonVariantTable(
               geneseq = geneseq,
               barcode_variant_file = table_path,
               substitutions_are_codon = True,
               substitutions_col = 'codon_substitutions',
               primary_target = primary_target)
    return variants, primary_target


def print_some_logs(args):
    sys.stdout.write('Reference: ' + args.reference_variant_counts + '\n')
    sys.stdout.write('Escape: ' + args.escape_variant_counts + '\n')
    sys.stdout.write('Experiment: ' + args.experiment + '\t')
    if args.experiment == 'FACS':
        sys.stdout.write('FACS result: ' + str(args.cells_sorted) + ' cells , and ' + str(args.frac_escape) + ' escaped\n')

def filter_expr_bind(escape_scores, mutbindexpr, variant_bind, variant_expr, min_expr_mut, min_expr_variant, min_bind_mut, min_bind_variant, filter_on_variants=False):
    print(mutbindexpr, variant_expr, variant_bind, filter_on_variants, min_expr_mut)
    # filter on mutations
    if mutbindexpr is not None:
        mut_bind_expr = pd.read_csv(mutbindexpr)
        assert mut_bind_expr['mutation_RBD'].nunique() == len(mut_bind_expr)
        for prop in ['bind', 'expr']:
            muts_adequate = set(mut_bind_expr
                                .query(f"{prop}_avg >= {(min_expr_mut if prop == 'expr' else min_bind_mut)}")
                                ['mutation_RBD']
                                )
            sys.stdout.write(str(len(muts_adequate))+" of "+str(len(mut_bind_expr))+" mutations have adequate "+prop+"\n")
            escape_scores[f"muts_pass_{prop}_filter"] = (
                escape_scores
                ['aa_substitutions']
                .map(lambda s: set(s.split()).issubset(muts_adequate))
                ) 
    else:
        sys.stdout.write("Skip mut_bind_expr filter\n")
        escape_scores["muts_pass_bind_filter"] = True
        escape_scores["muts_pass_expr_filter"] = True
    
    if filter_on_variants and os.path.exists(variant_bind):
        # filter on variants
        sys.stdout.write("Apply variant bind filter")
        col = 'delta_log10Ka'
        filter_name = "variant_pass_bind_filter"
        variant_pass_df = (
            pd.read_csv(variant_bind, keep_default_na=False, na_values=['NA'])
            .groupby(['library', 'target', 'barcode'])
            .aggregate(val=pd.NamedAgg(col, 'mean'))
            .reset_index()
            .assign(pass_filter=lambda x: x['val'] >= min_bind_variant)
            .rename(columns={'pass_filter': filter_name,
                            'val': "variant_bind"})
            )
        escape_scores = (
            escape_scores
            .drop(columns=filter_name, errors='ignore')
            .merge(variant_pass_df,
                how='left',
                validate='many_to_one',
                on=['library', 'target', 'barcode'],
                )
            )
        assert escape_scores[filter_name].notnull().all()
    else:
        sys.stdout.write("Skip variant bind filter\n")
        escape_scores["variant_pass_bind_filter"]=True


    if filter_on_variants and os.path.exists(variant_expr):
        # filter on variants
        sys.stdout.write("Apply variant expr filter")
        col = 'delta_ML_meanF'
        filter_name = "variant_pass_expr_filter"
        variant_pass_df = (
            pd.read_csv(variant_expr, keep_default_na=False, na_values=['NA'])
            .groupby(['library', 'target', 'barcode'])
            .aggregate(val=pd.NamedAgg(col, 'mean'))
            .reset_index()
            .assign(pass_filter=lambda x: x['val'] >= min_expr_variant)
            .rename(columns={'pass_filter': filter_name,
                            'val': "variant_expr"})
            )
        escape_scores = (
            escape_scores
            .drop(columns=filter_name, errors='ignore')
            .merge(variant_pass_df,
                how='left',
                validate='many_to_one',
                on=['library', 'target', 'barcode'],
                )
            )
        assert escape_scores[filter_name].notnull().all()
    else:
        sys.stdout.write("Skip variant expr filter\n")
        escape_scores["variant_pass_expr_filter"]=True
    
    # annotate as passing overall filter if passes all mutation and binding filters:
    escape_scores['pass_ACE2bind_expr_filter'] = (
            escape_scores['muts_pass_bind_filter'] &
            escape_scores['muts_pass_expr_filter'] &
            escape_scores['variant_pass_bind_filter'] &
            escape_scores['variant_pass_expr_filter']
            )
    
    return escape_scores.query('pass_ACE2bind_expr_filter == True')


def calc_epistatsis_model_and_output(df, suffix, lib, escape_name, ref_name, wtseqstart, name, output_dir, bunch, protein_chain = 'E', plot_max=1.0):
    sys.stdout.write(f'Fitting epistasis model {suffix} ...\n')
    binary_map = dms_variants.binarymap.BinaryMap(df, func_score_col='escape_score')
    model = dms_variants.globalepistasis.MonotonicSplineEpistasisGaussianLikelihood(binary_map)
    model.fit()
    sys.stdout.write('Done.\n')
    variant_counter = collections.Counter(model.binarymap.substitution_variants)
    muts_counter = collections.Counter(itertools.chain.from_iterable(s.split() for s in model.binarymap.substitution_variants))
    
    effects_df = (pd.DataFrame({'condition': escape_name+'_over_'+ref_name,
                    'library': lib,
                    'mutation': model.binarymap.all_subs})
        .assign(n_single_mut_measurements=lambda x: x['mutation'].map(variant_counter),
                n_any_mut_measurements=lambda x: x['mutation'].map(muts_counter),
                )
        .pipe(model.add_phenotypes_to_df, substitutions_col='mutation')
        .drop(columns='latent_phenotype')
        .rename(columns={'observed_phenotype': 'epistasis_model_score'})
        )
    sys.stdout.write("Removing mutations that do not have EITHER >= 1 single-mutant measurements or >= 3 any-mutant measurements.\n")
    effects_df = (
        effects_df
        .assign(sufficient_measurements=lambda x: (
                    (x['n_single_mut_measurements'] >= 1) |
                    (x['n_any_mut_measurements'] >= 3))
                )
        .query('sufficient_measurements == True')
        .drop(columns='sufficient_measurements')
    )
    raw_avg_single_mut_scores = (
        df.query('n_aa_substitutions == 1')
        .rename(columns={'aa_substitutions': 'mutation'})
        .groupby(['condition', 'library', 'mutation'])
        .aggregate(raw_single_mut_score=pd.NamedAgg('escape_score', 'mean'))
        .reset_index()
        )
    
    mn = max(0,np.nanquantile(effects_df['epistasis_model_score'], 0.01)) # 0
    mx = min(1,np.nanquantile(effects_df['epistasis_model_score'], 0.996)) # 1
    if mx < 10 * mn:
        mx = max(effects_df['epistasis_model_score'])
        if mx < 10 * mn:
            mx = 1
    print(mn, mx)
    effects_df = effects_df.merge(raw_avg_single_mut_scores, how='outer', validate='one_to_one').assign(site=lambda x: x['mutation'].str[1: -1].astype(int),
            wildtype=lambda x: x['mutation'].str[0],
            mutation=lambda x: x['mutation'].str[-1],
            ).assign(
                mut_escape=lambda x: [((i-mn)/(mx-mn) if (i >= mn and i <= mx) else (0 if (i < (mn+mx)/2) else 1)) for i in x['epistasis_model_score']], #x['epistasis_model_score'],
                # mut_escape=lambda x: [(i if (i >= mn and i <= mx) else (mn if (i < (mn+mx)/2) else mx)) for i in x['epistasis_model_score']], #x['epistasis_model_score'],
                protein_chain=protein_chain
            )

    mn = 0
    mx = 1 # np.nanquantile(effects_df['raw_single_mut_score'], 0.999)

    effects_df['site'] = effects_df['site'] + wtseqstart
    effects_df = effects_df.assign(
        single_mut_escape=lambda x: [((i-mn)/mx if (i >= mn and i <= mx) else (mn if (i < (mn+mx)/2) else (np.nan if np.isnan(i) else mx))) for i in x['raw_single_mut_score']],
        protein_site=effects_df['site'],
        label_site=effects_df['wildtype']+effects_df['site'].astype(str)
    )

    site_effects_df = (
        effects_df
        .groupby('site')
        .aggregate(
            site_avg_escape_frac_epistasis_model=pd.NamedAgg('mut_escape','mean'),
            site_total_escape_frac_epistasis_model=pd.NamedAgg('mut_escape', 'sum'),
            site_avg_escape_frac_single_mut=pd.NamedAgg('single_mut_escape', 'mean'),
            site_total_escape_frac_single_mut=pd.NamedAgg('single_mut_escape', 'sum'),
            )
        .reset_index()
        )
    site_effects_df.index = site_effects_df['site']

    effects_df = effects_df.assign(
            site_total_escape=site_effects_df['site_total_escape_frac_epistasis_model'][effects_df['site']].to_numpy(),
            site_mean_escape=site_effects_df['site_avg_escape_frac_epistasis_model'][effects_df['site']].to_numpy(),
            site_total_escape_single=site_effects_df['site_total_escape_frac_single_mut'][effects_df['site']].to_numpy(),
            site_mean_escape_single=site_effects_df['site_avg_escape_frac_single_mut'][effects_df['site']].to_numpy(),
        )

    effects_df.to_csv(os.path.join(output_dir, f"effects_df{suffix}.csv"), index=False)
    site_effects_df.to_csv(os.path.join(output_dir, f"site_effects_df{suffix}.csv"), index=False)
    site_effects_df = site_effects_df.assign(targetsite = site_effects_df['site'])
    
    plots = []

    plots.append(ggplot(site_effects_df, aes(x = 'targetsite', y = 'site_avg_escape_frac_epistasis_model')) + geom_point() + geom_text(aes(label = 'targetsite')) + scale_x_discrete(breaks=list(range(min(site_effects_df['targetsite']),max(site_effects_df['targetsite'])))) + ggtitle(name))
    plots.append(ggplot(site_effects_df, aes(x = 'targetsite', y = 'site_avg_escape_frac_single_mut')) + geom_point() + geom_text(aes(label = 'targetsite')) + scale_x_discrete(breaks=list(range(min(site_effects_df['targetsite']),max(site_effects_df['targetsite'])))))
    plots.append(ggplot(site_effects_df, aes(x = 'targetsite', y = 'site_total_escape_frac_epistasis_model')) + geom_point() + geom_text(aes(label = 'targetsite')) + scale_x_discrete(breaks=list(range(min(site_effects_df['targetsite']),max(site_effects_df['targetsite'])))) + ggtitle(name))
    plots.append(ggplot(site_effects_df, aes(x = 'targetsite', y = 'site_total_escape_frac_single_mut')) + geom_point() + geom_text(aes(label = 'targetsite')) + scale_x_discrete(breaks=list(range(min(site_effects_df['targetsite']),max(site_effects_df['targetsite'])))))
    plots.append(ggplot(site_effects_df, aes(x = 'targetsite', y = 'site_avg_escape_frac_epistasis_model')) + geom_point() + geom_text(aes(label = 'targetsite')) + scale_x_discrete(breaks=list(range(min(site_effects_df['targetsite']),max(site_effects_df['targetsite'])))) + ggtitle(name) + ylim(0,plot_max))
    plots.append(ggplot(site_effects_df, aes(x = 'targetsite', y = 'site_avg_escape_frac_single_mut')) + geom_point() + geom_text(aes(label = 'targetsite')) + scale_x_discrete(breaks=list(range(min(site_effects_df['targetsite']),max(site_effects_df['targetsite'])))) + ylim(0,plot_max))
    
    for i in range(len(plots)):
        plt = plots[i]
        bunch.add_plot(plt, 0, (i+1)*400,1500,300)

    ggsave(bunch, f'call_escape_score_plots{suffix}.html', path = output_dir)

def main():
    args = parser.parse_args()
    bunch = GGBunch()
    protein_chain = args.protchain
    
    assert (args.experiment == 'MACS') or (args.experiment == 'FACS' and (args.frac_escape is not None) and (args.cells_sorted is not None)), 'experiment argument error\n'

    ref = pd.read_csv(args.reference_variant_counts).query('variant_call_support >= @args.min_var_support_call')
    escape = pd.read_csv(args.escape_variant_counts).query('variant_call_support >= @args.min_var_support_call')

    ref_name = '_'.join(args.reference_variant_counts.split('/')[-1].split('_')[0:-2])
    escape_name = '_'.join(args.escape_variant_counts.split('/')[-1].split('_')[0:-2])

    lib = '2mutBA5Tmerged'

    primary_target = args.primary_target
    output_dir = os.path.join(args.output_dir, escape_name+'_over_'+ref_name)

    os.makedirs(output_dir, exist_ok=True)
    print_some_logs(args)

    ncounts_ref = np.sum(ref['count'].to_numpy())
    ncounts_escape = np.sum(escape['count'].to_numpy())
    frac_escape = args.frac_escape
    if frac_escape is None:
        frac_escape = 1.0

    variants, primary_target = build_variant_table(args.table, args.wildtype)

    ref = (ref.merge(escape[['sample','barcode','count']], on = ['barcode'])
              .rename(columns = {'count_x':'reference_count','count_y':'escape_count'}))
    ref = ref.assign(escape_score = lambda x: (x.escape_count/ncounts_escape) / (x.reference_count/ncounts_ref) * frac_escape)
    ref = ref[(ref['escape_score'] != np.inf) & (ref['escape_score'] >= 0) & (ref['reference_count'] > args.mincount)].sort_values('escape_score')
    
    # num_use = int(len(ref)*0.99)
    # ref = ref[0:num_use]
    df = ref.assign(aa_sub_type = lambda x: (['>1' if i >1 else str(i) for i in x.n_aa_substitutions]))
    df = df[df['target'] == primary_target].fillna('')
    df = (df.assign(norm_score = (df['escape_score'])/np.std(df['escape_score']), condition = escape_name+'_over_'+ref_name)
            .drop(columns = ['sample_x','sample_y'])
            .pipe(variants.classifyVariants,
                        primary_target=primary_target,
                        syn_as_wt=False))

    mn = np.nanquantile(df['escape_score'], 0.01)
    mx = np.nanquantile(df['escape_score'], 0.99)
    if mx == 0:
        _ = -3
        while mx == 0:
            mx = np.nanquantile(df['escape_score'], 1.0-(10**_))
            _ -= 1

    print(mn, mx)

    df = df.assign(raw_escape_score = df['escape_score'])
    
    df['escape_score'] = [(i if (i > mn and i < mx) else (0 if (i < (mn+mx)/2) else mx))/mx for i in df['raw_escape_score']]

    p1 = ggplot(df, aes(x = 'variant_class', y = 'escape_score')) + geom_boxplot()
    bunch.add_plot(p1, 0, 0)

    
    # For RBD libraries don't apply this filter
    args.variant_expr = None
    args.variant_bind = None
    
    df_filter = filter_expr_bind(df, args.mutbindexpr, args.variant_bind, args.variant_expr, args.exprmin, args.exprmin, args.bindmin, args.bindmin, args.usevarfilt).query('variant_class != "stop"')
    df = df.query('variant_class != "stop"')
            
    calc_epistatsis_model_and_output(df_filter, '', lib, escape_name, ref_name, args.wtseqstart, args.escape_variant_counts, output_dir, bunch)
    calc_epistatsis_model_and_output(df, '_no_filter', lib, escape_name, ref_name, args.wtseqstart, args.escape_variant_counts, output_dir, bunch)

if __name__ == '__main__':
    main()
