"""A trimed down distillation scripts for use with CAMPERS"""
import os
from os import path, mkdir
import click
import pandas as pd
from collections import Counter, defaultdict
import altair as alt
import networkx as nx
from itertools import tee
import re
import numpy as np
from datetime import datetime
from camper_dramkit.camper_annotate import CAMPER_NAME


# TODO: add RBH information to output
# TODO: add flag to output table and not xlsx
# TODO: add flag to output heatmap table

FRAME_COLUMNS = ['gene_id', 'gene_description', 'module', 'sheet', 'header', 'subheader']
RRNA_TYPES = ['5S rRNA', '16S rRNA', '23S rRNA']
HEATMAP_MODULES = ['M00001', 'M00004', 'M00008', 'M00009', 'M00012', 'M00165', 'M00173', 'M00374', 'M00375',
                   'M00376', 'M00377', 'M00422', 'M00567']
HEATMAP_CELL_HEIGHT = 10
HEATMAP_CELL_WIDTH = 10
KO_REGEX = r'^K\d\d\d\d\d$'
ETC_COVERAGE_COLUMNS = ['module_id', 'module_name', 'complex', 'genome', 'path_length', 'path_length_coverage',
                        'percent_coverage', 'genes', 'missing_genes', 'complex_module_name']
TAXONOMY_LEVELS = ['d', 'p', 'c', 'o', 'f', 'g', 's']

DEFAULT_CAMPER_DIST = os.path.join(os.path.dirname(__file__) , "data", "CAMPER_distillate.tsv")


# UTILS
def get_ordered_uniques(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def get_ids_from_row(row):
    id_list = list()
    # get kegg gene ids
    if 'kegg_genes_id' in row and not pd.isna(row['kegg_genes_id']):
        id_list += row['kegg_genes_id']
    # get kegg orthology ids
    if 'ko_id' in row and not pd.isna(row['ko_id']):
        id_list += [j for j in row['ko_id'].split(',')]
    # get ec numbers
    if 'kegg_hit' in row and not pd.isna(row['kegg_hit']):
        id_list += [i[1:-1] for i in re.findall(r'\[EC:\d*.\d*.\d*.\d*\]', row['kegg_hit'])]
    # get merops ids
    if 'peptidase_family' in row and not pd.isna(row['peptidase_family']):
        id_list += [j for j in row['peptidase_family'].split(';')]
    # get cazy ids
    if 'cazy_hits' in row and not pd.isna(row['cazy_hits']):
        id_list += [i[1:-1].split('_')[0] for i in re.findall(r'\[[A-Z]*\d*?\]', row['cazy_hits'])]
    # get pfam ids
    if 'pfam_hits' in row and not pd.isna(row['pfam_hits']):
        id_list += [j[1:-1].split('.')[0]
                    for j in re.findall(r'\[PF\d\d\d\d\d.\d*\]', row['pfam_hits'])]
    # custom campers id
    if f"{CAMPER_NAME}_id" in row:
        id_list += [row[f"{CAMPER_NAME}_id"]]
    return set(id_list)


def get_ids_from_annotation(frame):
    id_list = list()
    # get kegg gene ids
    if 'kegg_genes_id' in frame:
        id_list += [j.strip() for i in frame.kegg_genes_id.dropna() for j in i.split(',')]
    # get kegg orthology ids
    if 'ko_id' in frame:
        id_list += [j.strip() for i in frame.ko_id.dropna() for j in i.split(',')]
    # get kegg ec numbers
    if 'kegg_hit' in frame:
        for kegg_hit in frame.kegg_hit.dropna():
            id_list += [i[1:-1] for i in re.findall(r'\[EC:\d*.\d*.\d*.\d*\]', kegg_hit)]
    # get merops ids
    if 'peptidase_family' in frame:
        id_list += [j.strip() for i in frame.peptidase_family.dropna() for j in i.split(';')]
    # get cazy ids
    if 'cazy_id' in frame:
        id_list += [j for i in frame.cazy_id.dropna() for j in i.split('; ')]
    # get cazy ec numbers
    if 'cazy_hits' in frame:
        id_list += [f"{j[1:3]}:{j[4:-1]}" for i in frame.cazy_hits.dropna()
                    for j in re.findall(r'\(EC [\d+\.]+[\d-]\)', i)]
    # get pfam ids
    if 'pfam_hits' in frame:
        id_list += [j[1:-1].split('.')[0] for i in frame.pfam_hits.dropna()
                    for j in re.findall(r'\[PF\d\d\d\d\d.\d*\]', i)]
    if f"{CAMPER_NAME}_id" in frame:
        id_list += frame[f"{CAMPER_NAME}_id"].dropna().apply(lambda x: x.strip()).tolist()
    return Counter(id_list)


def fill_genome_summary_frame(annotations, genome_summary_frame, groupby_column):
    genome_summary_id_sets = [set([k.strip() for k in j.split(',')]) for j in genome_summary_frame['gene_id']]
    for genome, frame in annotations.groupby(groupby_column, sort=False):
        id_dict = get_ids_from_annotation(frame)
        counts = list()
        for i in genome_summary_id_sets:
            identifier_count = 0
            for j in i:
                if j in id_dict:
                    identifier_count += id_dict[j]
            counts.append(identifier_count)
        genome_summary_frame[genome] = counts
    return genome_summary_frame


def fill_genome_summary_frame_gene_names(annotations, genome_summary_frame, groupby_column):
    genome_summary_id_sets = [set([k.strip() for k in j.split(',')]) for j in genome_summary_frame['gene_id']]
    for genome, frame in annotations.groupby(groupby_column, sort=False):
        # make dict of identifiers to gene names
        id_gene_dict = defaultdict(list)
        for gene, row in frame.iterrows():
            ids = get_ids_from_row(row)
            for id_ in ids:
                id_gene_dict[id_].append(gene)
        # fill in genome summary_frame
        values = list()
        for id_set in genome_summary_id_sets:
            this_value = list()
            for id_ in id_set:
                this_value += id_gene_dict[id_]
            values.append(','.join(this_value))
        genome_summary_frame[genome] = values
    return genome_summary_frame


def make_genome_summary(annotations, genome_summary_frame, groupby_column='fasta'):
    summary_frames = list()
    # get ko summaries
    summary_frames.append(fill_genome_summary_frame(annotations, genome_summary_frame.copy(), groupby_column))

    # merge summary frames
    summarized_genomes = pd.concat(summary_frames, sort=False)
    return summarized_genomes


def write_summarized_genomes_to_xlsx(summarized_genomes, output_file):
    # turn all this into an xlsx
    with pd.ExcelWriter(output_file) as writer:
        for sheet, frame in summarized_genomes.groupby('sheet', sort=False):
            frame = frame.sort_values(['header', 'subheader', 'module', 'gene_id'])
            frame = frame.drop(['sheet'], axis=1)
            frame.to_excel(writer, sheet_name=sheet, index=False)


# TODO: add assembly stats like N50, longest contig, total assembled length etc
def make_genome_stats(annotations, rrna_frame=None, trna_frame=None, groupby_column='fasta'):
    rows = list()
    columns = ['genome']
    if 'scaffold' in annotations.columns:
        columns.append('number of scaffolds')
    if 'bin_taxonomy' in annotations.columns:
        columns.append('taxonomy')
    if 'bin_completeness' in annotations.columns:
        columns.append('completeness score')
    if 'bin_contamination' in annotations.columns:
        columns.append('contamination score')
    if rrna_frame is not None:
        columns += RRNA_TYPES
    if trna_frame is not None:
        columns.append('tRNA count')
    if 'bin_completeness' in annotations.columns and 'bin_contamination' in annotations.columns \
       and rrna_frame is not None and trna_frame is not None:
        columns.append('assembly quality')
    for genome, frame in annotations.groupby(groupby_column, sort=False):
        row = [genome]
        if 'scaffold' in frame.columns:
            row.append(len(set(frame['scaffold'])))
        if 'bin_taxonomy' in frame.columns:
            row.append(frame['bin_taxonomy'][0])
        if 'bin_completeness' in frame.columns:
            row.append(frame['bin_completeness'][0])
        if 'bin_contamination' in frame.columns:
            row.append(frame['bin_contamination'][0])
        has_rrna = list()
        if rrna_frame is not None:
            genome_rrnas = rrna_frame.loc[rrna_frame.fasta == genome]
            for rrna in RRNA_TYPES:
                sixteens = genome_rrnas.loc[genome_rrnas.type == rrna]
                if sixteens.shape[0] == 0:
                    row.append('')
                    has_rrna.append(False)
                elif sixteens.shape[0] == 1:
                    row.append('%s (%s, %s)' % (sixteens['scaffold'].iloc[0], sixteens.begin.iloc[0],
                                                sixteens.end.iloc[0]))
                    has_rrna.append(True)
                else:
                    row.append('%s present' % sixteens.shape[0])
                    has_rrna.append(False)
        if trna_frame is not None:
            row.append(trna_frame.loc[trna_frame[groupby_column] == genome].shape[0])  # TODO: remove psuedo from count?
        if 'assembly quality' in columns:
            if frame['bin_completeness'][0] > 90 and frame['bin_contamination'][0] < 5 and np.all(has_rrna) and \
               len(set(trna_frame.loc[trna_frame[groupby_column] == genome].Type)) >= 18:
                row.append('high')
            elif frame['bin_completeness'][0] >= 50 and frame['bin_contamination'][0] < 10:
                row.append('med')
            else:
                row.append('low')
        rows.append(row)
    genome_stats = pd.DataFrame(rows, columns=columns)
    return genome_stats


def build_module_net(module_df):
    """Starts with a data from including a single module"""
    # build net from a set of module paths
    num_steps = max([int(i.split(',')[0]) for i in set(module_df['path'])])
    module_net = nx.DiGraph(num_steps=num_steps, module_id=list(module_df['module'])[0],
                            module_name=list(module_df['module_name'])[0])
    # go through all path/step combinations
    for module_path, frame in module_df.groupby('path'):
        split_path = [int(i) for i in module_path.split(',')]
        step = split_path[0]
        module_net.add_node(module_path, kos=set(frame['ko']))
        # add incoming edge
        if step != 0:
            module_net.add_edge('end_step_%s' % (step - 1), module_path)
        # add outgoing edge
        module_net.add_edge(module_path, 'end_step_%s' % step)
    return module_net


def get_module_step_coverage(kos, module_net):
    # prune network based on what kos were observed
    pruned_module_net = module_net.copy()
    module_kos_present = set()
    for node, data in module_net.nodes.items():
        if 'kos' in data:
            ko_overlap = data['kos'] & kos
            if len(ko_overlap) == 0:
                pruned_module_net.remove_node(node)
            else:
                module_kos_present = module_kos_present | ko_overlap
    # count number of missing steps
    missing_steps = 0
    for node, data in pruned_module_net.nodes.items():
        if ('end_step' in node) and (pruned_module_net.in_degree(node) == 0):
            missing_steps += 1
    # get statistics
    num_steps = pruned_module_net.graph['num_steps'] + 1
    num_steps_present = num_steps - missing_steps
    coverage = num_steps_present / num_steps
    return num_steps, num_steps_present, coverage, sorted(module_kos_present)


def make_module_coverage_df(annotation_df, module_nets):
    kos_to_genes = defaultdict(list)
    for gene_id, ko_list in annotation_df['ko_id'].iteritems():
        if type(ko_list) is str:
            for ko in ko_list.split(','):
                kos_to_genes[ko].append(gene_id)
    coverage_dict = {}
    for i, (module, net) in enumerate(module_nets.items()):
        module_steps, module_steps_present, module_coverage, module_kos = \
            get_module_step_coverage(set(kos_to_genes.keys()), net)
        module_genes = sorted([gene for ko in module_kos for gene in kos_to_genes[ko]])
        coverage_dict[module] = [net.graph['module_name'], module_steps, module_steps_present, module_coverage,
                                 len(module_kos), ','.join(module_kos), ','.join(module_genes)]
    coverage_df = pd.DataFrame.from_dict(coverage_dict, orient='index',
                                         columns=['module_name', 'steps', 'steps_present', 'step_coverage', 'ko_count',
                                                  'kos_present', 'genes_present'])
    return coverage_df


def make_module_coverage_frame(annotations, module_nets, groupby_column='fasta'):
    # go through each scaffold to check for modules
    module_coverage_dict = dict()
    for group, frame in annotations.groupby(groupby_column, sort=False):
        module_coverage_dict[group] = make_module_coverage_df(frame, module_nets)
    module_coverage = pd.concat(module_coverage_dict)
    module_coverage.index = module_coverage.index.set_names(['genome', 'module'])
    return module_coverage.reset_index()


def make_module_coverage_heatmap(module_coverage, mag_order=None):
    num_mags_in_frame = len(set(module_coverage['genome']))
    c = alt.Chart(module_coverage, title='Module').encode(
        x=alt.X('module_name', title=None, sort=None, axis=alt.Axis(labelLimit=0, labelAngle=90)),
        y=alt.Y('genome', title=None, sort=mag_order, axis=alt.Axis(labelLimit=0)),
        tooltip=[alt.Tooltip('genome', title='Genome'),
                 alt.Tooltip('module_name', title='Module Name'),
                 alt.Tooltip('steps', title='Module steps'),
                 alt.Tooltip('steps_present', title='Steps present')
                 ]
    ).mark_rect().encode(color=alt.Color('step_coverage', legend=alt.Legend(title='% Complete'),
                                         scale=alt.Scale(domain=(0, 1)))).properties(
        width=HEATMAP_CELL_WIDTH * len(HEATMAP_MODULES),
        height=HEATMAP_CELL_HEIGHT * num_mags_in_frame)
    return c


def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def first_open_paren_is_all(str_):
    """Go through string and return true"""
    curr_level = 1
    for char in str_[1:-1]:
        if char == ')':
            curr_level -= 1
        elif char == '(':
            curr_level += 1
        if curr_level == 0:
            return False
    return True


def split_into_steps(definition, split_char=' '):
    """Very fancy split on string of chars"""
    curr_level = 0
    step_starts = [-1]
    for i, char in enumerate(definition):
        if char == '(':
            curr_level += 1
        if char == ')':
            curr_level -= 1
        if (curr_level == 0) and (char in split_char):
            step_starts.append(i)
    step_starts.append(len(definition))
    steps = list()
    for a, b in pairwise(step_starts):
        step = definition[a+1:b]
        if step.startswith('(') and step.endswith(')'):
            if first_open_paren_is_all(step):
                step = step[1:-1]
        steps.append(step)
    return steps


def is_ko(ko):
    return re.match(KO_REGEX, ko) is not None


def make_module_network(definition, network: nx.DiGraph = None, parent_nodes=('start',)):
    # TODO: Figure out how to add 'end' step to last step at end
    if network is None:
        network = nx.DiGraph()
    last_steps = []
    for step in split_into_steps(definition, ','):
        prev_steps = parent_nodes
        for substep in split_into_steps(step, '+'):
            if is_ko(substep):
                for prev_step in prev_steps:
                    network.add_edge(prev_step, substep)
                prev_steps = [substep]
            else:
                network, prev_steps = make_module_network(substep, network, prev_steps)
        last_steps += prev_steps
    return network, last_steps


def get_module_coverage(module_net: nx.DiGraph, genes_present: set):
    max_coverage = -1
    max_coverage_genes = list()
    max_coverage_missing_genes = list()
    max_path_len = 0
    for net_path in nx.all_simple_paths(module_net, source='start', target='end'):
        net_path = set(net_path[1:-1])
        overlap = net_path & genes_present
        coverage = len(overlap) / len(net_path)
        if coverage > max_coverage:
            max_coverage = coverage
            max_coverage_genes = overlap
            max_coverage_missing_genes = net_path - genes_present
            max_path_len = len(net_path)
    return max_path_len, len(max_coverage_genes), max_coverage, max_coverage_genes, max_coverage_missing_genes


def make_etc_coverage_df(etc_module_df, annotations, groupby_column='fasta'):
    etc_coverage_df_rows = list()
    for _, module_row in etc_module_df.iterrows():
        definition = module_row['definition']
        # remove optional subunits
        definition = re.sub(r'-K\d\d\d\d\d', '', definition)
        module_net, _ = make_module_network(definition)
        # add end node
        no_out = [node for node in module_net.nodes() if module_net.out_degree(node) == 0]
        for node in no_out:
            module_net.add_edge(node, 'end')
        # go through each genome and check pathway coverage
        for group, frame in annotations.groupby(groupby_column):
            # get annotation genes
            grouped_ids = set(get_ids_from_annotation(frame).keys())
            path_len, path_coverage_count, path_coverage_percent, genes, missing_genes = \
                get_module_coverage(module_net, grouped_ids)
            complex_module_name = 'Complex %s: %s' % (module_row['complex'].replace('Complex ', ''),
                                                      module_row['module_name'])
            etc_coverage_df_rows.append([module_row['module_id'], module_row['module_name'],
                                         module_row['complex'].replace('Complex ', ''), group, path_len,
                                         path_coverage_count, path_coverage_percent, ','.join(sorted(genes)),
                                         ','.join(sorted(missing_genes)), complex_module_name])
    return pd.DataFrame(etc_coverage_df_rows, columns=ETC_COVERAGE_COLUMNS)


def make_etc_coverage_heatmap(etc_coverage, mag_order=None, module_order=None):
    num_mags_in_frame = len(set(etc_coverage['genome']))
    charts = list()
    for i, (etc_complex, frame) in enumerate(etc_coverage.groupby('complex')):
        # if this is the first chart then make y-ticks otherwise none
        c = alt.Chart(frame, title=etc_complex).encode(
            x=alt.X('module_name', title=None, axis=alt.Axis(labelLimit=0, labelAngle=90),
                    sort=module_order),
            y=alt.Y('genome', axis=alt.Axis(title=None, labels=False, ticks=False), sort=mag_order),
            tooltip=[alt.Tooltip('genome', title='Genome'),
                     alt.Tooltip('module_name', title='Module Name'),
                     alt.Tooltip('path_length', title='Module Subunits'),
                     alt.Tooltip('path_length_coverage', title='Subunits present'),
                     alt.Tooltip('genes', title='Genes present'),
                     alt.Tooltip('missing_genes', title='Genes missing')
                     ]
        ).mark_rect().encode(color=alt.Color('percent_coverage', legend=alt.Legend(title='% Complete'),
                                             scale=alt.Scale(domain=(0, 1)))).properties(
            width=HEATMAP_CELL_WIDTH * len(set(frame['module_name'])),
            height=HEATMAP_CELL_HEIGHT * num_mags_in_frame)
        charts.append(c)
    concat_title = alt.TitleParams('ETC Complexes', anchor='middle')
    return alt.hconcat(*charts, spacing=5, title=concat_title)



def make_functional_df(annotations, function_heatmap_form, groupby_column='fasta'):
    # clean up function heatmap form
    function_heatmap_form = function_heatmap_form.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    function_heatmap_form = function_heatmap_form.fillna('')
    # build dict of ids per genome
    genome_to_id_dict = dict()
    for genome, frame in annotations.groupby(groupby_column, sort=False):
        id_list = get_ids_from_annotation(frame).keys()
        genome_to_id_dict[genome] = set(id_list)
    # build long from data frame
    rows = list()
    for function, frame in function_heatmap_form.groupby('function_name', sort=False):
        for bin_name, id_set in genome_to_id_dict.items():
            presents_in_bin = list()
            functions_present = set()
            for _, row in frame.iterrows():
                function_id_set = set([i.strip() for i in row.function_ids.strip().split(',')])
                present_in_bin = id_set & function_id_set
                functions_present = functions_present | present_in_bin
                presents_in_bin.append(len(present_in_bin) > 0)
            function_in_bin = np.all(presents_in_bin)
            row = frame.iloc[0]
            rows.append([row.category, row.subcategory, row.function_name, ', '.join(functions_present),
                         '; '.join(get_ordered_uniques(frame.long_function_name)),
                         '; '.join(get_ordered_uniques(frame.gene_symbol)), bin_name, function_in_bin,
                         '%s: %s' % (row.category, row.function_name)])
    return pd.DataFrame(rows, columns=list(function_heatmap_form.columns) + ['genome', 'present',
                                                                             'category_function_name'])


def make_functional_heatmap(functional_df, mag_order=None):
    # build heatmaps
    charts = list()
    for i, (group, frame) in enumerate(functional_df.groupby('category', sort=False)):
        # set variables for chart
        function_order = get_ordered_uniques(list(frame.function_name))
        num_mags_in_frame = len(set(frame['genome']))
        chart_width = HEATMAP_CELL_WIDTH * len(function_order)
        chart_height = HEATMAP_CELL_HEIGHT * num_mags_in_frame
        # if this is the first chart then make y-ticks otherwise none
        if i == 0:
            y = alt.Y('genome', title=None, sort=mag_order,
                      axis=alt.Axis(labelLimit=0, labelExpr="replace(datum.label, /_\d*$/gi, '')"))
        else:
            y = alt.Y('genome', axis=alt.Axis(title=None, labels=False, ticks=False), sort=mag_order)
        # set up colors for chart
        rect_colors = alt.Color('present',
                                legend=alt.Legend(title="Function is Present", symbolType='square',
                                                  values=[True, False]),
                                sort=[True, False],
                                scale=alt.Scale(range=['#2ca25f', '#e5f5f9']))
        # define chart
        # TODO: Figure out how to angle title to take up less space
        c = alt.Chart(frame, title=alt.TitleParams(group)).encode(
            x=alt.X('function_name', title=None, axis=alt.Axis(labelLimit=0, labelAngle=90), sort=function_order),
            tooltip=[alt.Tooltip('genome', title='Genome'),
                     alt.Tooltip('category', title='Category'),
                     alt.Tooltip('subcategory', title='Subcategory'),
                     alt.Tooltip('function_ids', title='Function IDs'),
                     alt.Tooltip('function_name', title='Function'),
                     alt.Tooltip('long_function_name', title='Description'),
                     alt.Tooltip('gene_symbol', title='Gene Symbol')]
        ).mark_rect().encode(y=y, color=rect_colors).properties(
            width=chart_width,
            height=chart_height)
        charts.append(c)
    # merge and return
    function_heatmap = alt.hconcat(*charts, spacing=5)
    return function_heatmap


# TODO: refactor this to handle splitting large numbers of genomes into multiple heatmaps here
def fill_liquor_dfs(annotations, module_nets, etc_module_df, function_heatmap_form, groupby_column='fasta'):
    module_coverage_frame = make_module_coverage_frame(annotations, module_nets, groupby_column)

    # make ETC frame
    etc_coverage_df = make_etc_coverage_df(etc_module_df, annotations, groupby_column)

    # make functional frame
    function_df = make_functional_df(annotations, function_heatmap_form, groupby_column)

    return module_coverage_frame, etc_coverage_df, function_df


def rename_genomes_to_taxa(function_df, labels, mag_order):
    function_df = function_df.copy()
    new_genome_column = [labels[i] for i in function_df['genome']]
    function_df['genome'] = pd.Series(new_genome_column, index=function_df.index)
    mag_order = [labels[i] for i in mag_order]
    return function_df, mag_order


def make_liquor_heatmap(module_coverage_frame, etc_coverage_df, function_df, mag_order=None, labels=None):
    module_coverage_heatmap = make_module_coverage_heatmap(module_coverage_frame, mag_order)
    etc_heatmap = make_etc_coverage_heatmap(etc_coverage_df, mag_order=mag_order)
    if labels is not None:
        function_df, mag_order = rename_genomes_to_taxa(function_df, labels, mag_order)
    function_heatmap = make_functional_heatmap(function_df, mag_order)

    liquor = alt.hconcat(alt.hconcat(module_coverage_heatmap, etc_heatmap), function_heatmap)
    return liquor


def make_liquor_df(module_coverage_frame, etc_coverage_df, function_df):
    liquor_df = pd.concat([module_coverage_frame.pivot(index='genome', columns='module_name', values='step_coverage'),
                           etc_coverage_df.pivot(index='genome', columns='complex_module_name',
                                                 values='percent_coverage'),
                           function_df.pivot(index='genome', columns='category_function_name', values='present')],
                          axis=1, sort=False)
    return liquor_df


def get_phylum_and_most_specific(taxa_str):
    taxa_ranks = [i[3:] for i in taxa_str.split(';')]
    phylum = taxa_ranks[1]
    most_specific_rank = TAXONOMY_LEVELS[sum([len(i) > 0 for i in taxa_ranks])-1]
    most_specific_taxa = taxa_ranks[sum([len(i) > 0 for i in taxa_ranks]) - 1]
    if most_specific_rank == 'd':
        return 'd__%s;p__' % most_specific_taxa
    if most_specific_rank == 'p':
        return 'p__%s;c__' % most_specific_taxa
    else:
        return 'p__%s;%s__%s' % (phylum, most_specific_rank, most_specific_taxa)


def make_strings_no_repeats(genome_taxa_dict):
    labels = dict()
    seen = Counter()
    for genome, taxa_string in genome_taxa_dict.items():
        final_taxa_string = '%s_%s' % (taxa_string, str(seen[taxa_string]))
        seen[taxa_string] += 1
        labels[genome] = final_taxa_string
    return labels


def summarize_genomes(annotations,  output_tsv, groupby_column='fasta',
                      camper_distillate=DEFAULT_CAMPER_DIST):
    start_time = datetime.now()

    # read in data
    annotations = pd.read_csv(annotations, sep='\t', index_col=0)
    if 'bin_taxnomy' in annotations:
        annotations = annotations.sort_values('bin_taxonomy')

    genome_summary_form = pd.read_csv(camper_distillate, sep='\t')
    summarized_genomes = make_genome_summary(annotations, genome_summary_form, groupby_column)
    summarized_genomes.drop(['sheet'], axis=1, inplace=True)
    summarized_genomes.to_csv(output_tsv, sep='\t', index=False)



@click.command()
@click.option('-a', '--annotations', help="annotations file, from dram annotate", 
              required=True)
@click.option('-o', '--output_tsv', help="output tsv loctation, it can overwright!", required=True)
@click.option('--groupby_column', type=str, default='fasta',
              help="past_annotations to append new annotations to.")
@click.option('--camper_distillate', default=DEFAULT_CAMPER_DIST,
              help="CAMPER Distillate file, if not specified the default will be used")
def summarize_genomes_cmd(annotations,  output_tsv, groupby_column='fasta',
                      camper_distillate=DEFAULT_CAMPER_DIST):
    summarize_genomes(annotations=annotations,  
                      output_tsv=output_tsv, 
                      groupby_column=groupby_column,
                      camper_distillate=camper_distillate)

if __name__ == "__main__":
    summarize_genomes_cmd()

"""
import os
os.system("python workflow/scripts/camper_distill.py -a output/t1/annotations.tsv -o results/CAMPER_distillate.tsv --custom_distillate CAMPER/test_CAMPER_distillate.tsv")
os.system("rm -rf output/t1/dist")


"""
