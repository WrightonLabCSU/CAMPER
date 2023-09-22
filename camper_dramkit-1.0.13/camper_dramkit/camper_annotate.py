"""A cut down of dram annotator that only adds CAMPER Annotations"""
import re
import os
import warnings
import subprocess
from glob import glob
from functools import partial
from datetime import datetime
from os import path, mkdir, stat
from typing import Callable
from skbio.io import read as read_sequence
from shutil import rmtree
from skbio.io import write as write_sequence
import pandas as pd
import numpy as np
import click
# Remove

CAMPER_NAME = "CAMPER"
DEFAULT_CUSTOM_FA_DB_LOC = os.path.join(os.path.dirname(__file__),  "data", "CAMPER_blast.faa")
DEFAULT_CUSTOM_HMM_LOC = os.path.join(os.path.dirname(__file__), "data", "CAMPER.hmm")
DEFAULT_CUSTOM_HMM_CUTOFFS_LOC = os.path.join(os.path.dirname(__file__), "data", "CAMPER_hmm_scores.tsv")
DEFAULT_CUSTOM_FA_DB_CUTOFFS_LOC = os.path.join(os.path.dirname(__file__), "data", "CAMPER_blast_scores.tsv")

BOUTFMT6_COLUMNS = ['qId', 'tId', 'seqIdentity', 'alnLen', 'mismatchCnt',
                    'gapOpenCnt', 'qStart', 'qEnd', 'tStart', 'tEnd', 'eVal',
                    'bitScore']
HMMSCAN_ALL_COLUMNS = ['query_id', 'query_ascession', 'query_length', 
                       'target_id', 'target_ascession', 'target_length',
                       'full_evalue', 'full_score', 'full_bias', 
                       'domain_number', 'domain_count', 'domain_cevalue',
                       'domain_ievalue', 'domain_score', 'domain_bias',
                       'target_start', 'target_end', 'alignment_start',
                       'alignment_end', 'query_start', 'query_end', 'accuracy',
                       'description']
HMMSCAN_COLUMN_TYPES = [str, str, int, str, str, int, float, float, float, int,
                        int, float, float, float, float, int, int, int, int, 
                        int, int, float, str]


def remove_suffix(text, suffix):
    if text.endswith(suffix):
        return text[:-1*len(suffix)]
    return text  # or whatever


def merge_files(files_to_merge, outfile, has_header=False):
    """It's in the name, if has_header assumes all files have the same header"""
    with open(outfile, 'w') as outfile_handle:
        if has_header:
            outfile_handle.write(open(files_to_merge[0]).readline())
        for file in files_to_merge:
            with open(file) as f:
                if has_header:
                    _ = f.readline()
                outfile_handle.write(f.read())


def run_process(command, shell=False, capture_stdout=True, check=True, verbose=False):
    """Standardization of parameters for using subprocess.run, provides verbose mode and option to run via shell"""
    stdout = subprocess.DEVNULL
    stderr = subprocess.DEVNULL
    if verbose:
        stdout = None
        stderr = None
    if capture_stdout:
        return subprocess.run(command, check=check, shell=shell, stdout=subprocess.PIPE,
                              stderr=stderr).stdout.decode(errors='ignore')
    else:
        subprocess.run(command, check=check, shell=shell, stdout=stdout, stderr=stderr)


def make_mmseqs_db(fasta_loc, output_loc, create_index=True, threads=10, verbose=False):
    """Takes a fasta file and makes a mmseqs2 database for use in blast searching and hmm searching with mmseqs2"""
    run_process(['mmseqs', 'createdb', fasta_loc, output_loc], verbose=verbose)
    if create_index:
        tmp_dir = path.join(path.dirname(output_loc), 'tmp')
        run_process(['mmseqs', 'createindex', output_loc, tmp_dir, '--threads', str(threads)], verbose=verbose)


def multigrep(search_terms, search_against, split_char='\n', output='.'):
    # TODO: multiprocess this over the list of search terms
    """Search a list of exact substrings against a database, takes name of mmseqs db index with _h to search against"""
    hits_file = path.join(output, 'hits.txt')
    with open(hits_file, 'w') as f:
        f.write('%s\n' % '\n'.join(search_terms))
    results = run_process(['grep', '-a', '-F', '-f', hits_file, search_against], capture_stdout=True, verbose=False)
    processed_results = [i.strip() for i in results.strip().split(split_char)
                         if len(i) > 0]
    # remove(hits_file)
    return {i.split()[0]: i for i in processed_results if i != ''}



def filter_fasta(fasta_loc, min_len=5000, output_loc=None):
    """Removes sequences shorter than a set minimum from fasta files, outputs an object or to a file"""
    kept_seqs = (seq for seq in read_sequence(fasta_loc, format='fasta') if len(seq) >= min_len)
    if output_loc is None:
        return list(kept_seqs)
    else:
        write_sequence(kept_seqs, format='fasta', into=output_loc)


def run_prodigal(fasta_loc, output_dir, mode='meta', trans_table='11', verbose=False):
    """Runs the prodigal gene caller on a given fasta file, outputs resulting files to given directory"""
    output_gff = path.join(output_dir, 'genes.gff')
    output_fna = path.join(output_dir, 'genes.fna')
    output_faa = path.join(output_dir, 'genes.faa')

    run_process(['prodigal', '-i', fasta_loc, '-p', mode, '-g', trans_table, '-f', 'gff', '-o', output_gff, '-a',
                 output_faa, '-d', output_fna], verbose=verbose)
    return output_gff, output_fna, output_faa


def get_best_hits(query_db, target_db, output_dir='.', query_prefix='query', target_prefix='target',
                  bit_score_threshold=60, threads=10, verbose=False):
    """Uses mmseqs2 to do a blast style search of a query db against a target db, filters to only include best hits
    Returns a file location of a blast out format 6 file with search results
    """
    # make query to target db
    tmp_dir = path.join(output_dir, 'tmp')
    query_target_db = path.join(output_dir, '%s_%s.mmsdb' % (query_prefix, target_prefix))
    run_process(['mmseqs', 'search', query_db, target_db, query_target_db, tmp_dir, '--threads', str(threads)],
                verbose=verbose)
    # filter query to target db to only best hit
    query_target_fa_db_top = path.join(output_dir, '%s_%s.tophit.mmsdb' % (query_prefix, target_prefix))
    run_process(['mmseqs', 'filterdb', query_target_db, query_target_fa_db_top, '--extract-lines', '1'], verbose=verbose)
    # filter query to target db to only hits with min threshold
    query_target_fa_db_top_filt = path.join(output_dir, '%s_%s.tophit.minbitscore%s.mmsdb'
                                         % (query_prefix, target_prefix, bit_score_threshold))
    run_process(['mmseqs', 'filterdb', '--filter-column', '2', '--comparison-operator', 'ge', '--comparison-value',
                 str(bit_score_threshold), '--threads', str(threads), query_target_fa_db_top, query_target_fa_db_top_filt],
                verbose=verbose)
    # convert results to blast outformat 6
    forward_output_loc = path.join(output_dir, '%s_%s_hits.b6' % (query_prefix, target_prefix))
    run_process(['mmseqs', 'convertalis', query_db, target_db, query_target_db, forward_output_loc,
                 '--threads', str(threads)], verbose=verbose)
    return forward_output_loc


def get_reciprocal_best_hits(query_db, target_db, output_dir='.', query_prefix='query', target_prefix='target',
                             bit_score_threshold=60, rbh_bit_score_threshold=350, threads=10, verbose=False):
    """Take results from best hits and use for a reciprocal best hits search"""
    # TODO: Make it take query_target_db as a parameter
    # create subset for second search
    query_target_fa_db_top_filt = path.join(output_dir, '%s_%s.tophit.minbitscore%s.mmsdb'
                                         % (query_prefix, target_prefix, bit_score_threshold))  # I DON'T LIKE THIS
    query_target_fa_db_filt_top_swapped = path.join(output_dir, '%s_%s.minbitscore%s.tophit.swapped.mmsdb'
                                                 % (query_prefix, target_prefix, bit_score_threshold))
    # swap queries and targets in results database
    run_process(['mmseqs', 'swapdb', query_target_fa_db_top_filt, query_target_fa_db_filt_top_swapped, '--threads',
                 str(threads)], verbose=verbose)
    target_fa_db_filt = path.join(output_dir, '%s.filt.mmsdb' % target_prefix)
    # create a subdatabase of the target database with the best hits as well as the index of the target database
    run_process(['mmseqs', 'createsubdb', query_target_fa_db_filt_top_swapped, target_db, target_fa_db_filt], verbose=verbose)
    run_process(['mmseqs', 'createsubdb', query_target_fa_db_filt_top_swapped, '%s_h' % target_db,
                 '%s_h' % target_fa_db_filt], verbose=verbose)

    return get_best_hits(target_fa_db_filt, query_db, output_dir, target_prefix, query_prefix, rbh_bit_score_threshold,
                         threads, verbose)


def process_reciprocal_best_hits(forward_output_loc, reverse_output_loc, target_prefix='target'):
    """Process the forward and reverse best hits results to find reverse best hits
    Returns the query gene, target gene, if it was a reverse best hit, % identity, bit score and e-value
    """
    forward_hits = pd.read_csv(forward_output_loc, sep='\t', header=None, names=BOUTFMT6_COLUMNS)
    forward_hits = forward_hits.set_index('qId')
    reverse_hits = pd.read_csv(reverse_output_loc, sep='\t', header=None, names=BOUTFMT6_COLUMNS)
    reverse_hits = reverse_hits.set_index('qId')
    def check_hit(row:pd.Series):
        rbh = False
        if row.tId in reverse_hits.index:
            rbh = row.name == reverse_hits.loc[row.tId].tId
        return {'%s_hit' % target_prefix:      row.tId,
                '%s_RBH' % target_prefix:      rbh,
                '%s_identity' % target_prefix: row.seqIdentity,
                '%s_bitScore' % target_prefix: row.bitScore,
                '%s_eVal' % target_prefix:     row.eVal,
                'index':                       row.name
                }
    hits = forward_hits.apply(check_hit, axis=1, result_type='expand')
    # NOTE these lines may not be necessary
    hits.set_index('index', drop=True, inplace=True)
    hits.index.name = None
    return hits


def get_kegg_description(kegg_hits, header_dict):
    """Gets the KEGG IDs, and full KEGG hits from list of KEGG IDs for output in annotations"""
    gene_description = list()
    ko_list = list()
    for kegg_hit in kegg_hits.kegg_hit:
        header = header_dict[kegg_hit]
        gene_description.append(header)
        kos = re.findall(r'(K\d\d\d\d\d)', header)
        if len(kos) == 0:
            ko_list.append('')
        else:
            ko_list.append(','.join(kos))
    # TODO: change kegg_id to kegg_genes_id so that people get an error and not the wrong identifier
    new_df = pd.DataFrame([kegg_hits['kegg_hit'].values, ko_list, gene_description],
                          index=['kegg_genes_id', 'ko_id', 'kegg_hit'], columns=kegg_hits.index)
    return pd.concat([new_df.transpose(), kegg_hits.drop('kegg_hit', axis=1)], axis=1, sort=False)


def get_uniref_description(uniref_hits, header_dict):
    """Gets UniRef ID's, taxonomy and full string from list of UniRef IDs for output in annotations"""
    gene_description = list()
    uniref_list = list()
    gene_taxonomy = list()
    for uniref_hit in uniref_hits.uniref_hit:
        header = header_dict[uniref_hit]
        gene_description.append(header)
        uniref_list.append(header[header.find('RepID=') + 6:])
        gene_taxonomy.append(re.search(r'Tax=(.*?) (\S*?)=', header).group(1))
    new_df = pd.DataFrame([uniref_list, gene_description, gene_taxonomy],
                          index=['uniref_id', 'uniref_hit', 'uniref_taxonomy'],
                          columns=uniref_hits.index)
    return pd.concat([new_df.transpose(), uniref_hits.drop('uniref_hit', axis=1)], axis=1, sort=False)


def get_basic_description(hits, header_dict, db_name='viral', db_info_path=None):
    """Filter hits, and get descriptions
    (text before first space)"""
    hit_list = list()
    if db_info_path is not None:
        db_info = pd.read_csv(db_info_path, sep='\t', index_col=0)
        db_info.iloc[0]
        hits.iloc[0]
        hits[f"{db_name}_hits"]
        #hits = 
        hits[f"{db_name}_bitScore"]
        hits[f"{db_name}_bitScore"]
        hits = hits.merge(db_info, how='left',left_on=f"{db_name}_hit", 
                           right_index=True)
        hits.assign()
    description = list()
    for hit in hits['%s_hit' % db_name]:
        header = header_dict[hit]
        hit_list.append(hit)
        description.append(header)
    new_df = pd.DataFrame([hit_list, description],
                          index=['%s_id' % db_name, '%s_hit' % db_name],
                          columns=hits.index)
    return pd.concat([new_df.transpose(), hits.drop('%s_hit' % db_name, axis=1)], axis=1, sort=False)


def get_peptidase_description(peptidase_hits, header_dict):
    peptidase_list = list()
    peptidase_family = list()
    peptidase_descirption = list()
    for peptidase_hit in peptidase_hits.peptidase_hit:
        header = header_dict[peptidase_hit]
        peptidase_list.append(peptidase_hit)
        peptidase_family.append(re.search(r'#\w*.#', header).group()[1:-1])
        peptidase_descirption.append(header)
    new_df = pd.DataFrame([peptidase_list, peptidase_family, peptidase_descirption],
                          index=['peptidase_id', 'peptidase_family', 'peptidase_hit'], columns=peptidase_hits.index)
    return pd.concat([new_df.transpose(), peptidase_hits.drop('peptidase_hit', axis=1)], axis=1, sort=False)


def get_sig_row(row, evalue_lim:float=1e-15):
    """Check if hmm match is significant, based on dbCAN described parameters"""
    tstart, tend, tlen, evalue = row[['target_start', 'target_end', 'target_length', 'full_evalue']].values
    perc_cov = (tend - tstart)/tlen
    if perc_cov >= .35 and evalue <= evalue_lim:
        return True
    else:
        return False


# TODO: refactor following to methods to a shared run hmm step and individual get description steps
def parse_hmmsearch_domtblout(file):
    df_lines = list()
    for line in open(file):
        if not line.startswith('#'):
            line = line.split()
            line = line[:22] + [' '.join(line[22:])]
            df_lines.append(line)
    hmmsearch_frame = pd.DataFrame(df_lines, columns=HMMSCAN_ALL_COLUMNS)
    for i, column in enumerate(hmmsearch_frame.columns):
        hmmsearch_frame[column] = hmmsearch_frame[column].astype(HMMSCAN_COLUMN_TYPES[i])
    return hmmsearch_frame


def find_best_dbcan_hit(genome:str, group:pd.DataFrame):
    group['perc_cov'] = group.apply(
        lambda x: (x['target_end'] - x['target_start']) / x['target_length'],
        axis=1)
    group.sort_values('perc_cov', inplace=True)
    group.columns
    group.sort_values('full_evalue', inplace=True)
    return group.iloc[0]['target_id']


def run_hmmscan(genes_faa:str, db_loc:str, db_name:str, output_loc:str, formater:Callable,
                threads:int=2, verbose:bool=False):
    output = path.join(output_loc, f'{db_name}_results.unprocessed.b6')
    run_process(['hmmsearch', '--domtblout', output, '--cpu', str(threads), db_loc, genes_faa], verbose=verbose)
    #definition Parse hmmsearch output
    if not (path.isfile(output) and stat(output).st_size > 0):
        return pd.DataFrame()
    hits = parse_hmmsearch_domtblout(output)
    return formater(hits)


def get_gene_data(fasta_loc):
    """Take the prodigal gene headers and get the scaffold that it came from
    Based on idba_ud 'scaffold_#' scaffold names with gene name after
    """
    df_dict = dict()
    for seq in read_sequence(fasta_loc, format='fasta'):
        split_label = seq.metadata['id'].split('_')
        scaffold = '_'.join(split_label[:-1])
        gene_position = split_label[-1]
        start_position, end_position, strandedness = seq.metadata['description'].split('#')[1:4]
        df_dict[seq.metadata['id']] = [scaffold, int(gene_position), int(start_position), int(end_position),
                                       int(strandedness)]
    return pd.DataFrame.from_dict(df_dict, orient='index', columns=['scaffold', 'gene_position', 'start_position',
                                                                    'end_position', 'strandedness'])


def get_unannotated(fasta_loc, annotations):
    """Get the genes from the fasta which did not get any annotations"""
    return [seq.metadata['id'] for seq in read_sequence(fasta_loc, format='fasta')
            if seq.metadata['id'] not in annotations]


def get_minimum_bitscore(info_db):
    bit_score_threshold = min(info_db[['A_rank', 'B_rank']].min().values)
    return bit_score_threshold


def camper_blast_search(query_db, target_db, working_dir, info_db_path, start_time,
                          db_name, threads=10, verbose=False):
    """A convenience function to do a blast style reciprocal best hits search"""
    # Get kegg hits
    info_db = pd.read_csv(info_db_path, sep='\t', index_col=0)
    print('%s: Getting forward best hits from %s' % (str(datetime.now() - start_time), db_name))
    bit_score_threshold = get_minimum_bitscore(info_db)
    hits_path = get_best_hits(query_db, target_db, working_dir, 'gene', db_name, 
                         bit_score_threshold, threads, verbose=verbose)
    return camper_blast_search_formater(hits_path, db_name, info_db, start_time)

def bitScore_per_row(row):
    if row['score_type'] == 'domain':
        return row.domain_score
    elif row['score_type'] == 'full':
        return row.full_score
    elif row['score_type'] == '-':
        return None
    else:
        raise ValueError("The score_type must be 'domain', 'full', or 'j")

def rank_per_row(row):
    r_a = row['A_rank']
    r_b = row['B_rank']
    score = row['bitScore']
    if score is None:
        return None
    if float(score) >= float(r_a):
         return 'A'
    if pd.isnull(r_b):
        return None
    if float(score) >= float(r_b):
         return 'B'
    return None

#TODO decide if we need use_hmmer_thresholds:bool=False
def camper_hmmscan_formater(hits:pd.DataFrame,  db_name:str, 
                             hmm_info_path:str, top_hit:bool=True):
    if len(hits) < 1:
        return pd.DataFrame(columns = [
            f"{db_name}_id", f"{db_name}_rank", f"{db_name}_bitScore", 
            f"{db_name}_hits", f"{db_name}_search_type"])
    if hmm_info_path is None:
        hmm_info = None
        hits = hits[hits.apply(get_sig_row, axis=1)]
    else:
        hmm_info = pd.read_csv(hmm_info_path, sep='\t', index_col=0)
        hits = hits.merge(hmm_info, how='left',left_on="target_id", right_index=True)
        hits['bitScore'] = hits.apply(bitScore_per_row, axis=1)
        hits['score_rank'] = hits.apply(rank_per_row, axis=1)
        hits.dropna(subset=['score_rank'], inplace=True)
    if len(hits) == 0:
        # if nothing significant then return nothing, don't get descriptions
        return pd.DataFrame()
    if top_hit:
        # Get the best hits
        # TODO check we want top hit
        hits = hits.sort_values('full_evalue').drop_duplicates(subset=["query_id"])
    hits.set_index('query_id', inplace=True, drop=True)
    hits.rename_axis(None, inplace=True)
    if 'definition' in hits.columns:
        hits = hits[['target_id', 'score_rank', 'bitScore', 'definition']]
        hits.columns = [f"{db_name}_id", f"{db_name}_rank", 
                        f"{db_name}_bitScore", f"{db_name}_hits"]
    else:
        hits = hits[['target_id', 'score_rank', 'bitScore']]
        hits.columns = [f"{db_name}_id", f"{db_name}_rank", 
                        f"{db_name}_bitScore"]
    # Rename
    hits[f"{db_name}_search_type"] = 'hmm'
    return hits



def camper_blast_search_formater(hits_path, db_name, info_db, start_time):
    if stat(hits_path).st_size == 0:
        return pd.DataFrame()
    hits = pd.read_csv(hits_path, sep='\t', header=None, 
                       names=BOUTFMT6_COLUMNS, index_col='qId')
    hits = hits.merge(info_db, how='left',left_on="tId", right_index=True)
    rank_col = f"{db_name}_rank"
    hits[rank_col] = hits.apply(rank_per_row, axis=1)
    hits.dropna(subset=[rank_col], inplace=True)
    print('%s: Getting descriptions of hits from %s' % (str(datetime.now() - start_time), db_name))
    hits = hits[['tId', rank_col, 'bitScore', 'ID_for_distillate', 'definition']]
    hits.rename(columns={
            'tId': f"{db_name}_hits",
            'ID_for_distillate': f"{db_name}_id", 
            'bitScore': f"{db_name}_bitScore",  
            'definition': f"{db_name}_definition"
        },
        inplace=True)
    hits[f"{db_name}_search_type"] = 'blast'
    return hits


def count_motifs(gene_faa, motif='(C..CH)'):
    motif_count_dict = dict()
    for seq in read_sequence(gene_faa, format='fasta'):
        motif_count_dict[seq.metadata['id']] = len(list(seq.find_with_regex(motif)))
    return motif_count_dict


def strip_endings(text, suffixes: list):
    for suffix in suffixes:
        if text.endswith(suffix):
            text = text[:len(text) - len(suffix)]
    return text


def process_camper(output_dir, custom_fa_db_loc=(), custom_hmm_loc=(), 
                       custom_name=(), threads=1, verbose=False):
    mkdir(output_dir)
    if (len(custom_fa_db_loc) != len(custom_name)) or (len(custom_hmm_loc) != len(custom_name)):
        raise ValueError('Lengths of custom location list and custom name list must be the same.')
    custom_locs = dict()
    for name, fa_loc, hmm_loc in zip(custom_name, custom_fa_db_loc, custom_hmm_loc):
        new_fa_db_loc = path.join(output_dir, '%s.custom.mmsdb' % name)
        make_mmseqs_db(fa_loc, new_fa_db_loc, threads=threads, verbose=verbose)
        run_process(['hmmpress', '-f', hmm_loc], verbose=verbose)  # all are pressed just in case
        custom_locs[name] = {'fa': new_fa_db_loc, 'hmm': hmm_loc} 
    return custom_locs


def process_custom_hmm_cutoffs(custom_hmm_cutoffs_loc, custom_name, verbose=False):
    if custom_hmm_cutoffs_loc is None:
        return {}
    if len(custom_hmm_cutoffs_loc) < 0:
        return {}
    if custom_name is None:
        raise ValueError("You can't use the custom_hmm_cutoffs_loc argument without the custom_name and"
                         " custom_hmm_locs aguments specified.")
    if len(custom_hmm_cutoffs_loc) != len(custom_name):
        warnings.warn(f"Custom hmm cutoffs and descriptions were only provided to the first {len(custom_hmm_cutoffs_loc)}."
                      " The rest of the custom hmms will use standard cutoffs and have no descriptions.")
    return {custom_name[i]:j for i, j in enumerate(custom_hmm_cutoffs_loc)}

def process_custom_fa_db_cutoffs(custom_fa_db_cutoffs_loc, custom_name, verbose=False):
    if custom_fa_db_cutoffs_loc is None: 
        return {}
    if len(custom_fa_db_cutoffs_loc) < 0:
        return {}
    if custom_name is None:
        raise ValueError("You can't use the custom_fa_db_cutoffs_loc argument without the custom_name and"
                         " custom_fa_db_locs aguments specified.")
    if len(custom_fa_db_cutoffs_loc) != len(custom_name):
        warnings.warn(f"Custom fasta db cutoffs and descriptions were only provided to the first {len(custom_fa_db_cutoffs_loc)}."
                      " The rest of the custom fast dbs will use standard cutoffs and have no descriptions.")
    return {custom_name[i]:j for i, j in enumerate(custom_fa_db_cutoffs_loc)}


def annotate_orf(gene_faa:str, tmp_dir:str, start_time, 
                 custom_locs=(), 
                 custom_hmm_cutoffs_locs=(), 
                 custom_fa_db_cutoffs_locs=(), 
                 bit_score_threshold=60, rbh_bit_score_threshold=350,
                 threads=10, verbose=False):

    # Return a base dataframe so all genes are present
    yield pd.DataFrame( index=[seq.metadata['id'] for seq in read_sequence(gene_faa, format='fasta')])

    # Get kegg hits
    print('%s: Turning genes from prodigal to mmseqs2 db' % str(datetime.now() - start_time))
    query_db = path.join(tmp_dir, 'gene.mmsdb')
    make_mmseqs_db(gene_faa, query_db, create_index=True, threads=threads, verbose=verbose)

    for name, locs in custom_locs.items():
        print('%s: Getting hits from %s' % (str(datetime.now() - start_time), name))
        fasta = camper_blast_search(query_db=query_db, 
                                      target_db=locs['fa'], 
                                      working_dir=tmp_dir, 
                                      info_db_path=custom_fa_db_cutoffs_locs[name],
                                      start_time=start_time, 
                                      db_name=name, 
                                      threads=threads,
                                      verbose=verbose,)
        hmm = run_hmmscan(genes_faa=gene_faa,
                                           db_loc=locs['hmm'],
                                           db_name=name,
                                           threads=threads,
                                           output_loc=tmp_dir,
                                           formater=partial(
                                               camper_hmmscan_formater,
                                               db_name=name,
                                               hmm_info_path=custom_hmm_cutoffs_locs.get(name),
                                               top_hit=True
                                           ))
        full = pd.concat([fasta, hmm])
        yield (full
               .groupby(full.index)
               .apply(
                   lambda x: (x
                              .sort_values(f'{name}_search_type', ascending=True) # make sure hmm is first
                              .sort_values(f'{name}_bitScore', ascending=False)
                              .iloc[0])
                  )
               )



def check_fasta(input_faa):
    fasta_locs = glob(input_faa)
    if len(fasta_locs) == 0:
        raise ValueError('Given fasta locations returns no paths: %s' % input_faa)
    print('%s fastas found' % len(fasta_locs))
    return fasta_locs


def get_fa_db_name(fasta_loc):
    return path.splitext(path.basename(remove_suffix(fasta_loc, '.gz')))[0]

def merge_in_new_annotations(new_annotations:pd.DataFrame, past_annotations:pd.DataFrame):
    past_cols = set(past_annotations.columns)
    new_cols = set(new_annotations.columns)
    if new_cols.issubset(past_cols):
        raise Warning("You have passed an annotations file that containes"
                      " columns matching the new annotations. You most"
                      " likely are annotating with the same database again."
                      " The falowing columns will be replaced in the new"
                      " annotations file:\n%s" %
                      past_cols.intersection(new_cols))
    past_annotations = past_annotations[list(past_cols - new_cols)]
    while set(new_annotations.index) != set(past_annotations.index):
        if np.all([new_annotations.index.str.startswith('genes_')]): 
            new_annotations.index = new_annotations.index.str[6:]
            continue
        raise ValueError("The name in the genes.faa file dose not match the"
                         " annotations provided. Thus the old annotations can't"
                         " be merged to.")
    new_annotations = pd.merge(past_annotations, new_annotations, how="left", left_index=True, right_index=True)


def annotate_genes(input_faa, output_dir='.', bit_score_threshold=60, rbh_bit_score_threshold=350,
                   past_annotations_path:str=None, 
                   camper_fa_db_loc=(DEFAULT_CUSTOM_FA_DB_LOC), 
                   camper_fa_db_cutoffs_loc=(DEFAULT_CUSTOM_FA_DB_CUTOFFS_LOC,), 
                   camper_hmm_loc=(DEFAULT_CUSTOM_HMM_LOC),
                   camper_hmm_cutoffs_loc=(DEFAULT_CUSTOM_HMM_CUTOFFS_LOC), 
                   use_uniref=False, use_vogdb=False, 
                   rename_genes=True, keep_tmp_dir=True, threads=10, verbose=True):
    fasta_locs = glob(input_faa)
    if len(fasta_locs) == 0:
        raise ValueError('Given fasta locations returns no paths: %s' % input_faa)
    print('%s fastas found' % len(fasta_locs))
    # set up
    start_time = datetime.now()
    print('%s: Annotation started' % str(datetime.now()))

    if len(camper_fa_db_loc + camper_hmm_loc) < 1:
        raise ValueError("For some reason there are no database selected to"
                         " annotate against. This is most likely a result of"
                         " bad arguments.")

    mkdir(output_dir)
    tmp_dir = path.join(output_dir, 'working_dir')
    mkdir(tmp_dir)

    # setup camper databases to be searched
    camper_locs = process_camper(output_dir=path.join(tmp_dir, 'custom_dbs'), 
                                       custom_fa_db_loc=camper_fa_db_loc,
                                       custom_hmm_loc=camper_hmm_loc,
                                       custom_name=[CAMPER_NAME],  
                                       threads=threads,
                                       verbose=verbose)
    camper_fa_db_cutoffs_locs = process_custom_fa_db_cutoffs(camper_fa_db_cutoffs_loc, [CAMPER_NAME])
    camper_hmm_cutoffs_locs = process_custom_hmm_cutoffs(camper_hmm_cutoffs_loc, [CAMPER_NAME])
    print('%s: Retrieved database locations and descriptions' % (str(datetime.now() - start_time)))

    # annotate
    annotation_locs = list()
    faa_locs = list()
    # TODO IO is so slow!!! We need to use more memory and make less files
    for fasta_loc in fasta_locs:
        # set up
        fasta_name = get_fa_db_name(fasta_loc)
        fasta_dir = path.join(tmp_dir, fasta_name)
        mkdir(fasta_dir)


        # annotate
        annotations = pd.concat([
            orf for orf in
            annotate_orf(fasta_loc, fasta_dir,
                         start_time, camper_locs,
                         camper_hmm_cutoffs_locs, camper_fa_db_cutoffs_locs, 
                         bit_score_threshold, rbh_bit_score_threshold, 
                         threads, verbose)
                       ], axis=1)

        # add fasta name to frame and index, write file
        annotations.insert(0, 'fasta', fasta_name)
        if rename_genes:
            annotations.index = annotations.fasta + '_' + annotations.index
        annotation_loc = path.join(fasta_dir, 'annotations.tsv')
        annotations.to_csv(annotation_loc, sep='\t')
        annotation_locs.append(annotation_loc)
        
        # List files in the directory and delete those matching the pattern "*mmsdb*"
        for file_name in os.listdir(fasta_dir):
            if "mmsdb" in file_name:
                file_path = os.path.join(fasta_dir, file_name)
                if os.path.isfile(file_path):
                    os.remove(file_path)
                    print(f"Removed: {file_name}")

    # merge
    all_annotations = pd.concat([pd.read_csv(i, sep='\t', index_col=0) for i in annotation_locs], sort=False)
    all_annotations = all_annotations.sort_values('fasta')
    if past_annotations_path is not None:
        past_annotations = pd.read_csv(past_annotations_path, sep='\t', index_col=0)
        merge_in_new_annotations(new_annotations = all_annotations, 
                                 past_annotations = past_annotations)
    all_annotations.to_csv(path.join(output_dir, 'annotations.tsv'), sep='\t')
    if len(faa_locs) > 0:
        merge_files(faa_locs, path.join(output_dir, 'genes.faa'))


    # clean up
    if not keep_tmp_dir:
        rmtree(tmp_dir)

    print("%s: Completed annotations" % str(datetime.now() - start_time))


@click.command("annotate_genes")
@click.option('-i', '--input_faa', help="fasta file, optionally with wildcards to point to "
                                   "individual MAGs", required=True)
@click.option('-o', '--output_dir', help="output directory", required=True)
@click.option('-a', '--past_annotations_path',
              help="past_annotations to append new annotations to.")
# @click.option('--bit_score_threshold', type=int, default=60,
#          help='minimum bitScore of search to retain hits')
# @click.option('--rbh_bit_score_threshold', type=int, default=350,
#          help='minimum bitScore of reverse best hits to retain hits')
@click.option('--camper_fa_db_loc', default=[DEFAULT_CUSTOM_FA_DB_LOC], multiple=True,
         help="Location of CAMPER fastas to annotate against")
@click.option('--camper_fa_db_cutoffs_loc', default=[DEFAULT_CUSTOM_FA_DB_CUTOFFS_LOC], multiple=True,
              help="Location of file with camper fasta cutoffs and descriptions")
@click.option('--camper_hmm_loc',   default=[DEFAULT_CUSTOM_HMM_LOC], multiple=True,
              help="Location of camper hmms to annotate against")
@click.option('--camper_hmm_cutoffs_loc', default=[DEFAULT_CUSTOM_HMM_CUTOFFS_LOC], multiple=True,
              help="Location of file with CAMPER HMM cutoffs and descriptions")
# @click.option('--use_uniref', default=False,
#               help='Annotate these fastas against UniRef, drastically increases run time and '
#                    'memory requirements')
@click.option('--keep_tmp_dir', default=False)
@click.option('--threads', type=int, default=10, help='number of processors to use')
def annotate_genes_cmd(input_faa, output_dir='.', bit_score_threshold=60, rbh_bit_score_threshold=350,
                   past_annotations_path:str=None, 
                   camper_fa_db_loc=(DEFAULT_CUSTOM_FA_DB_LOC), 
                   camper_fa_db_cutoffs_loc=(DEFAULT_CUSTOM_FA_DB_CUTOFFS_LOC,), 
                   camper_hmm_loc=(DEFAULT_CUSTOM_HMM_LOC),
                   camper_hmm_cutoffs_loc=(DEFAULT_CUSTOM_HMM_CUTOFFS_LOC), 
                   use_uniref=False, use_vogdb=False, 
                   rename_genes=True, keep_tmp_dir=True, threads=10, verbose=True):
    annotate_genes(
        input_faa=input_faa,
        output_dir=output_dir,
        bit_score_threshold=bit_score_threshold,
        rbh_bit_score_threshold=rbh_bit_score_threshold,
        past_annotations_path=past_annotations_path,
        camper_fa_db_loc=camper_fa_db_loc,
        camper_fa_db_cutoffs_loc=camper_fa_db_cutoffs_loc,
        camper_hmm_loc=camper_hmm_loc,
        camper_hmm_cutoffs_loc=camper_hmm_cutoffs_loc,
        use_uniref=use_uniref,
        use_vogdb=use_vogdb,
        rename_genes=rename_genes,
        keep_tmp_dir=keep_tmp_dir,
        threads=threads,
        verbose=verbose)

if __name__ == "__main__":
    annotate_genes_cmd()

