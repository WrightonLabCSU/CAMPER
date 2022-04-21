import pytest

import os
from io import StringIO
from datetime import datetime
from functools import partial
from filecmp import cmp
import time
from shutil import copy

import pandas as pd
from skbio.io import read as read_sequence

from camper_dramkit.camper_annotate import filter_fasta, run_prodigal, get_best_hits, \
    get_reciprocal_best_hits, process_reciprocal_best_hits, get_kegg_description, get_uniref_description, \
    get_basic_description, get_peptidase_description, get_sig_row, get_gene_data, get_unannotated, \
    count_motifs, strip_endings, process_camper, parse_hmmsearch_domtblout,  camper_hmmscan_formater, \
    make_mmseqs_db, camper_blast_search_formater, annotate_genes, DEFAULT_CUSTOM_FA_DB_LOC, \
    DEFAULT_CUSTOM_FA_DB_CUTOFFS_LOC, DEFAULT_CUSTOM_HMM_LOC, DEFAULT_CUSTOM_HMM_CUTOFFS_LOC

ALT_CUSTOM_FA_DB_LOC = os.path.join("..", "CAMPER_blast.faa")
ALT_CUSTOM_HMM_LOC = os.path.join("..", "CAMPER.hmm")
ALT_CUSTOM_HMM_CUTOFFS_LOC = os.path.join("..", "CAMPER_hmm_scores.tsv")
ALT_CUSTOM_FA_DB_CUTOFFS_LOC = os.path.join("..", "CAMPER_blast_scores.tsv")


@pytest.fixture()
def prodigal_dir(fasta_loc, tmpdir):
    prodigal_output = tmpdir.mkdir('prodigal_output')
    gff, fna, faa = run_prodigal(fasta_loc, str(prodigal_output))
    return gff, fna, faa

@pytest.fixture()
def fasta_loc():
    return os.path.join('tests', 'test_data', 'NC_001422.fasta')


def test_camper_hmmscan_formater():
    hits = pd.read_csv(
        'tests/test_data/test_camper_hmmscan_formater_hits.csv', 
        index_col=0)
    received = camper_hmmscan_formater(
        hits=hits, db_name='CAMPER', 
        hmm_info_path="tests/test_data/test_camper_hmmscan_formater_scores.tsv", 
        top_hit=False)
    expected = pd.read_csv(
        'tests/test_data/test_camper_hmmscan_formater_out_notop_hit.csv',
        index_col=0)
    assert received.equals(expected), "Something wrong with camper hmm"
    received = camper_hmmscan_formater(
        hits=hits, db_name='CAMPER', 
        hmm_info_path="tests/test_data/test_camper_hmmscan_formater_scores.tsv", 
        top_hit=True)
    expected = pd.read_csv(
        'tests/test_data/test_camper_hmmscan_formater_out.csv',
        index_col=0)
    assert received.equals(expected), "Something wrong with camper hmm"


def test_camper_blast_search_formater():
    received = camper_blast_search_formater(
        hits_path='tests/test_data/test_camper_blast_search_formater_hits.tsv',
        db_name='CAMPER', 
        info_db = pd.read_csv(
            './tests/test_data/test_camper_blast_search_formater_scores.tsv', 
            sep='\t', 
            index_col=0) ,
        start_time = datetime.now())
    expected = pd.read_csv(
         './tests/test_data/test_camper_blast_search_formater_out.csv',
         index_col=0)
    assert received.equals(expected), "Something wrong with camper blast search"


@pytest.fixture()
def prodigal_gff(prodigal_dir):
    return prodigal_dir[0]


def test_annotate_genes(tmpdir):
    tmp_out = tmpdir.mkdir('annotate_genes')
    input_faa = os.path.join('tests', 'test_data', 'camper_test_genes.faa')
    input_faa_empty = os.path.join('tests', 'test_data', 'camper_test_genes_empty.faa')
    output_dir=os.path.join(tmp_out, 'annotations')
    output_dir_empty=os.path.join(tmp_out, 'annotations_empty')
    camper_fa_db_loc = DEFAULT_CUSTOM_FA_DB_LOC \
        if os.path.exists(DEFAULT_CUSTOM_FA_DB_LOC) else ALT_CUSTOM_FA_DB_LOC
    camper_fa_db_cutoffs_loc = DEFAULT_CUSTOM_FA_DB_CUTOFFS_LOC \
        if os.path.exists(DEFAULT_CUSTOM_FA_DB_CUTOFFS_LOC) else ALT_CUSTOM_FA_DB_CUTOFFS_LOC
    camper_hmm_loc = DEFAULT_CUSTOM_HMM_LOC \
        if os.path.exists(DEFAULT_CUSTOM_HMM_LOC) else ALT_CUSTOM_HMM_LOC
    camper_hmm_cutoffs_loc = DEFAULT_CUSTOM_HMM_CUTOFFS_LOC \
        if os.path.exists(DEFAULT_CUSTOM_HMM_CUTOFFS_LOC) else ALT_CUSTOM_HMM_CUTOFFS_LOC
    annotate_genes(input_faa=input_faa, output_dir=output_dir,
                   camper_fa_db_loc=(camper_fa_db_loc,), 
                   camper_fa_db_cutoffs_loc=(camper_fa_db_cutoffs_loc,), 
                   camper_hmm_loc=(camper_hmm_loc,),
                   camper_hmm_cutoffs_loc=(camper_hmm_cutoffs_loc,),
                   keep_tmp_dir=False
                   )
    assert (
        pd.read_csv(os.path.join(output_dir, "annotations.tsv"), 
                    sep='\t', index_col=0)
        .equals(
            pd.read_csv(os.path.join('tests', 'test_data', 'camper_test_annotations.tsv'), 
                        sep='\t', index_col=0)))
    # TODO test if fa file is empty
    # annotate_genes(input_faa=input_faa_empty, output_dir=output_dir_empty,
    #                camper_fa_db_loc=(camper_fa_db_loc,), 
    #                camper_fa_db_cutoffs_loc=(camper_fa_db_cutoffs_loc,), 
    #                camper_hmm_loc=(camper_hmm_loc,),
    #                camper_hmm_cutoffs_loc=(camper_hmm_cutoffs_loc,),
    #                keep_tmp_dir=False
    #                )
    # assert (
    #     pd.read_csv(os.path.join(output_dir_empty, "annotations.tsv"), 
    #                 sep='\t', index_col=0)
    #     .equals(
    #         pd.read_csv(os.path.join('tests', 'test_data', 'camper_test_annotations_empty.tsv'), 
    #                     sep='\t', index_col=0)))

@pytest.fixture()
def prodigal_fna(prodigal_dir):
    return prodigal_dir[1]


@pytest.fixture()
def prodigal_faa(prodigal_dir):
    return prodigal_dir[2]


def test_run_prodigal(prodigal_gff, prodigal_fna, prodigal_faa):
    assert os.path.isfile(prodigal_gff)
    assert os.path.isfile(prodigal_fna)
    assert os.path.isfile(prodigal_faa)


@pytest.fixture()
def mmseqs_db_dir(tmpdir):
    output_loc = tmpdir.mkdir('make_mmseqs_db_test')
    return output_loc


@pytest.fixture()
def mmseqs_db(prodigal_faa, mmseqs_db_dir):
    output_file = str(mmseqs_db_dir.join('mmseqs_db.mmsdb'))
    make_mmseqs_db(prodigal_faa, output_file, True, 1)
    return output_file


@pytest.fixture()
def phix_proteins():
    return os.path.join('tests', 'test_data', 'NC_001422.faa')


@pytest.fixture()
def target_mmseqs_db(mmseqs_db_dir, phix_proteins):
    output_file = str(mmseqs_db_dir.join('target.mmsdb'))
    make_mmseqs_db(phix_proteins, output_file, True, 1)
    return output_file


@pytest.fixture()
def best_hits_loc(mmseqs_db, target_mmseqs_db, mmseqs_db_dir):
    best_hits_loc = get_best_hits(mmseqs_db, target_mmseqs_db, mmseqs_db_dir, threads=1, verbose=False)
    return best_hits_loc


def test_get_best_hits(best_hits_loc):
    assert os.path.isfile(best_hits_loc)


@pytest.fixture()
def reverse_best_hits_loc(best_hits_loc, mmseqs_db, target_mmseqs_db, mmseqs_db_dir):
    reverse_best_hits_loc = get_reciprocal_best_hits(mmseqs_db, target_mmseqs_db, mmseqs_db_dir, threads=1,
                                                     verbose=False)
    return reverse_best_hits_loc


def test_get_reciprocal_best_hits(reverse_best_hits_loc):
    assert os.path.isfile(reverse_best_hits_loc)


@pytest.fixture()
def processed_hits():
    forward = os.path.join('tests', 'test_data', 'query_target_hits.b6')
    reverse = os.path.join('tests', 'test_data', 'target_query_hits.b6')
    processed_hits = process_reciprocal_best_hits(forward, reverse)
    return processed_hits


def test_process_reciprocal_best_hits(processed_hits):
    assert processed_hits.shape == (7, 5)
    assert set(processed_hits.loc[processed_hits.target_RBH].index) == {'NC_001422.1_5', 'NC_001422.1_4',
                                                                        'NC_001422.1_7', 'NC_001422.1_6'}


def test_get_kegg_description():
    header_dict = {'aad:TC41_2367': 'aad:TC41_2367  ABC-type molybdate transport system periplasmic component-like '
                                    'protein',
                   'aar:Acear_0854': 'aar:Acear_0854  hypothetical protein; K05810 conserved hypothetical protein',
                   'aar:Acear_1520': 'aar:Acear_1520  hypothetical protein'}
    kegg_hits_data = [['aad:TC41_2367', 10e-5],
                      ['aar:Acear_0854', 10e-6],
                      ['aar:Acear_1520', 10e-10]]
    kegg_hits = pd.DataFrame(kegg_hits_data, index=['gene1', 'gene2', 'gene3'], columns=['kegg_hit', 'eVal'])
    kegg_hits_add_description = get_kegg_description(kegg_hits, header_dict)
    print(kegg_hits_add_description.head())
    assert kegg_hits_add_description.shape == (3, 4)
    assert kegg_hits_add_description.loc['gene2', 'ko_id'] == 'K05810'
    assert kegg_hits_add_description.loc['gene2', 'kegg_genes_id'] == 'aar:Acear_0854'
    assert kegg_hits_add_description.loc['gene1', 'ko_id'] == ''
    assert kegg_hits_add_description.loc['gene1', 'kegg_genes_id'] == 'aad:TC41_2367'


def test_get_uniref_description():
    header_dict = {'UniRef90_A0A139CGD2': 'UniRef90_A0A139CGD2 Phosphate transport system permease protein PstA n=2 '
                                          'Tax=Candidatus Frackibacter TaxID=2017975 RepID=A0A139CGD2_9FIRM',
                   'UniRef90_2642661139': 'UniRef90_2642661139 Ga0073286_10147 conserved hypothetical protein n=1 '
                                          'Tax=Fuchsiella alkaliacetigena WG11 RepID=Ga0073286_10147'}
    uniref_hits_data = [['UniRef90_A0A139CGD2', 10e-5],
                        ['UniRef90_2642661139', 10e-20]]
    uniref_hits = pd.DataFrame(uniref_hits_data, index=['gene1', 'gene2'], columns=['uniref_hit', 'eVal'])
    uniref_hits_add_description = get_uniref_description(uniref_hits, header_dict)
    assert uniref_hits_add_description.shape == (2, 4)
    assert uniref_hits_add_description.loc['gene1', 'uniref_id'] == 'A0A139CGD2_9FIRM'
    assert uniref_hits_add_description.loc['gene1', 'uniref_taxonomy'] == 'Candidatus Frackibacter'


def test_get_basic_description():
    header_dict = {'YP_009015653.1': 'YP_009015653.1 gp350 [Bacillus virus G]',
                   'NP_077550.1': 'NP_077550.1 EsV-1-65 [Ectocarpus siliculosus virus 1]'}
    viral_hits_data = [['NP_077550.1'],
                       ['YP_009015653.1']]
    viral_hits = pd.DataFrame(viral_hits_data, index=['gene1', 'gene2'], columns=['viral_hit'])
    viral_hits_add_description = get_basic_description(viral_hits, header_dict)
    assert viral_hits_add_description.shape == (2, 2)
    assert viral_hits_add_description.loc['gene1', 'viral_id'] == 'NP_077550.1'


def test_get_peptidase_description():
    header_dict = {'MER0025711': 'MER0025711 - family S12 unassigned peptidases (Cytophaga hutchinsonii) [S12.UPW]#S12#'
                                 '{peptidase unit: 588-980}~source ZP_00119535~',
                   'MER0068848': 'MER0068848 - family C56 non-peptidase homologues (Silicibacter sp. TM1040) [C56.UNW]#'
                                 'C56#{peptidase unit: 31-167}~source YP_611907~'}
    peptidase_hits_data = [['MER0025711'], ['MER0068848']]
    peptidase_hits = pd.DataFrame(peptidase_hits_data, index=['gene1', 'gene2'], columns=['peptidase_hit'])
    peptidase_hits_add_description = get_peptidase_description(peptidase_hits, header_dict)
    assert peptidase_hits_add_description.shape == (2, 3)
    assert peptidase_hits_add_description.loc['gene1', 'peptidase_id'] == 'MER0025711'
    assert peptidase_hits_add_description.loc['gene1', 'peptidase_family'] == 'S12'
    assert peptidase_hits_add_description.loc['gene2', 'peptidase_hit'] == 'MER0068848 - family C56 non-peptidase ' \
                                                                           'homologues (Silicibacter sp. TM1040) [C56' \
                                                                           '.UNW]#C56#{peptidase unit: 31-167}~source' \
                                                                           ' YP_611907~'


def test_get_sig_row():
    names = ['target_start', 'target_end', 'target_length', 'full_evalue']
    assert not get_sig_row(pd.Series([1, 85, 100, 1], index=names))
    assert get_sig_row(pd.Series([1, 86, 100, 1e-16], index=names))
    assert not get_sig_row(pd.Series([1, 29, 100, 1e-20], index=names))

@pytest.fixture()
def phix_prodigal_genes():
    phix_seq = ">NC_001422.1_1 # 51 # 221 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-" \
               "10bp;gc_cont=0.404\n" \
               "MSRKIILIKQELLLLVYELNRSGLLAENEKIRPILAQLEKLLLCDLSPSTNDSVKN*\n" \
               ">NC_001422.1_2 # 390 # 848 # 1 # ID=1_2;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=" \
               "5-10bp;gc_cont=0.468\n" \
               "MSQVTEQSVRFQTALASIKLIQASAVLDLTEDDFDFLTSNKVWIATDRSRARRCVEACVY\n" \
               "GTLDFVGYPRFPAPVEFIAAVIAYYVHPVNIQTACLIMEGAEFTENIINGVERPVKAAEL\n" \
               "FAFTLRVRAGNTDVLTDAEENVRQKLRAEGVM*\n" \
               ">NC_001422.1_3 # 848 # 964 # 1 # ID=1_3;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;" \
               "gc_cont=0.496\n" \
               "MSKGKKRSGARPGRPQPLRGTKGKRKGARLWYVGGQQF*\n" \
               ">NC_001422.1_4 # 1001 # 2284 # 1 # ID=1_4;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5" \
               "-10bp;gc_cont=0.449\n" \
               "MSNIQTGAERMPHDLSHLGFLAGQIGRLITISTTPVIAGDSFEMDAVGALRLSPLRRGLA\n" \
               "IDSTVDIFTFYVPHRHVYGEQWIKFMKDGVNATPLPTVNTTGYIDHAAFLGTINPDTNKI\n" \
               "PKHLFQGYLNIYNNYFKAPWMPDRTEANPNELNQDDARYGFRCCHLKNIWTAPLPPETEL\n" \
               "SRQMTTSTTSIDIMGLQAAYANLHTDQERDYFMQRYHDVISSFGGKTSYDADNRPLLVMR\n" \
               "SNLWASGYDVDGTDQTSLGQFSGRVQQTYKHSVPRFFVPEHGTMFTLALVRFPPTATKEI\n" \
               "QYLNAKGALTYTDIAGDPVLYGNLPPREISMKDVFRSGDSSKKFKIAEGQWYRYAPSYVS\n" \
               "PAYHLLEGFPFIQEPPSGDLQERVLIRHHDYDQCFQSVQLLQWNSQVKFNVTVYRNLPTT\n" \
               "RDSIMTS*\n" \
               ">NC_001422.1_5 # 2395 # 2922 # 1 # ID=1_5;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10" \
               "bp;gc_cont=0.420\n" \
               "MFQTFISRHNSNFFSDKLVLTSVTPASSAPVLQTPKATSSTLYFDSLTVNAGNGGFLHCI\n" \
               "QMDTSVNAANQVVSVGADIAFDADPKFFACLVRFESSSVPTTLPTAYDVYPLNGRHDGGY\n" \
               "YTVKDCVTIDVLPRTPGNNVYVGFMVWSNFTATKCRGLVSLNQVIKEIICLQPLK*\n" \
               ">NC_001422.1_6 # 2931 # 3917 # 1 # ID=1_6;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=" \
               "5-10bp;gc_cont=0.451\n" \
               "MFGAIAGGIASALAGGAMSKLFGGGQKAASGGIQGDVLATDNNTVGMGDAGIKSAIQGSN\n" \
               "VPNPDEAAPSFVSGAMAKAGKGLLEGTLQAGTSAVSDKLLDLVGLGGKSAADKGKDTRDY\n" \
               "LAAAFPELNAWERAGADASSAGMVDAGFENQKELTKMQLDNQKEIAEMQNETQKEIAGIQ\n" \
               "SATSRQNTKDQVYAQNEMLAYQQKESTARVASIMENTNLSKQQQVSEIMRQMLTQAQTAG\n" \
               "QYFTNDQIKEMTRKVSAEVDLVHQQTQNQRYGSSHIGATAKDISNVVTDAASGVVDIFHG\n" \
               "IDKAVADTWNNFWKDGKADGIGSNLSRK*\n" \
               ">NC_001422.1_7 # 3981 # 5384 # 1 # ID=1_7;partial=01;start_type=ATG;rbs_motif=GGAGG;rbs_spacer=5-10" \
               "bp;gc_cont=0.455\n" \
               "MVRSYYPSECHADYFDFERIEALKPAIEACGISTLSQSPMLGFHKQMDNRIKLLEEILSF\n" \
               "RMQGVEFDNGDMYVDGHKAASDVRDEFVSVTEKLMDELAQCYNVLPQLDINNTIDHRPEG\n" \
               "DEKWFLENEKTVTQFCRKLAAERPLKDIRDEYNYPKKKGIKDECSRLLEASTMKSRRGFA\n" \
               "IQRLMNAMRQAHADGWFIVFDTLTLADDRLEAFYDNPNALRDYFRDIGRMVLAAEGRKAN\n" \
               "DSHADCYQYFCVPEYGTANGRLHFHAVHFMRTLPTGSVDPNFGRRVRNRRQLNSLQNTWP\n" \
               "YGYSMPIAVRYTQDAFSRSGWLWPVDAKGEPLKATSYMAVGFYVAKYVNKKSDMDLAAKG\n" \
               "LGAKEWNNSLKTKLSLLPKKLFRIRMSRNFGMKMLTMTNLSTECLIQLTKLGYDATPFNQ\n" \
               "ILKQNAKREMRLRLGKVTVADVLAAQPVTTNLLKFMRASIKMIGVSNL\n"
    return StringIO(phix_seq)


# Do this better
def test_get_gene_data(phix_prodigal_genes):
    scaffold_gene_df = get_gene_data(phix_prodigal_genes)
    assert scaffold_gene_df.shape == (7, 5)


def test_get_unannotated(phix_proteins):
    annotated_genes = ['NP_040704.1', 'NP_040703.1', 'NP_040713.1', 'NP_040712.1', 'NP_040711.1', 'NP_040710.1',
                       'NP_040709.1', 'NP_040707.1']
    unannotated_genes = ['NP_040705.1', 'NP_040706.1', 'NP_040708.1']
    test_unannotated_genes = get_unannotated(phix_proteins, annotated_genes)
    assert set(unannotated_genes) == set(test_unannotated_genes)


@pytest.fixture()
def phix_annotations():
    return pd.DataFrame([['A', 'K1', 'U1', None, 'a_bug1'],
                         ['B', 'K2', 'U2', 'P2', 'a_bug2'],
                         ['C', None, 'U3', 'P3', 'a_bug3'],
                         ['D', 'K4', 'U4', 'P4', 'a_bug4'],
                         ['A', 'K5', 'U5', 'P5', 'a_bug5'],
                         ['E', 'K6', 'U6', 'P6', 'a_bug6'],
                         ['C', 'K7', 'U7', 'P7', 'a_bug7']],
                        index=['NC_001422.1_1', 'NC_001422.1_2', 'NC_001422.1_3', 'NC_001422.1_4', 'NC_001422.1_5',
                               'NC_001422.1_6', 'NC_001422.1_7'],
                        columns=['rank', 'kegg_hit', 'uniref_hit', 'pfam_hits', 'bin_taxonomy'])


@pytest.fixture()
def phix_annotations_no_kegg():
    return pd.DataFrame([['B', 'U1', None, 'a_bug1'],
                         ['B', 'U2', 'P2', 'a_bug2'],
                         ['C', 'U3', 'P3', 'a_bug3'],
                         ['D', 'U4', 'P4', 'a_bug4'],
                         ['B', 'U5', 'P5', 'a_bug5'],
                         ['E', 'U6', 'P6', 'a_bug6'],
                         ['C', 'U7', 'P7', 'a_bug7']],
                        index=['NC_001422.1_1', 'NC_001422.1_2', 'NC_001422.1_3', 'NC_001422.1_4', 'NC_001422.1_5',
                               'NC_001422.1_6', 'NC_001422.1_7'],
                        columns=['rank', 'uniref_hit', 'pfam_hits', 'bin_taxonomy'])


def test_count_motifs(phix_proteins):
    motif_counts = count_motifs(phix_proteins, motif='(A.A)')
    assert len(motif_counts) == 11
    assert motif_counts['NP_040713.1'] == 7
    assert motif_counts['NP_040709.1'] == 0


def test_strip_endings():
    assert strip_endings('abc.efg', ['.efg', '.jkl']) == 'abc'
    assert strip_endings('abc.jkl', ['.efg', '.jkl']) == 'abc'
    assert strip_endings('123456', ['.efg', '.jkl']) == '123456'


def test_parse_hmmsearch_domtblout():
    parsed_hit = parse_hmmsearch_domtblout(os.path.join('tests', 'test_data', 'hmmsearch_hit.txt'))
    assert parsed_hit.shape == (1, 23)
    assert parsed_hit.loc[0, 'query_id'] == 'NP_040710.1'
    assert parsed_hit.loc[0, 'query_length'] == 38
    assert parsed_hit.loc[0, 'target_id'] == 'Microvir_J'
    assert parsed_hit.loc[0, 'target_ascession'] == 'PF04726.13'
    assert parsed_hit.loc[0, 'full_evalue'] == 6.900000e-31


@pytest.fixture()
def fake_phix_annotations():
    return pd.DataFrame([[pd.np.NaN],
                         ['GH13'],
                         [pd.np.NaN],
                         [pd.np.NaN],
                         [pd.np.NaN],
                         [pd.np.NaN],
                         [pd.np.NaN]],
                        index=['NC_001422.1_1', 'NC_001422.1_2', 'NC_001422.1_3', 'NC_001422.1_4',
                               'NC_001422.1_5', 'NC_001422.1_6', 'NC_001422.1_7'],
                        columns=['cazy_id'])


