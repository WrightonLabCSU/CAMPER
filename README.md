# BETA RELEASE

The CAMPER software and data is currently in a testing state and needs more validations before it should be used for any publication or production task. This RELEASE is for early adopters to get a taste of what CAMPER can do and how it can be incorporated into a workflow. An official release will be available soon.

# License Pending

This software remains unlicensed and the Wrighton Lab reserves all rights for the time being. Contributions will be welcomed once an appropriate license or licenses can be written. This should serve as another reminder that this project is in beta.

# CAMPER (BETA)

<p align="center">
<img width="660" alt="camper_logo" src="https://user-images.githubusercontent.com/95941779/176459947-c5497a73-2b77-42f0-99ee-6a836f23d6a7.png">
</p>


Table of Contents
=================
  * [Overview](https://github.com/WrightonLabCSU/CAMPER#overview)
  * [Installation and Usage](https://github.com/WrightonLabCSU/CAMPER#installation-and-usage)
	* [Using CAMPER within DRAM](https://github.com/WrightonLabCSU/CAMPER#1-using-camper-within-dram)
	* [Standalone CAMPER](https://github.com/WrightonLabCSU/CAMPER#2-camper-standalone-tool-camper_dramkit)
	* [I just want your profiles to search on my own](https://github.com/WrightonLabCSU/CAMPER#3-i-just-want-to-run-your-blast-and-hmm-searches-on-my-own)
  * [CAMPER Outputs](https://github.com/WrightonLabCSU/CAMPER#camper-outputs)
  * [CAMPER Road Map](https://github.com/WrightonLabCSU/CAMPER#camper-map)


# Overview
**C**urated **A**nnotations for **M**icrobial (**P**oly)phenol **E**nzymes and **R**eactions (**CAMPER**) is a tool that annotates genes likely involved in transforming polyphenols, and provides chemical context for these transformations in a summarized form. 

To facilitate the inference of polyphenol metabolism from genomes, CAMPER includes 12 custom Hidden-Markov Model (HMM) profiles and 32 Basic Local Alignment Search Tool (BLAST) searches for (poly)phenol-active genes. We also provide recommended score cut-offs for searches using two ranks: a more stringent, trusted rank (A) and a more relaxed, exploratory rank (B). Beyond these 44 profiles, nearly 300 HMMs from other databases (KEGG/kofamscan, dbCAN2) are included in the CAMPER annotation library.

CAMPER summarizes the gene annotations into 102 modules representing different polyphenol transformations. These modules are classified by the family and sub-family of polyphenols used as substrates (following [Phenol-Explorer](http://phenol-explorer.eu/compounds/classification) Ontology) and by the oxygen requirements for the genes involved. These modules can be as small as a single gene, up to a maximum of 12 genes in the largest module.
![camper_pathways_tree_for_github-01](https://user-images.githubusercontent.com/95941779/171468538-3f2cc169-2170-4612-880b-22ad11d7c9e9.png)
**Figure 1. CAMPER consists of 102 polyphenol transformation modules, organized by substrates (Family and Subfamily) and oxygen requirements. Modules can be composed of 1-12 reaction steps.**

For more detailed information on the organization and outputs, see the [CAMPER Outputs](https://github.com/WrightonLabCSU/CAMPER#camper-outputs) section below, and for module visuals, see the [CAMPER Map](https://github.com/WrightonLabCSU/CAMPER#camper-map).

## CAMPER DATA
The CAMPER data set consists of 5 files, each serving a key role in enabling reproducible annotation of gene data.
  - CAMPER_blast.faa: A fasta file of CAMPER genes used as a target in a BLAST style search provided by mmseqs search.
  - CAMPER.hmm: A HMM file used as the target in an HMM profile search provide by MMseqs profilesearch
  - CAMPER_blast_scores.tsv: Provides the minimum cut off scores for search results and quality ranks with BLAST style searches.
  - CAMPER_hmm_scores.tsv: Provides the minimum cut off scores for search results and quality ranks with HMM Profile searches.
  - CAMPER_distillate.tsv: A custom distillate, for use with DRAM or with CAMPER_DRAMKit, to summarize the annotation results.

To support interpretation, the file `compound_key.txt` can be used to easily identify which pathways are needed for particular compounds. `CAMPER_RoadMap_v1.pdf` can be used to visualize pathways.

---
## 1. CAMPER as a standalone tool

### Set up and Installation
#### A. Setup with Conda

The simplest way to get started with CAMPER is with Conda, using the enviroment.yaml provided in this repository.

CAMPER comes with the latest version of the CAMPER annotation library preloaded, so if all you want to do is annotate and distill called genes with CAMPER, you only need the following commands.

```
  wget https://raw.githubusercontent.com/WrightonLabCSU/CAMPER/main/CAMPER_DRAMKit/environment.yaml
  conda env create --name CAMPER -f ./environment.yaml
```
You can activate the newly made environment with the command:
```
conda activate CAMPER
```

#### B. Setup with pip
If you are not able to use Conda, you can still install CAMPER_DRAMKit with pip using the command below. Note that first you will need to manually install [scikit-bio](http://scikit-bio.org/), and [MMseqs2](https://github.com/soedinglab/mmseqs2), as these tools can't be installed with the other pip dependencies.
```
pip install camper_dramkit
```

### Using CAMPER
CAMPER runs in two steps: first `camper_annotate` and second `camper_distill`. The annotate step takes as input amino acid fasta files, and creates and output directory with a file call `annotations.tsv`. Then, the distill step takes this annotations file and outputs a tsv file. 

```
camper_annotate -i <input faa file> -o <name of output directory> --threads <number of processors>
camper_distill  -a <path to annotations.tsv> -o <name of output.tsv>
```
To give multiple faa files as input, use the wildcard form:
```
camper_annotate -i 'path/to/*faa' -o <name of output directory>
```
For descriptions of the output files, see the [CAMPER Outputs](https://github.com/WrightonLabCSU/CAMPER#camper-outputs) section below.

---
## 2. Using CAMPER with DRAM

If your goal is to integrate CAMPER into your regular genome annotation pipeline, we recommend running it as part of [DRAM](https://github.com/WrightonLabCSU/DRAM). This will provide curated annotation and summarization of polyphenol transformation genes in addition to the regular DRAM databases. **CAMPER will be available as a callable database in DRAM v1.5+**.

There are two steps to running CAMPER in DRAM: (1) annotation and (2) summarization (distillation).
Supply the `--use_camper` flag during the annotation step, like so:

```
DRAM.py annotate --use_camper -i 'my_bins/*.fa' -o DRAM_wCAMPER
DRAM.py distill -i DRAM_wCAMPER/annotations.tsv -o DRAM_wCAMPER_distilled
```
The difference in outputs between this and default DRAM is that you will find CAMPER-specific columns added to the `annotations.tsv` and you will find a CAMPER tab in your `metabolism_summary.xlsx` output.

For descriptions of the content in output files, see the [CAMPER Outputs](https://github.com/WrightonLabCSU/CAMPER#camper-outputs) section below.

#### Installing and running CAMPER in the DRAM environment
If you want to add CAMPER annotations to existing DRAM annotations, it may make sense to install CAMPER in the same Conda environment. This is easy to do if you have already made a DRAM Conda environment with the [instructions in the DRAM README](https://github.com/WrightonLabCSU/DRAM), then you can add CAMPER with the following commands:

```
wget https://raw.githubusercontent.com/WrightonLabCSU/CAMPER/main/CAMPER_DRAMKit/environment.yaml
conda env update --name DRAM -f ./environment.yaml
```

Note: If you install CAMPER, you will get the latest version of the CAMPER database with it. If you want more control over the database, you can override the default data with the instructions in [Other Tools and flags](https://github.com/WrightonLabCSU/CAMPER#other-tools-and-flags).

### DRAM Combination Workflow
It is possible to add CAMPER annotations to genomes annotated with DRAM v1.3+, with one additional command. First follow the instructions above to update your DRAM environment with CAMPER_DRAMKit. Then with that environment activated, you should be able to run the following commands to make a new raw annotations file with all the DRAM data you expect, and the CAMPER data added in.

If you are not able to update your DRAM environment for whatever reason, you will simply need to switch environments mid-workflow.

```
DRAM.py annotate -i 'my_bins/*.fa' -o dram_output
camper_annotate -i 'my_bins/*.fa' -a <path to DRAM/annotations.tsv> -o <output directory to create>
camper_distill  -a <path to camper output/annotations.tsv> -o <name of output.tsv>
```
For descriptions of the `annotations.tsv` and summary file, see the [CAMPER Outputs](https://github.com/WrightonLabCSU/CAMPER#camper-outputs) section below.

---
### Manually Specifying the Location of CAMPER Files
The behavior of the `camper_annotate` and `camper_distill` commands is controlled by the latest version of the CAMPER dataset. If you want to use an older version of CAMPER, it is suggested you install the older version of the CAMPER_DRAMKit tool, as they will be released together and be mutually compatible. However, if you must, you can also specify the files to use with  `camper_annotate` and `camper_distill` using the appropriate arguments. An example is shown below.
```
camper_annotate -i my_genes.faa -o my_output \
	--camper_fa_db_loc CAMPER_blast.faa \
	--camper_fa_db_cutoffs_loc CAMPER_blast_scores.tsv \
	--camper_hmm_loc CAMPER.hmm  \
   --camper_hmm_cutoffs_loc CAMPER_hmm_scores.tsv
camper_distill  -a my_output/annotations.tsv -o my_output/distillate.tsv \
	--camper_distillate CAMPER_distillate.tsv
```
---
# CAMPER Outputs
Approaches [1](https://github.com/WrightonLabCSU/CAMPER#1-using-camper-within-dram) and [2](https://github.com/WrightonLabCSU/CAMPER#2-camper-standalone-tool-camper_dramkit) output two files: the raw information for given searches (`annotations.tsv`) and the summarized information across searches (the distillate, either the `metabolism_summary.xlsx` if run through DRAM or the `distillate.tsv` from CAMPER_DRAMKit).

**Raw annotations**: This is either a standalone file, or columns added to a file, depending on the search approach. This file tells you the genes in your dataset that pass CAMPER annotation thresholds, what they are annotated as, and the scores. It includes the following columns:
  - `camper_hits`, A longer ID giving the CAMPER ID, gene abbreviation, and gene description.
  - `camper_rank`, A match quality rank based on the value of the bit score (A or B). For BLAST-style searches, an A rank is a bitscore >=200 and B >=120. For HMM-style searches, scores are specific to each profile (see `CAMPER_hmm_scores.tsv`).
  - `camper_bitScore`, The bitscore from the best search result. If more than one search meets at least a B-rank for a given gene, the search with the higher score is reported.
  - `camper_id`: Unique CAMPER ID used in the distillation step, of the form D000XX.
  - `camper_definition`: A short description of the CAMPER match in the database.
  - `camper_search_type`: Tells you if a HMM profile or blast search found this match.

**Distillate**: This is either a single file, or the **CAMPER** tab in the `metabolism_summary.xlsx` file. Each row in this file corresponds to a gene in a CAMPER module. This file gives you gene counts of genes in CAMPER modules. It includes the following columns:
 - `gene_id`, the database IDs assigned to this gene. These can be from CAMPER (D000XX), KEGG (KXXXX), dbCAN (AAX), or EC numbers. Note, some IDs are included more than once in the sheet if they are involved in more than one module!
- `gene_description`, A more informative description of the gene in the step, including gene abbreviation and gene name.
- `module`, The CAMPER module that the given gene belongs to. There are 101 modules in CAMPER.
- `header`, The classification for the polyphenol substrate following [Phenol-Explorer](http://phenol-explorer.eu/compounds/classification) Ontology. In the form: Polyphenol;Family;Sub-Family;Compound.
- `subheader`, This contains information about routes, steps, and subunits. Sometimes, a given transformation can be accomplished in more than one sequence of steps: these are termed 'Routes'. Steps indicate the sequential transformations in the module. Subunits denote if the given gene encodes a subunit of a larger complex that carries out a step. Sometimes steps are labeled as "optional" if they are not required.
- `specifc_reaction`, This gives examples of reactions when possible.
- `oxygen`, This is either "oxic", "anoxic"," or "both" for reactions that require oxygen, don't require oxygen, or can function with or without, respectively. Note: these are largely based on literature reporting and the systems they were characterized in, and should be used as guidelines.
- `EC`, The EC number (if known) for a reaction.
- `Notes`, Any important information to know about the genes, for example: manual curation to do, note on gene clusters, should they be extracellular etc.
**The remaining columns will be counts of each gene in your input files.**

# CAMPER Map
This is also provided as a PDF file.
![CAMPER_RoadMap_v1-01](https://user-images.githubusercontent.com/95941779/171663024-4e53989c-f88b-4689-98b6-5f4df2c03b45.png)

# Happy CAMPER-ing!
Annotations, organization, and conceptualization by [Bridget McGivern](https://github.com/bmcgivern13). Coding and implementation by [Rory Flynn](https://github.com/rmFlynn).
