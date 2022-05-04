
# BETA RELEASE

The CAMPER software and data is currently in a testing state and needs more validations before it should be used for any publication or production task. This RELEASE is for early adopters to get a taste of what CAMPER can do and how it can be incorporated into a workflow. A official release will be available soon.

# License Pending

This software remains unlicensed and the Wrighton Lab reserves all rights for the time being. Contributions will be welcomed once an appropriate license or licenses can be written. This should serve as another reminder that this project is in beta.

# CAMPER (BETA)
 
  A revolutionary new data set, to enable the understanding of polyphenol metabolism in its true complexity.

<img width="741" alt="Screen Shot 2022-03-07 at 1 47 10 PM" src="https://user-images.githubusercontent.com/95941779/157345312-27679138-c32c-4e76-8923-a2c776bccbe9.png">


## CAMPER DATA
The CAMPER data set consists of 4 files, each serving a key role in enabling reproducible, intelligent annotation of gene data.
  - CAMPER_blast.fa: A fa file of CAMPER genes used as a target in a BLAST style search provided by mmseqs search.
  - CAMPER.hmm: A HMM file used as the target in an HMM profile search provide by MMseqs profilesearch
  - CAMPER_blast_scores.tsv: Determines minimum cut off scores for search results and quality ranks with blast style searches.
  - CAMPER_hmm_scores.tsv: Determines minimum cut off scores for search results and quality ranks with HMM Profile searches.
  - CAMPER_distillate.tsv: A custom distillate, for use with DRAM or with CAMPER_DRAMKit
## Using CAMPER as part of DRAM

To facilitate the use of the CAMPER data set, we are integrating it into [DRAM](https://github.com/WrightonLabCSU/DRAM).

With the release of DRAM1.4.0 you will be able to use the CAMPER data set without any additional data sets simply by supplying the `--use_camper` flag using the annotate command, like so:

```
DRAM.py annotate --use_camper -i 'my_bins/*.fa' -o annotation
```

When you run `DRAM.py distill` it should detect that CAMPER data was included, and you will find a CAMPER tab in your metabolic summary.

If you want to use CAMPER as part of your annotation pipeline along with an older version of DRAM see the (With DRAM) section below.

## CAMPER Standalone Tool: CAMPER_DRAMKit

In addition to the CAMPER module in DRAM, we have also created a standalone tool, named CAMPER_DRAMKit, that can be used by itself or in conjunction with older versions of DRAM. In Fact, CAMPER_DRAMKit is really a much smaller version of DRAM that follows much the same workflow as DRAM and has similar capabilities. The setup will be a bit different for each use case, so you may want to review the options before you begin. If you are wondering what this tool can do for you, you can jump to the [Usage](https://github.com/WrightonLabCSU/CAMPER/edit/main/README.md#usage) section and get a better idea for how the tool works.

#### Setup With Conda

If you are only interested in CAMPER, or don't want to wait for the DRAM1.4 release, you can use the standalone tool CAMPER_DRAMKit. The simplest way to get started is with Conda, using the enviroment.yaml provided in this repository.

CAMPER_DRAMKit comes with the latest version of CAMPER preloaded, so if all you want to do is annotate and distill called genes with CAMPER, you only need the following commands.

```
wget https://github.com/WrightonLabCSU/CAMPER/main/CAMPER_DRAMKit/environment.yaml
conda create --name CAMPER -f ./environment.yaml
```

However, if you want to have all of CAMPER at your fingertips you can use the following command to both download this repository, including all the CAMPER data and install the CAMPER_DRAMKit tool.

```
git clone https://github.com/WrightonLabCSU/CAMPER.git
cd CAMPER/CAMPER_DRAMKit
conda create --name CAMPER -f ./environment.yaml
```

In both cases, you can activate the newly made environments with the command:

```
conda activate CAMPER
```

Provided all things have gone smoothly, you will be able to activate this environment at any time and use any of the commands outlined in the Usage section below. If there are any problems, please open an issue in the GitHub repository.

#### With pip

If you are not able to use Conda, you can still install CAMPER_DRAMKit with pip using the command below. Note that first you will need to install manually install [scikit-bio](http://scikit-bio.org/), and [MMseqs2](https://github.com/soedinglab/mmseqs2), as these tools can't be installed with the other pip dependencies.

```
pip install camper_dramkit
```

#### Installing With DRAM

If you intend to use CAMPER_DRAMKit with DRAM it may be expedient to install them in the same Conda environment. This is easy to do if you have already made a DRAM Conda environment with the [instructions in the README](https://github.com/WrightonLabCSU/DRAM) then you can add CAMPER with the following commands:

```
wget https://github.com/WrightonLabCSU/CAMPER/main/CAMPER_DRAMKit/environment.yaml
conda env update --name DRAM -f ./environment.yaml
```

#### Note

If you install CAMPER_DRAMKit you will get the latest version of the CAMPER database with it. If you want more control over the database, you can override the default data with the instructions in [Other Tools and flags](https://github.com/WrightonLabCSU/CAMPER#other-tools-and-flags).
## Usage

Once installed, CAMPER_DRAMKit will provide 3 commands, `camper_annotate`, `camper_distill` and `combine_annotations_lowmem`.  These commands along-side DRAM can enable a variety of workflows.

### Standalone Workflow

The simplest workflow is the annotation of a single gene faa file. An example of such a workflow is shown below.

```
camper_annotate -i my_genes.faa -o my_output
camper_distill  -i my_output/annotations.tsv -o my_output/distillate.tsv
```

These commands will make 2 files and the output folder. In the command above I named it `my_output` but you can give it any name appropriate to your project. The output folder gives a place for your files to live and provides a place for temporary files to be stored when using the annotation file. The user will be most interested in the data produced by the program, aka the raw annotations file, and the distillate.

Raw annotations: The `camper_annotate` command makes the output folder and the raw annotations file within it, always named annotations.tsv. With this tab separated file we are trying to provide a consistent, thoughtful, and user-friendly version of the annotations against the CAMPER database. The annotations.tsv has the following columns.
An unnamed index, with the gene ID from the faa file.
  - `fasta` holding the name of the faa file itself.
  - `camper_hits`, A longer CAMPER ID where applicable
  - `camper_rank`, A match quality rank based on the value of the bit score. This is determined by the  .
  - `camper_bitScore`, The bit score from the best search result.
  - `camper_id`: Unique CAMPER ID used in the distillation step.
  - `camper_definition`: A short description of the CAMPER match in the database.
  - `camper_search_type`: Tells you if a HMM profile or blast search found this match.
This script can also be used with a previously run DRAM annotations  file in order to add the CAMPER related fields, above, to the pre-existing annotations. See the next section for more details. If you don’t have DRAM and you need to call your genes [Prodigal](https://github.com/hyattpd/Prodigal) is the suggested tool.

Note also that the gene ids may be slightly different in the annotations file made without a DRAM annotations file to match to, as the verbatim names from the faa file are used instead DRAMs formatted names.

Distillate: Once you have the annotations file, you can use `camper_distill` in order to make a mini distillate. This distillate is a single tab separated file that can be used in the same way as the metabolism summary. It provides a wealth of information like detailed gene descriptions for all CAMPER genes, and a count of hits to each gene in the CAMPER DB. The `camper_distill` script can also leverage dbCAN ids and KO ids from KOfam or KEGG in order to provide more insight into polyphenol metabolism.

### DRAM Combination Workflow

As previously stated, DRAM1.4.0 will include CAMPER by default as an easy tool, but it is possible to use CAMPER with any version of DRAM after 1.3, with one additional command. First follow the instructions above to update your DRAM environment with CAMPER_DRAMKit. Then with that environment activated, you should be able to run the following commands to make a new raw annotations file with all the DRAM data you expect and the CAMPER data added in.

If you are not able to update your DRAM environment for whatever reason, you will simply need to switch environments mid-workflow.

```
DRAM.py annotate -i 'my_bins/*.fa' -o dram_output
camper_annotate -i my_genes.faa -a dram_output_annotation -o camper_dram_output
```

This will create a new set of raw annotations with CAMPER data added, in this case the path of the new file will be `camper_dram_output/annotations.tsv`. You now have two options, if you are only interested in polyphenol metabolism the fastest way to get a distilled summary of related genes is to use the `camper_distill` command to get a distillate with all the key genes from both DRAM and CAMPER data.

```
camper_distill  -i camper_dram_output/annotations.tsv -o camper_dram_output/camper_distillate.tsv
```
### Other Tools and flags

In order to further customize your workflow, you can take advantage of a few more options in the CAMPER_DRAMKit package.
**Combine Annotations With low memory:** You may want to reannotate many annotation files, possibly from more than one version of DRAM. To this end we include the `combine_annotations_lowmem` which should combine many annotations quickly and with a small memory footprint, even if they come from different versions of DRAM. The command is used like so:
```
Combine_annotations_lowmem -i /path/to/many/dramfolders/*/annotations.tsv -o combined_annotation.tsv
```
The input path needs to be a wild card pointing to a set of DRAM annotation files, this is passed to the python glob command, but the format should be familiar to anyone who uses bash and you can test it with the `ls` command.
**Manually Specifying the Location of CAMPER Files:** The behavior of the `camper_annotate` and `camper_distill` commands is controlled by the latest version of the CAMPER dataset. If you want to use an older version of CAMPER, it is suggested you install the older version of the CAMPER_DRAMKit tool, as they will be released together and be mutually compatible. However, if you must, you can also specify the files to use with  `camper_annotate` and `camper_distill` using the appropriate arguments. An example is shown below.
```
camper_annotate -i my_genes.faa -o my_output \
    --camper_fa_db_loc CAMPER_blast.fa \
    --camper_fa_db_cutoffs_loc CAMPER.hmm \
    --camper_hmm_loc CAMPER_blast_scores.tsv  \
   --camper_hmm_cutoffs_loc CAMPER_hmm_scores.tsv
camper_distill  -i my_output/annotations.tsv -o my_output/distillate.tsv \
    --camper_distillate CAMPER_distillate.tsv
```
If at any time you forget these arguments, remember that running any script with the `--help` flag will provide more information. Also note that if you do not specify one or more arguments, the default data will be used.
