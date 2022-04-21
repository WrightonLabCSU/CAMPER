# CAMPER

## Abstract
A prevailing paradigm in arctic peatlands these ecosystems is that polyphenols and phenolic compounds contribute to soil carbon sequestrationcarbon storage by binding microbial extracellular enzymes, limiting soil carbon decomposition and carbon dioxide (CO2) production. In this model, coined the enzyme latch, polyphenol degradation is controlled by a single oxygen-requiring enzyme known as a polyphenol oxidase. We hypothesized this model oversimplifies oversimplified both polyphenol chemistry and microbial metabolism and persists due to an absence of high-resolution, molecular evidence. Using a permafrost peatland thaw gradient in Sweden, we coupled field multi-omics with soil geochemical data to audit (poly)phenolic chemistry and microbial gene expression in the thawed palsa, bog, and fen habitats.  Using both traditional and ‘omics methods, oOur data failed to support the enzyme latch using both traditional and ‘omics methods, and instead showed showing positive relationships between soil polyphenolss, microbial activity, and CO2 concentrations. To uncover the microbial processes contributing to these observations, we developed an an annotation tool for (poly)phenol degradation genes and applied it to peat metatranscriptomes. Beyond polyphenol oxidase, several expressed enzymes genes could modulated polyphenol concentrations across peatland habitatshabitats and depths. Microorganisms from Microbial genomes from 10 phyla expressed genes for transforming polyphenols or their degradation products, with the majority assigned to novel lineages within the Acidobacteriota, Actinobacteriota, and Chloroflexota. Paired to metabolite data, we positPaired metabolite data indicated that flavonoids, lignans, and phenolic acids are important (poly)phenolic carbon sources, that when  in this peatlanddecomposed that can lead to CO2 production, rather than constrain it. This chemical and biological framework defines a previously ignored aspect oenigmatic aspectf peatland soil microbial metabolismcarbon cycle, revealing a a prevalent metabolic and complex decomposition process network active in aArctic soilscarbon cycles.

<img width="741" alt="Screen Shot 2022-03-07 at 1 47 10 PM" src="https://user-images.githubusercontent.com/95941779/157345312-27679138-c32c-4e76-8923-a2c776bccbe9.png">


## CAMPER DATA



## Using CAMPER as part of DRAM

To facilitate the use of the CAMPERS data set we integrating it into [DRAM](https://github.com/WrightonLabCSU/DRAM).

With the releas of DRAM1.4.0 you will be able to use the CAMPERS data set without any aditional data sets simply by suplying the `--use_campers` flag useing the annotate comand, like so.

```
DRAM.py annotate --use_campers -i 'my_bins/*.fa' -o annotation
```

When you run `DRAM.py distill` it should detect that campers data was included and you will find a campers tab in your metaboic summary.  

If you want to use campers as part of your annotation pipline along with an older version of DRAM see the (Whith DRAM) section bellow.

## Campers Standalone Tool: CAMPER_DRAMKit

In adition to the CAMPER module in DRAM, we have also created a stand alone tool, named CAMPER_DRAMKit that can be used by its self or in conjunction with older versions of DRAM. Infact CAMPER_DRAMKit is realy a much smaller verson of DRAM that falows much the same work flow as dram and has simulare capabilitys. The set up will be a bit difrent for each use case so you may want to revew the options before you begine. If you are wondering what this tool can do for you you can jump to the Useage section and get a beter idea for how the tool works.

#### Setup With Conda

If you are only interestd in CAMPER, or don't want to wait for the DRAM1.4 release you can use the stand alone tool included in this repository. The simplest whay to get started is too download this repository and use conda to install CAMPER_DRAMKit. Using the firts set of comands below. You will need all the campers data to be downloaded so if you don't wan't to use git you will need to download the data with the link above. If you already have the data and just want to install or upgrade the packege the second set of comands will work for you. 

```
git clone https://github.com/WrightonLabCSU/CAMPER.git
cd CAMPER/CAMPER_DRAMKit
conda create --name CAMPER -f ./environment.yaml
```

or if you have already downloaded the campers data,

```
wget https://github.com/WrightonLabCSU/CAMPER/main/CAMPER_DRAMKit/environment.yaml
conda create --name CAMPER -f ./environment.yaml
```

In both cases you can activate the newly made enviroments with the comand:

```
conda activate CAMPER
```

Provided all things have gone smothly you will be able to activate this enviroment at any time and use any of the comands outlined in the usige section below. If there are any problems pleas open an ishue in the github repository.

#### With pip

If you are not able to use conda you can still install CAMPER_DRAMKit with pip using the comand below. Note that first you will need to install manualy install (scikit-bio)[http://scikit-bio.org/], and (MMseqs2)[https://github.com/soedinglab/mmseqs2], and these tools cant be install with the other pip depenacies. 

```
pip install camper_dramkit
```

#### Installing With DRAM

If you intend to use CAMPER_DRAMKit with DRAM it may be expediant to install them in the same conda enviroment. This is easy to do if you have already made a DRAM Conda enviroment with the [instructions in the README](https://github.com/WrightonLabCSU/DRAM) then you can add CAMPERS with the falowing comands:

```
wget https://github.com/WrightonLabCSU/CAMPER/main/CAMPER_DRAMKit/environment.yaml
conda env update --name DRAM -f ./environment.yaml
```

## Usage

Once installed it is easy to use camper dram
