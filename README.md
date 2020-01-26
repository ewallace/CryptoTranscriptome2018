# CryptoTranscriptome2018

Cryptococcus neoformans transcriptome and translatome analysis 2018-2020.

Edward Wallace (Edward.Wallace@ed.ac.uk), Corinne Maufrais, Guilhem Janbon. 

This work is the analysis, including most figure construction, for the paper:

> Quantitative global studies reveal differential translational control by start codon context across the fungal kingdom
> Wallace, Edward; Maufrais, Corinne ; Sales-Lee, Jade; Tuck, Laura; de Oliveira, Luciana; Feuerbach, Frank; Moyrand, Frédérique; Natarajan, Prashanthi; Madhani, Hiten; Janbon, Guilhem.
> Nucleic Acids Research, 2020 (in press)

A preprint of this work is on biorxiv: https://doi.org/10.1101/654046

**Please cite the paper if you refer to or re-use any of this code.**

This repository contains the code used for figures in the paper. It also contains additional figures that didn't make it into the paper, f incomplete or unsuccessful analyses, and so on. We are sharing it here as-is because it is better to work open.

Code is reusable under an Apache 2.0 Licence https://www.apache.org/licenses/LICENSE-2.0

# Contents

## ATGContext

ATG Context and translation analysis in Cryptococcus neoformans H99 and Cryptococcus deneoformans JEC21.

This contains R markdown files .Rmd format, making cryptococcus-centric figures in the paper.

## ATGContextOthers

ATG Context and translation analysis in other model fungi: Saccharomyces cerevisiae, Candida albicans, Schizosaccharomyces pombe, Neurospora crassa.

This contains R markdown files .Rmd format, making species-comparison files in the paper.

## CDSFasta

Coding sequence fasta files for Cryptococcus

## CDSFastaOthers	

Coding sequence fasta files for other model fungi

## CnGFF	

Annotation/ genome feature files gff3 for Cryptococcus

## GFFOthers

Annotation/ genome feature files gff3 for other model fungi

## DiffExp	

DESeq2 differential expression data on delta-upf1 vs wildtype JEC21.

## FigSVG

Some svg-format files of sequence logos etc.

## ProteinFasta	

Protein sequence fasta files for Cryptococcus

## ProteinFastaOthers	

Protein sequence fasta files for other model fungi

## SharedScripts

SharedScripts/sharedScripts.R contains R functions that are used in analyses, collected together into a single source file.


## CnRiboCode	

Code to run RiboCode for ORF analysis on Cryptococcus. Not completely successful and didn't reach the paper.

## uORF

Upstream ORF analysis, both sequence fasta and genome feature files. Not completely successful and didn't reach the paper.
