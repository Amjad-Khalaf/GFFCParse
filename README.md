# GFFCParse

<p align="center">

<img src ="https://github.com/user-attachments/assets/b6aa0d8c-3f8f-49b3-8805-da0177e5d9fc" width = "500">

</p>

# Table of Contents

- [GFFCParse](#gffcparse)
- [About](#about)
- [Installation](#installation)
- [Usage](#usage)
- [Text File Output](#text-file-output)
- [Plot Output](#plot-output)
- [Miscellaneous](#miscellaneous)

## About

<div align="justify">
  
**GFF-C**<sub>hinchilla</sub>**-Parse** is a tool written in C (with the plotting carried out in R) to parse a [General Feature Format (GFF) file](https://www.ensembl.org/info/website/upload/gff.html), and: 

1. Summarise it.
2. Extract intron, exon & mRNA sizes, locations, and counts into TSV-formatted text files.
3. Make png and svg summary plots.

</div>

## Installation

<div align="justify">

To run GFFCParse, you will need an environment with the following installed:

</div>

```
ggplot2==3.5.1
dplyr==1.1.4
scales==1.3.0
```
<div align="justify">
  
To install GFFCParse:

</div>

```
git clone https://github.com/Amjad-Khalaf/GFFCParse.git
cd GFFCParse
make
export PATH=$PATH:path_to_GFFCParse_directory #The compiled C binary will call "plot.r", so make sure it is found in the path and can be called without the abolute path
```

## Usage

<div align="justify">
  
To run GFFCParse, all you need to do it provide it with a tab-separated GFF file and an output prefix of your choice.

</div>

```
GFFCParse <gff_file> <output_prefix>
```

<div align="justify">

Once run, GFFCParse will produce 4 text files and 4 plots (two copies of each, one in png and one in svg).

</div>

## Text File Output

### 1. <output_prefix>_**mRNA**.tsv

| Gene      | Num Copies | Gene Start | Gene End | Gene Len | Num mRNAs | mRNA                  | mRNA Start | mRNA End | mRNA Len | Num Exons | Sum Exon Len | Avg Exon Len | Max Exon Len | Num Introns | Sum Intron Len | Avg Intron Len | Max Intron Len |
|-----------|------------|------------|----------|----------|-----------|-----------------------|------------|----------|----------|-----------|--------------|--------------|--------------|-------------|----------------|----------------|----------------|
| 52at6029  | 4          | 12885      | 15212    | 2327     | 1         | 52at6029_mRNA_1       | 12885      | 15212    | 2327     | 1         | 2327         | 2327         | 2327         | 0           | 0              | 0              | 0              |
| 592at6029 | 6          | 39244      | 40320    | 1076     | 1         | 592at6029_mRNA_1      | 39244      | 40320    | 1076     | 1         | 1076         | 1076         | 1076         | 0           | 0              | 0              | 0              |
| 978at6029 | 8          | 13687      | 15240    | 1553     | 1         | 978at6029_mRNA_1      | 13687      | 15240    | 1553     | 1         | 1553         | 1553         | 1553         | 0           | 0              | 0              | 0              |

   
### 2. <output_prefix>_**exon**.tsv

| Gene      | Gene Start | Gene End | Gene Len | mRNA                  | mRNA Start | mRNA End | mRNA Len | Exon                  | Exon Start | Exon End | Exon Len |
|-----------|------------|----------|----------|-----------------------|------------|----------|----------|-----------------------|------------|----------|----------|
| 52at6029  | 12885      | 15212    | 2327     | 52at6029_mRNA_1       | 12885      | 15212    | 2327     | 52at6029_mRNA_1_exon_1 | 12885      | 15212    | 2327     |
| 592at6029 | 39244      | 40320    | 1076     | 592at6029_mRNA_1      | 39244      | 40320    | 1076     | 592at6029_mRNA_1_exon_1 | 39244      | 40320    | 1076     |
| 978at6029 | 13687      | 15240    | 1553     | 978at6029_mRNA_1      | 13687      | 15240    | 1553     | 978at6029_mRNA_1_exon_1 | 13687      | 15240    | 1553     |
  
### 3. <output_prefix>_**intron**.tsv

| Gene      | Gene Start | Gene End | Gene Len | mRNA                  | mRNA Start | mRNA End | mRNA Len | Intron                  | Intron Start | Intron End | Intron Len |
|-----------|------------|----------|----------|-----------------------|------------|----------|----------|-------------------------|--------------|------------|------------|
| 3507at6029| 221708     | 224246   | 2538     | 3507at6029_mRNA_1      | 221708     | 224246   | 2538     | 3507at6029_mRNA_1_intron_1 | 221753       | 224042     | 2289       |
| 2773at6029| 26803      | 42583    | 15780    | 2773at6029_mRNA_1      | 26803      | 42583    | 15780    | 2773at6029_mRNA_1_intron_1 | 26911        | 42406      | 15495      |
| 2773at6029| 10723      | 26627    | 15904    | 2773at6029_mRNA_1      | 10723      | 26627    | 15904    | 2773at6029_mRNA_1_intron_1 | 10966        | 26531      | 15565      |

   
### 4. <output_prefix>_**summary**.tbl

```
GFF SUMMARY
=================================================================
Description                                        Count / Length / Percentage
-----------------------------------------------------------------
Total number of genes:                                    179
Total number of all isoforms:                             179
Percentage of single-isoform genes:                    100.00%
Percentage of multi-isoform genes:                       0.00%
Percentage of mono-exonic isoforms (intronless):        89.39%
Percentage of multi-exonic isoforms (with introns):      10.61%
-----------------------------------------------------------------
Total length of all isoforms:                          370356bp
Total number of exons detected:                           206
Total length of exons detected:                        166903bp
Total number of introns detected:                          27
Total length of introns detected:                      203399bp
Total number of all mono-exonic isoforms:                 160
Total length of all mono-exonic isoforms:              149555bp
Total number of all multi-exonic isoforms:                 19
Total length of all multi-exonic isoforms:             220801bp
=================================================================
Single-Copy Summary
-----------------------------------------------------------------
Total number of single-copy isoforms:                      18
Total length of single-copy isoforms:                   38280bp
Total number of exons in single-copy isoforms:             30
Total length of exons in single-copy isoforms:          14352bp
Total number of introns in single-copy isoforms:           12
Total length of introns in single-copy isoforms:        23904bp
Total number of mono-exonic single-copy isoforms:          11
Total length of mono-exonic single-copy isoforms:        7912bp
Total number of multi-exonic single-copy isoforms:          7
Total length of multi-exonic single-copy isoforms:      30368bp
=================================================================
Multi-Copy Summary
-----------------------------------------------------------------
Total number of multi-copy isoforms:                      161
Total length of multi-copy isoforms:                   332076bp
Total number of exons in multi-copy isoforms:             176
Total length of exons in multi-copy isoforms:          152551bp
Total number of introns in multi-copy isoforms:            15
Total length of introns in multi-copy isoforms:        179495bp
Total number of mono-exonic multi-copy isoforms:          149
Total length of mono-exonic multi-copy isoforms:       141643bp
Total number of multi-exonic multi-copy isoforms:          12
Total length of multi-exonic multi-copy isoforms:      190433bp
=================================================================
```

## Plot Output

### 1. **mRNA count**

<p align="center">

<img src = "https://github.com/user-attachments/assets/314eb00e-3c2d-4d5d-b42b-d2e702d1b207" width = 600>

</p>

### 2. **mRNA length**

<p align="center">

<img src = "https://github.com/user-attachments/assets/35f32137-6826-46b4-8baa-b3a2a155497a" width = 600>

</p>

### 3. **Exon count**

<p align="center">

<img src = "https://github.com/user-attachments/assets/805bd270-e036-4259-80b4-e006de44acdd" width = 600>

</p>

### 4. **Intron count**

<p align="center">

<img src = "https://github.com/user-attachments/assets/d09392f9-e54b-4c81-9b50-e19e91088944" width = 600>

</p>

## Miscellaneous

<div align="justify">

1. Plot.r may not work for your specific gff file (if for example you have a huge number of genes) or would like to analyse different categories, so please feel free to edit it!
  
2. With **BUSCO output**, a gene's gff may not necessarily reflect the final reported genes in the .fna and the .faa (there is a fitering step after the individual gffs' are made). To get a gff which reflects the final called genes, you can do this:

</div>

```
grep -Ff <(grep '^>' 433at6029.fna | sed 's/>//; s/:.*//') 433at6029.gff
```
