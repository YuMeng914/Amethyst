# Amethyst
Amethyst is a pipeline for analysis of Deep Mutational Scanning (DMS) sequencing data. This release features a GUI which allows users to easily install and operate, and generates a heatmap which visualizes the enrichment of each variant in the DMS library. Also by integrating NGmerge (Gaspar JM, 2018) into this pipeline, users can directly upload the paired-end sequencing data acquired from NGS. The enrichment score for each variant is calculated as the logfold change of its percentage in experiment group compared to that in control group. Variants with read counts below 20 are filtered out to reduce bias. <br>

## Installation
Amethyst v1.0.0 is a software designed for MacOS. Download AmethystPipeline.dmg from release, open it, and drag the Amethyst app icon to **Applications** folder.<br>
<br>
Note that this software is developed for macOS but is not signed or notarized by Apple. Hence, when running the software for the first time, macOS may display a security warning. This is expected for open-source research software. To allow execution:<br>
1. Open **System Settings → Privacy & Security**<br>
2. Scroll down to **Security**<br>
3. Click **"Allow Anyway"**<br>
4. Re-run the software<br>

## Usage
Amethyst requires several parameters for the analysis process:<br>
<br>
**Fixed sequence (DNA)**: This should be the sequence of usually 7 to 10 nucleotides immediately before your mutagenesis ORF. For example, if your DMS library involves Site Saturation Mutagenesis (SSM) on amino acids 101-200 on protein X, then you can input the wild-type DNA sequence of amino acids 98-100. **DNA sequence should be in lower case.** Note that fixed sequence should be singularly present in the sequencing reads.<br>
<br>
**AA length**: This is the length of your mutagenesis ORF. For example, if your DMS library involves Site Saturation Mutagenesis (SSM) on amino acids 101-200 on protein X, you should put 100. <br>
<br>
**Reference Sequence (Protein)**: This is the wild-type protein sequence of your mutagenesis ORF. **Protein sequence should be in upper case.** Note that you should only input the sequence of the segment in which the mutagenesis is introduced rather than the whole protein.<br>
<br>
**Input Files**: These should be the paired-end sequencing files before merging, usually acquired directly from NGS. This pipeline aims to identify the enrichment of each variant between the control group and the experiment group. For example, if you are screening for Gain-Of-Function (GOF) variants, the control group could be the entirety of cells before sorting or the low bin collected during sorting, and the experiment group should be the high bin collected during sorting. The files can be in either .fastq or .fastq.gz format. <br>
<br>
**Output PDF Filename**: Input the desired name for the heatmap that is to be generated. It will be in the same folder as the input files.<br>

## Citations
This release integrates the following third-party tools:<br>
**NGmerge**: Gaspar JM. BMC Bioinformatics. 2018 Dec 20;19(1):536.

## License
Copyright (c) 2025 Yu Meng (yu.meng@epfl.ch)<br>
<br>
This project is licensed under the MIT License.<br>
See the LICENSE file for details.<br>
