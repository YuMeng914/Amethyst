# Amethyst
Amethyst is a pipeline for analysis of Deep Mutational Scanning (DMS) sequencing data. This release features a GUI which allows users to easily install and operate, and generates a heatmap which visualizes the enrichment of each variant in the DMS library. Also by integrating NGmerge (Gaspar JM, 2018) into this pipeline, users can directly upload the paired-end sequencing data acquired from NGS.

## Installation
Amethyst v1.0.0 is a software designed for MacOS. Download AmethystPipeline.dmg from release, open it, and drag the Amethyst app icon to **Applications** folder. When opening the Amethyst app for the first time, **Right-click** the icon and choose **Open**.

## Usage
Amethyst requires several parameters for the analysis process:<br>
<br>
**Fixed sequence (DNA)**: This should be the sequence of usually 7 to 10 nucleotides immediately before your mutagenesis ORF. For example, if your DMS library involves Site Saturation Mutagenesis (SSM) on amino acids 101-200 on protein X, then you can input the wild-type DNA sequence of amino acids 98-100. **DNA sequence should be in lower case.** Note that fixed sequence should be singularly present in the sequencing reads.<br>
<br>
**AA length**: This is the length of your mutagenesis ORF. For example, if your DMS library involves Site Saturation Mutagenesis (SSM) on amino acids 101-200 on protein X, you should put 100. <br>
<br>
**Reference Sequence (Protein)**: This is the wild-type protein sequence of your mutagenesis ORF. **Protein sequence should be in upper case.** Note that you should only input the sequence of the segment in which the mutagenesis is introduced rather than the whole protein.<br>
<br>

## Citations
This release integrates the following third-party tools:
**NGmerge**: Gaspar JM. BMC Bioinformatics. 2018 Dec 20;19(1):536.
