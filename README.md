# Project: Ordering the Evolution and Function of Intrinsically Disordered Regions in Neurodegenerative Proteins

Repositories used:
+ BioComputingUP/AlphaFold-disorder
+ GSLBiotech/clustal-omega
+ agemagician/ProtTrans

For dependencies please see requirements.txt

## 1. Calculating Relative Solvent Accessibility (RSA) using Alphafold-disorder

Run the following command to calculate RSA within Alphafold-disorder subdirectory:
```bash
python3 alphafold_disorder.py -i pdbs/ -o out.tsv 
```
Alphafold-disorder outputs two files by default, *out_data.tsv* and *out_pred.tsv*. *out_pred.tsv* is used for IDR prediction as it is the final calculated prediction, *out_data.tsv* is an intermediate calculation (DSSP output). All PDB files used for this project are saved within the 'pdbs' subfolder housed within the alphafold-disorder parent folder. Previous calculations are saved in 'predictions' subfolder.

## 2. Predict Intrinstically Disordered Regions (IDRs) Using Gaussian Filtering

Run the following command to predict IDRs:

```bash
python idr_prediction.py --i path/to/out_pred.tsv --o output.png
```

Script using 1-D Guassian ($\sigma=7$) and 0.581 RSA_threshold criterion to label IDRs within sequence data from `out_pred.tsv`. Plotting logic is hard coded, please edit labeling accordingly.

## 2.1. Benchmarking IDR Prediction

prediction_benchmark.py contains static data from various \sigma value idr_prediction.py iteriations (Human - FUS, TDP-43, and Axin-1). $\sigma$ value indexed by i in protein_species_{i}. Jaccard index calculation implemented in script.

## 3. Generating Multiple Sequence Alignments and Phylogeny Trees using Clustal Omega

Run the following commands within Clustal Omega subdirectory:
*Sequences for FUS, TDP-43, and Axin-1 are saved within the `sequences` subdirectory*

```bash
clustalo -i sequences.fasta -o aligned_sequences.aln --outfmt=clu
clustalo -i sequences.fasta -o phylogeny_tree.dnd --guidetree-out
```

plot_tree.py plots a phylogenetic tree with branch annotations from a Newick format file using Bio.Phylo.

```bash
python plot_tree,py --i path/to/dnd_file --o outpng.png
```
MSA_disorder.py creates MSA graphic using Clustal Omega output (`aligned_sequences.aln`) and IDR ranges of all selected species of that protein *(Can be accessed from sequences subdirectory)*. The path to alignment file and IDR ranges are hard coded.

## 4. ProtTrans-Based Embedding and Clustering for Functional Analysis of IDRs

ProtTrans embeddings, IDR non-conserved residues, and ELM labeling were performed within `prottrans/prottrans.ipynb`. Required sequence and IDR range inputs can be found in the `sequences` subdirectory. Data processing utilized t-SNE, DBSCAN, and enrichment analysis without parameter modifications from the Jupyter notbeook.
