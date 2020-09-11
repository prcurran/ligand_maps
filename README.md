# Ligand Maps: Generating a Pharmacophoric Models from Superimposed PDB Data.


This workflow was developed during my PhD and was used as a reference to compare hotspot pharmacophore models (generated with the Hotspots API) against.

There are 4 main steps to the workflow: 
1. Downloading data
2. Clustering ligands
3. Superimposing selected structures
4. Pharmacophore generation


## Download the PDB data

Using the new RCSB web services, all structures for a given UniProt ID and structure resolution < 3.0. 

## Clustering the Ligands

The central assumption in this approach is that is that more abundant pharmacophoric features are more important for ligand binding. However, often chemical scaffolds are overrepresented in PDB collections (many entries with the putative ligand (ATP in the case of CDK2) or perhaps many compounds from a single lead series published on a crystallography enabled project).  Therefore, care must be taken to ensure features aren't detected as a consequence of an overrepresented scaffolds but rather due to chemically diverse ligands utilising common intermolecule interactions. 

To cluster the detected ligands, I used fingerprints generated with RDKit and sklearn TSNE dimensionality reduction with hdbscan to cluster.

`CDK2.png` visualises the compounds in 2D space, the points are coloured by cluster ID and the 'X' denote the selected representatives.

## Superimposition

The selected proteins are superimposed with the CSD Python API. Although the BioPython based superimposer built in part 1 of this task works well, I have found mixing PDB file parsers can be problematic so I have used the CSD Python API for this demo.

## Pharmacophore Generation

Pharmacophore features are detected using CrossMiner feature definitions. A concensus pharmacophore is found by adding feature points to a grid and detected local maxima.

Run `pymol_file.py` in PyMOL to see the out visualisation.

