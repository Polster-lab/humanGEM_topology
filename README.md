# ihuman topology

MATLAB script for extracting the topology of the connectivity between molecular components in humanGEM. As a result the folowing files are produced:

| File name | Description|
|----------|----------|
| ihuman_rxns.txt | File listing all reaction IDs, reaction names, formulas and encoding genes (grRules) |
| ihuman_mets.txt | File listing all metabolites IDs, metabolite names, and chemical formulas |
| ihuman_subSystems.txt| File listing all metabolic subsystems (pathways) represented in humanGEM |
| ihuman_rxnGeneMatrix.txt | Incidence matrix connecting reactions in the model to the genes that encode for their enzmes|
| ihuman_genesMetMatrix.txt | Incidence matrix connecting genes in the model to the metabolites present in their corresponding reactions |
| ihuman_geneSubSystemsMatrix.txt | Incidence matrix connecting genes in the model to metabolic subsystems (presence or their corresponding reactions in a given pathway) |
| ihuman_metSubSstemsMatrix.txt| Incidence matrix connecting metabolites to metabolic subsystems in the model |

An additional text file storing the model version and the date of last results update is also saved for version control.

## Required Software
- MATLAB 2020 (or later)
- The RAVEN toolbox (see installation instructions [here](https://github.com/SysBioChalmers/RAVEN/wiki/Installation#installation))

## Installation
Clone master branch from Polster-Lab GitHub.

## Instructions
To update results in this repository, set the subfolder 'ihuman_topology/code' in this repo as a working directory in MATLAB, then type
´extractModelTopology´ in the MATLAB command window. Results are stored as .txt files in the subfolder 'ihuman_topology/results'.

Polster Lab, Division of Systems and Synthetic Biology, Chalmers University of Technology.

Last update: 29th of May 2024.