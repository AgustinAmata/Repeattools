# Repeattools: analysis of transposable elements across species
## Introduction
Transposable elements (TEs) or transposons are DNA mobile elements with a great
ability to move across the genome and intervene in its structure and evolution, representing a
source of genetic variability. Several environmental and genetic factors regulate their
diversification and activity, so these structures can be useful to study evolutionary relationships
between species. There exist different types of software to identify and classify TEs, such as
RepeatModeler or TESorter. Programs such as RepeatMasker can recognize the location of the
TEs in the genome. Nonetheless, there is an important lack of tools that allow to analyze the
diversity of TEs between species, compare their transposon profiles and analyze the divergence
of different copies of a given TE in different genomes.

Repeattools is a package created in Python that tries to address this problem.
It is designed to perform statistical analysis and visualization on the output of the programs
RepeatMasker and TESorter. Simply put, Repeattools processes the transposable element landscape
of the selected species and provides several forms of visualization. Specifically, it can output the
following:

- A heatmap comparing the ET profile of the different species. For a given species, this profile is composed of the number
of reads for each type of ET found in its genome.
- A PCA that also provides a form of comparing the ET profile.
- Box plots that compare the divergence distribution of a given ET across several species and groups.
- Violin plots that compare the divergence distribution of all ETs in a certain [classification level](##classification-system-of-Repeattools)

## Dependencies and installation
### Dependencies
Dependencies for Repeattools are already specified in the [setup.py](./setup.py) of the repository. As a remainder:
- Matplotlib version 3.7.3
- Pandas version 2.0.3
- Numpy version 1.25.1
- Scipy version 1.11.1
- Scikit-learn version 1.3.0
- Seaborn version 0.12.2
- PyQt5 version 5.15.9
- Ete3 version 3.1.3

## Classification system of Repeattools

