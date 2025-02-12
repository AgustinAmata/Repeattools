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

- A heatmap comparing the TE profile of the different species. For a given species, this profile is composed of the number
of reads for each type of TE found in its genome.
- A PCA that also provides a form of comparing the TE profile.
- Box plots that compare the divergence distribution of a given TE across several species and groups.
- Violin plots that compare the divergence distribution of all TEs in a certain [classification level](#classification-system-of-Repeattools)

## Dependencies and installation
### Dependencies
Dependencies for Repeattools are already specified in the [setup.py](./setup.py) of the repository. As a remainder:

+ **At least Python version 3.10.11 (not tested on Python 3.12 onwards, so I strongly suggest using
Python 3.10 or 3.11)**

| Third-party library | Version |
|---------------------|---------|
| Matplotlib          | 3.7.3   |
| Pandas              | 2.0.3   |
| Numpy               | 1.25.1  |
| Scipy               | 1.11.1  |
| Scikit-learn        | 1.3.0   |
| Seaborn             | 0.12.2  |
| PyQt5               | 5.15.9  |
| Ete3                | 3.1.3   |

### Installation

The package can be installed in a virtual environment using [venv](https://docs.python.org/3.10/library/venv.html)

To create a virtual environment, use the following command:
```
$ python -m venv path/to/myenv

#To activate the virtual environment
#POSIX
$ source <myenv>/bin/activate

#Windows
$ <myenv>\Scripts\activate

#To deactivate the virtual environment
$ <myenv> deactivate
```

Once the virtual environment is set up and **activated**, proceed to download the repository in your preferred directory and execute setup.py:
```
$ git clone https://github.com/AgustinAmata/Repeattools.git

#Go to the Repeattools directory
$ cd path/to/Repeattools
#Install dependencies from setup.py
$ pip install .

# (If there is any change in the repository) Updates and merges your local work with the online repository
$ git pull Repeattools
```

To test that the installation was done correctly, use the following command in the Repeattools directory:
```
$ python -m unittest discover
```

## Classification system of Repeattools
Repeattools uses a hierarchical classification system of four levels: Class, Subclass, Superfamily, and Element.
This classification system was created considering the classification system of RepeatMasker and TESorter.
RECollector makes the conversion to this classification for each TE found in the RepeatMasker and TESorter files
for each species.

Examples:

| RepeatMasker/TESorter class/family column |    Class    |    Subclass    |    Superfamily    |    Element    |
|-------------------------------------------|-------------|----------------|-------------------|---------------|
|                  Unknown                  |   Unknown   |     Unknown    |      Unknown      |    Unknown    |
|                  DNA                      |   Class_II  |     Unknown    |      Unknown      |    Unknown    |
|                  DNA/Maverick             |   Class_II  | DNA_Polymerase |      Maverick     |    Unknown    |
|                  DNA/Casposons            |   Class_II  | DNA_Polymerase |      Casposon     |    Casposon   |

Note: some TESorter detected TEs are converted into a RepeatMasker equivalent before being properly classified. Example:
TIR/hAT is converted into RepeatMasker's DNA/hAT.

## Analysis of RepeatMasker and TESorter with RECollector
RECollector is the first program to be used in the Repeattools package to analyze TEs. It outputs the following:
- A TE count matrix for a certain classification level (that can be specified by the user - default is Superfamily)
in which rows correspond to the different TEs in that level and columns correspond to the diffretent analyzed species.
So, each element in the matrix correspond to the number of copies of a certain TE in a certain species.
- A directory containing the divergence data files for each of the TEs and their copies in the genome of the different
species analyzed in a certain classification level.

### Directory data structure
Before proceding with RECollector, the user must create the directory that is going to be analyzed. It must follow the
following structure:
```
MainSpeciesDir
├── Species1Dir
│   ├── Species1_RepeatMasker_file.out
│   └── Species1_TESorter_file.cls.tsv
└── Species1Dir
    ├── Species2_RepeatMasker_file.out
    └── Species2_TESorter_file.cls.tsv

#The name of the file does not matter, what matters is the file extension, so be sure to use the files from RepeatMasker
#and TESorter that contain the correct extension (.out for RM, and .cls.tsv for TES)
```
## ET profile comparison with REPlotCounts


## ET divergence distribution comparison with REPlotDivergence

