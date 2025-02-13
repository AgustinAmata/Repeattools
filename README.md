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
following figures:

- A heatmap comparing the TE profile of the different species. For a given species, this profile is composed of the number
of reads for each type of TE found in its genome.
- A PCA that also provides a form of comparing the TE profile.
- Box plots that compare the divergence distribution of a given TE across several species and groups.
- Violin plots that compare the divergence distribution of all TEs in a certain [classification level](#classification-system-of-Repeattools)

Repeattools is composed of three programs: RECollector, REPlotCounts, and REPlotDivergence.
They are explained later in the document.

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
or the [conda package manager](https://docs.conda.io/projects/conda/en/stable/)

To create a virtual environment using venv, use the following commands:
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
for each species. Conversions can be consulted in [config.py](Repeattools/src/config.py).

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
RECollector is the first program to be used in the Repeattools package to analyze TEs.

### Usage
For a complete description of the program and its options use:
```
$ python RECollector.py --help
```

Simple usage of the program is:
```
$ python RECollector.py -i MainSpeciesDir/ -n names_file -o out_file
```

The user can include the `--override` option so that RECollector can
add to the classification information of Unknown copies of RepeatMasker the information retrieved by
TESorter for these same copies. `--override` also affects RepeatMasker TE copies that were only classified
as Class_II if the correspondent information of TESorter matches.

Example:

Before `--override`:
| TE Copy number |    Class    |    Subclass    |    Superfamily    |    Element    | TESorter Class/family |
|----------------|-------------|----------------|-------------------|---------------|-----------------------|
|        1       |   Unknown   |     Unknown    |      Unknown      |    Unknown    |        LINE/L1        |
|        2       |   Class_II  |     Unknown    |      Unknown      |    Unknown    |    Maverick/Unknown   |
|        3       |   Class_II  |     Unknown    |      Unknown      |    Unknown    |     LTR/Retrovirus    |
|        4       |   Class_II  | DNA_Polymerase |      Casposon     |    Casposon   |       TIR/Merlin      |

After `--override`:
| TE Copy number |    Class    |    Subclass    |    Superfamily    |    Element    | TESorter Class/family |                   'overriden'                   |
|----------------|-------------|----------------|-------------------|---------------|-----------------------|-------------------------------------------------|
|        1       |   Class_I   |     LINE       | DFAM-LINE_group_II|    L1         |        LINE/L1        |            Yes (Initial Unknown copy)           |
|        2       |   Class_II  | DNA_Polymerase |      Maverick     |    Unknown    |    Maverick/Unknown   |             Yes (Class_II-compatible)           |
|        3       |   Class_II  |     Unknown    |      Unknown      |    Unknown    |     LTR/Retrovirus    |    No (Class_II not compatible with Class_I)    |
|        4       |   Class_II  | DNA_Polymerase |      Casposon     |    Casposon   |       TIR/Merlin      | No (Copy already classified beyond Class level) |

Additionally, the user can also include the `--depth` option followed by one of the choices included in the program
(e.g., `--depth class`), which can be known by consulting `python RECollector.py --help`.
If `--depth` is not included, RECollector, by default, will create its outputs based on the Superfamily classification level.

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

#Be sure to use the files from RepeatMasker and TESorter that contain the correct extension
#(.out for RM, and .cls.tsv for TES)
```

### Names file structure
Additionally, RECollector uses another tab-separated file to name each of the species for
their representation in REPlotCounts and REPlotDivergence. Rows correspond to the different
species, the first column corresponds to the directory name for that species and the second
column corresponds to the chosen name for that species for its graphic representation.

```
dirname_species_1	species_1_name
dirname_species_2	species_2_name
```
### Outputs
- A TE count matrix for a certain classification level (that can be specified by the user with `--depth`)
in which rows correspond to the different TEs in that level and columns correspond to the diffretent analyzed species.
So, each element in the matrix correspond to the number of copies of a certain TE in a certain species.
- A directory containing the divergence data files for each of the TEs and their copies in the genome of the different
species analyzed in a certain classification level.
- A log file containing information on the command that was run, the process flow of the program, and additional stuff.
In case something went wrong during the data processing, the log file will also contain information of the error that caused
the abortion of the program.
- A tab-separated file meant to be used with REPlotDivergence.

## TE profile comparison with REPlotCounts


### Usage
For a complete description of the program and its options use:
```
$ python REPlotCounts.py --help
```

Simple usage of the program is:
```
$ python REPlotCounts.py -i TECountMatrix.csv -g groups_file -o out_directory/
```

The user can also specify the size of the heatmap with the `--hsize` option followed by two numbers
separated by a space, the first number specifies the width whereas the second specifies the height
(e.g., `--hsize 20 13`). By default, heatmaps are 15 inches wide and 15 inches tall.

### Groups file structure
Additionally, REPlotCounts uses a tab-separated file to group the different species.
Rows correspond to the different species, the first column corresponds to the name of the species
(which must be the same name given to the species when processed with RECollector)
and the second column corresponds to the chosen group for that species.

```
species_1_name	group_1
species_1_name	group_2
species_1_name	group_2
```

### Outputs
- A heatmap (normalized by z-score) that allows to observe the TE profile of species and
compare them with each other and between the same group/different groups. A dendrogram can
be added with `--dendro`.
- A PCA plot that can also help with inferring relationships between species/groups
(percentage of variance explained by both principal components is provided to help with interpretation).
The name of each species can also be plotted with `--names`.

Groups are color-coded.

## TE divergence distribution comparison with REPlotDivergence


### Usage 

For a complete description of the program and its options use:
```
$ python REPlotDivergence.py --help
```

Simple usage of the program is:
```
# If the user only wants to create violin plots
$ python REPlotDivergence.py -v DivergenceDir/ -n RECollector_names_file -o out_directory

#If the user only wants to create box plots
$ python REPlotDivergence.py -b DivergenceDir/<TE_to_analyze>_divergence.csv -g groups_file -o out_directory
```
**Notes: the groups file can be the same one used for REPlotCounts. The user can also create
violin and box plots at the same time if they consider it.**

### Outputs
