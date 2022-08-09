# ESSdb

This repository contains programs for the creation of a Database of **E**nzymatic **S**tep
**S**equences (**ESS**s). An Enzymatic Step is represented by the Enzyme Commission number (EC number)
that describes a given reaction. These scripts creates 2 versions of ESSs, one using the 4 levels (ess4)
of each EC number and the other using only the first 3 levels. The first 3 levels can describe
the type of reaction of an enzyme.

The Database created can be used for the comparison of of the ESS using the software
[alignESS](https://github.com/acpooth/alignESS).

## Dependencies.

This version of ESSdb is written in Python 2.7 and work with the latest version of the following
packages:

 - sqlite3
 - numpy
 - matplotlib
 - networkx
 
The preferred way to get the dependencies is creating a conda environment as shown here. It is
necessary to install anaconda, miniconda or mambaforge distributions.

### Conda environment.

The conda environment can be created with the following command:

```bash
conda create -n essdbpy2 -c conda-forge python=2 sqlite numpy matplotlib  networkx
```

Also the *.yml* file included in this repository can be used to create the environment, as follows:

```bash
conda env create -f essdbpy2.yml
```

The *.yml* file contains the working and tested version of all dependencies.


## Folder structure

The program for the creation of the Database of ESS requires a specific directory structure
to properly work.



The initial directory structure contains 3 main directories:
 - **Scripts**: All the scripts are stored here.
 - **Enzyme**: The enzyme (EC number) to gene relation lists for each species is stored here.
	 Each file must be named using the KEGG letter code for each organism.
 - **Maps**: The metabolic map KGML files are stored here. This directory must contain
	 a directory for each species represented in `Enzymes/` directory using the same KEGG 
	 letter code for organisms. The KGML files for each organism must be inside the corresponding
	 sub directory.
	 **IMPORTANT**: Only KEGG maps representing metabolism are allowed by the program. These maps
	 are those withe codes bellow 01000. For example eco00010 or hsa00540

The following is an example of the initial directory structure before the creation of the database.
Please see the [KEGG REST API](https://www.kegg.jp/kegg/rest/keggapi.html) for examples of how to 
obtain the files (list files and KGML) necessary for this program.

```
essdb
│
├── Enzymes
│   ├── eco.list
│   ...
│   └── hsa.list
├── Maps
│   ├── eco
│   │   ├── eco00010.xml
│   │   ...
│   │   └── eco00020.xml
│   └── hsa
│       ├── hsa00010.xml
│       ...
│       └── hsa00540.xml
└── Scripts
    ├── kegg2seq.py
    ├── make_seq.py
    ├── make_seq.pyc
    ├── metabolismTypeAss.py
    ├── noRedundantDB.py
    ├── parse_kgml.py
    ├── parse_kgml.pyc
    ├── traducir.py
    └── traducir.pyc
```

Final folder structure after execution. The *../Db* folder is created and within, the
*seq.db* file is created. This file is a *SQLite* database file.

```
essdb
│
├── Db
│   └── seqs.db
├── Enzymes
│   ├── eco.list
│   ...
│   └── hsa.list
├── Maps
│   ├── eco
│   │   ├── eco00010.xml
│   │   ...
│   │   └── eco00020.xml
│   └── hsa
│       ├── hsa00010.xml
│       ...
│       └── hsa00540.xml
└── Scripts
    ├── kegg2seq.py
    ├── make_seq.py
    ├── make_seq.pyc
    ├── metabolismTypeAss.py
    ├── noRedundantDB.py
    ├── parse_kgml.py
    ├── parse_kgml.pyc
    ├── traducir.py
    └── traducir.pyc
```

## Execution.

**The scripts must be executed from ./Scripts folder**

1. Get enzyme list files for each organism.
2. Get Metaboilic maps for each organism.
3. Enter the ./Scripts folder.
4. Activate conda environment.

```bash
conda activate essdbpy2
```

5. Execute kegg2seq.py

```bash
python kegg2seq.py
```

6. Create non redundant table of ESS.

```bash
python metabolismTypeAss.py ../Db/seqs.db
```

7. Add a metabolism type column to seqs table

```bash
python noRedundantDB.py ../Db/seqs.db
```

## Database structure


## Important.

- Only creates the ESS from those organisms that has a Enzyme list in ./Enzymes folder
  - Those organisms with enzyme list without forlder in ./Maps will raise an error.

## Papers.
The programs presented here were used all or in parts in the following papers.

1. Comparison of Metabolic Pathways in Escherichia coli by Using Genetic Algorithms
P Ortegon, AC Poot-Hernández, E Perez-Rueda, K Rodriguez-Vazquez
Computational and structural biotechnology journal 13, 277-285. 2015.

2. The alignment of enzymatic steps reveals similar metabolic pathways and probable recruitment events in Gammaproteobacteria
AC Poot-Hernandez, K Rodriguez-Vazquez, E Perez-Rueda
BMC genomics 16 (1), 957. 2015.

3. Identification of functional signatures in the metabolism of the three cellular domains of life
P Escobar-Turriza, R Hernandez-Guerrero, AC Poot-Hernández, ...
PloS one 14 (5). 2019.


