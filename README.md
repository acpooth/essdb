# ESSdb

This repository contains programs for the creation of a Database of **E**nzymatic **S**tep
**S**equences (**ESS**s). An ESS is defined as a set of consecutive enzymatic reactions
in a metabolic network. Each enzymatic reaction, also called enzymatic step is represented
by the Enzyme Commission number (EC number) that describes a given reaction.
These scripts creates 2 versions of ESSs, one using the 4 levels (ess4)
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

This script has the option -g (--graphtype) that specifies the type of graph used
for the construction of the ESS database: 'directed' or 'nodirected'. When the graphs are 
directed, BFS algorithm takes into account the reversibility of the reactions and 
only follows forward edges. This modifies the number and length of the ESSs generated. 
By default the program uses no directed graphs.

8. Add a metabolism type column to seqs table.

```bash
python noRedundantDB.py ../Db/seqs.db
```


7. Create non redundant table of ESS.

```bash
python metabolismTypeAss.py ../Db/seqs.db
```


## Database structure

The database of ESSs contains 3 tables that contains the sequences and relevant information about them. The tables are the
following:

 1. **seqs**: Contains the primary ESSs generated from KEGG *kgml* files. The sequences are stored in different formats, i.e. protein sequences,
 3 level EC ESS, 4 level EC ESS, species and metabolic map of origin, among other information.
 
 2. **map_stats**: Contains statistics about the process of generation of ESSs for each metabolic map.
 
 3. **nrseqs**: Contains the non-redundant set of ESSs obtained from **seqs** table.

The structure of the database is shown in the following figure:

![Structure of ESS database](,/db_structure.png)


### seqs table.

This table is the direct result of the parsing and processing of *kgml* files. The table contains the following columns:

  - **id** : Unique identifier for ESS.
  
  - **gen_seq** : Sequences of genes identifiers. Each enzymatic step  may be represented for one or more gene
  identifiers. If more than one identifier is present, they are separated with a  white space, and may represent an
  enzymatic complex or isozymes. Each enzymatic step is delimited by a '>'. These are the basic sequences from
  whom the other sequences are generated.
  		
  - **ec3**: ESS using the first 3 levels of EC number classification. The enzymatic steps are
  separated with a colon ':'. The first 3 levels or classification may be used to do a general descripction of a reaction.
  **These are the sequences used for the creation of the non redundant database (see bellow) and for ESS alignment**.
  
    - **ec4**: ESS using the 4 levels of classification. 
  
  - **arch**: Currently this column is empty. Originally, this column was used for the storing of a sequences of
  protein architectures. The architectures are defined according the protein family assignation in SuperFamily
  database.
  
  - **len**: Length of the ESS, measured as the number of enzymatic steps.
  
  - **map**: Name of the metabolic map from where the sequence was obtained.
  
  - **mapid**: KEGG metabolic map identifier from where the sequence was obtained. The code is created joining
  the species code (3 or 4 letters) and the map code (5 digit number). Example, for the Glycolisis/Gluconeogenesis
  of *E. coli* is eco00010.
  
  - **metabolism**: KEGG general type of metabolism from where the sequence was obtained.
  
  - **sp**: KEGG species identifier from where the sequence was obtained. It consist of a 3 to 4 letter code that
  abbreviates the species name.
  
  - **nrid**: Identifier of the non redundant database that represent the ESS.
  
  
### map_stat table.

  - **mapid**: KEGG metabolic map identifier from where the sequence was obtained. The code is created joining
  the species code (3 or 4 letters) and the map code (5 digit number). Example, for the Glycolisis/Gluconeogenesis
  of *E. coli* is eco00010. Unique identifier that represents a unique *kgml* file.
  
  - **name**: Full name of the metabolic map.
  
  - **sp**: KEGG species identifier of the metabolic map.
  
  - **reactions**: Number of reactions parsed from a metabolic map (*kgml* file).
  
  - **reversibles**: Number of reversible reactions in a metabolic map.
  
  - **genes**: Number of genes annotated in a metabolic map.
  
  - **ec3s**: Number of different EC numbers at 3 level of classification.
  
  - **ec4s**: Number of different EC numbers at 4 level of classification.
  
  - **no_EC**: Number of genes with no EC number assigned.
  
  - **nodes**: Number of nodes of the gene graph generated from the *kgml* file.
  
  - **starts**: Number of *start nodes* identified from in the metabolic map. Each *start node*
  is used for the creation of the *Breath First Search* tree.
  
  - **maplinks**: Number links to other maps present in each map.
  
  - **nseqs**: Number of ESSs generated from a map.
  
  - **metabolism**: KEGG general type of metabolism type category for a map.
  
### nrseqs tabel.

  - **nrid**: Unique identifier for the non redundant ESS table.
  
  - **ec3**: ESS using the first 3 levels of EC number classification. The enzymatic steps are
  separated with a colon ':'. **These are the sequences used for the ESS alignment**.
  
  - **len**: Length of the ESS, measured as the number of enzymatic steps.


## Creation of Enzymatic Step Sequences.




## Important.

- Only creates the ESS from those organisms that has a Enzyme list in ./Enzymes folder
  - Those organisms with enzyme list without folder in ./Maps will raise an error.
- Currently, the program must be executed completely. If a error is raised, the program
must be executed from the beginning.


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


