



# Conda environment.

```bash
conda create -n essdbpy2 -c conda-forge python=2 sqlite numpy matplotlib  networkx
```

# Execution.

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

# Important.

- Only creates the ESS from those organisms that has a Enzyme list in ./Enzymes folder
  - Those organisms with enzyme list without forlder in ./Maps will raise an error.
  
  




