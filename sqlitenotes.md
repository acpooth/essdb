

# Select from one table into other
select nrseqs.nrid, nrseqs.ec3 from nrseqs INNER  JOIN seqs ON seqs.nrid =  nrseqs.nrid WHERE seqs.sp = 'hsa';
