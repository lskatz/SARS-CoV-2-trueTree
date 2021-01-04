# Synopsis

Creating a dataset for a "true" tree for SARS-CoV-2.
Where we know the actual tree and phylogeny and create the genome assemblies to go with it.
Tree-making pipelines can test their algorithm against the true tree
and benchmark how well they reconstructed it.

# Methods

* TreeToReads
  * Dependencies
  * TreeToReads
* download GISAID MSA
* prepare it for random sampling
  * rename seqids to integers
  * compress with bgzip
  * index with `samtools faidx`
  
        samtools faidx anonDeflines.fasta.gz

* Evolutionary model - some evolutionary parameters are required for TreeToReads.
  * rate matrix for 100 reps, 100 random seqs each

        (set -e; for rep in `seq 1 100`; do echo "rep $rep" >&2; rm -f RAxML_*.T2 >&2; shuf -i 1-291834 -n 100 | xargs -n 1 samtools faidx anonDeflines.fasta.gz > 100.fasta; raxmlHPC -p $RANDOM -m GTRGAMMA -n T2 -s 100.fasta >& /dev/null; grep alpha RAxML_info.T2; done) > raxml.gtrgamma.rates.txt
        cat raxml.gtrgamma.rates.txt | perl -lane 'if(/ct gt: (.+)/){$d=$1; $d=~s/\s+/\t/g; print $d;}' | datamash mean 1 sstdev 1 mean 2 sstdev 2 mean 3 sstdev 3 mean 4 sstdev 4 mean 5 sstdev 5 mean 6 sstdev 6 | perl -lane '$l=""; for(my $i=0;$i<@F;$i+=2){$l.=sprintf("%0.2f+-%0.2f ", $F[$i], $F[$i+1]);} print $l;'
        # => 0.23+-0.05 0.84+-0.11 0.14+-0.05 0.21+-0.05 2.80+-0.27 1.00+-0.00

  * Find the number of variable sites with more random sampling

        for rep in `seq 1 100`; do echo "rep $rep" >&2; shuf -i 1-291834 -n 100 | xargs -n 1 -P 1 samtools faidx anonDeflines.fasta.gz | goalign stats | grep variable; done | datamash mean 2 sstdev 2
        # => 7523.82        1308.1794534483
