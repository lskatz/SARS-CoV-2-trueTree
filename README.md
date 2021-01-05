# Synopsis

Creating a dataset for a "true" tree for SARS-CoV-2.
Where we know the actual tree and phylogeny and create the genome assemblies to go with it.
Tree-making pipelines can test their algorithm against the true tree
and benchmark how well they reconstructed it.

# Methods

* TreeToReads
  * Dependencies
  * TreeToReads
* download
  * GISAID MSA and tree (deprecated)
  * Nextstrain
    * Go to nCov page, go to the bottom, click on download for tree and metadata.  Can't do this on CLI.
    * mv the tree and metadata to `data/nextstrain-2020-01-04/`
    * Get the wuhan-1 genome (first line of the metadata field) for the anchor genome in TreeToReads
    
    ```
        wget 'https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=LR757998.1&rettype=fasta' -O data/nextstrain-2020-01-04/wuhan-1.fasta
        sed -i.bak 's|>.*|>Wuhan/WH01/2019|' wuhan-1.fasta # so the fasta matches the tree
        # Download the other sequences to a file.
        # Not all ncbi accessions are available.
        # Need to get NCBI_API_KEY in the env if possible.
        #  => will get an error of 429: too many requests more often if you don't have this in env.
        cut -f 9 nextstrain_ncov_global_metadata.tsv | tail -n +2 | grep . | xargs -n 1 -P 1 bash -c 'numtries=0; while [ $numtries -lt 3 ]; do wget "https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$0&rettype=fasta&ncbi_api_key=$NCBI_API_KEY" -O - && break; sleep 1; numtries=$(($numtries+1)); done' > in.fasta
        # Make an alignment
        mafft in.fasta > mafft.fasta
        # Prepare the alignment for the next steps
        cat data/nextstrain-2020-01-04/mafft.fasta | perl -plane 'if(/>/){$i++; $_=">$i";}' | bgzip -c > anonDeflines.fasta.gz
        samtools faidx anonDeflines.fasta.gz
    ```

    * Fix the tree so that it is fit for `seq-gen`.

          perl -MBio::TreeIO -e '
            $tree=Bio::TreeIO->new(-file=>"nextstrain_ncov_global_tree.resolved.nwk")->next_tree; 
            $tree->force_binary; 
            $tree->contract_linear_paths; 
            for my $node($tree->get_nodes){
              $node->branch_length || $node->branch_length(rand(1e-7)); 
              if(!$node->is_Leaf){
                $node->id("");
              } 
            } 
            print $tree->as_text("newick")."\n";
          ' > anonymized.nwk

* prepare it for random sampling
  * rename seqids to integers
  * compress with bgzip
  * index with `samtools faidx`
  
        samtools faidx anonDeflines.fasta.gz

* Evolutionary model - some evolutionary parameters are required for TreeToReads.
  * rate matrix for 100 reps, 100 random seqs each. Set the integer on `shuf` to the number of sequences in the input alignment file.

```
        (set -e; for rep in `seq 1 100`; do echo "rep $rep" >&2; rm -f RAxML_*.T2 >&2; shuf -i 1-291834 -n 100 | xargs -n 1 samtools faidx anonDeflines.fasta.gz > 100.fasta; raxmlHPC -p $RANDOM -m GTRGAMMA -n T2 -s 100.fasta >& /dev/null; grep alpha RAxML_info.T2; done) > raxml.gtrgamma.rates.txt
        cat raxml.gtrgamma.rates.txt | perl -lane 'if(/ct gt: (.+)/){$d=$1; $d=~s/\s+/\t/g; print $d;}' | datamash mean 1 sstdev 1 mean 2 sstdev 2 mean 3 sstdev 3 mean 4 sstdev 4 mean 5 sstdev 5 mean 6 sstdev 6 | perl -lane '$l=""; for(my $i=0;$i<@F;$i+=2){$l.=sprintf("%0.2f+-%0.2f ", $F[$i], $F[$i+1]);} print $l;'
        #gisaid     => 0.23+-0.05 0.84+-0.11 0.14+-0.05 0.21+-0.05 2.80+-0.27 1.00+-0.00
        #nextstrain => 0.25+-0.04 0.82+-0.10 0.15+-0.03 0.27+-0.05 2.99+-0.29 1.00+-0.00
```

  * Find the number of variable sites with more random sampling

        for rep in `seq 1 100`; do echo "rep $rep" >&2; shuf -i 1-316 -n 100 | xargs bash -c 'samtools faidx anonDeflines.fasta.gz $@' | goalign stats | grep variable; done | datamash mean 2 sstdev 2
        #gisaid     => 7523.82        1308.1794534483
        #nextstrain => 6399.42        1519.2121441688

* Run TreeToReads. Use parameters as found above.  I decided to exclude ART in my path to avoid read simulation.

    #REQUIRED PARAMETERS
    # tree must be newick or Nexus format, and include branch lengths
    treefile_path = data/nextstrain-2020-01-04/anonymized.nwk
    number_of_variable_sites = 6400 
    # base_genome_name should be the label of a tip in your tree
    base_genome_name = Wuhan/WH01/2019
    base_genome_path = data/nextstrain-2020-01-04/wuhan-1.fasta
    output_dir = TTR.out

    #parameters of evolutionary model (comma seperated), in order ac, ag, at, cg, ct, gc (gc = 1)
    rate_matrix = 0.25,0.82,0.15,0.27,2.99,1.00

    #parameters for read simulation
    coverage = 100 #either an integer or a file name of a comma delimited file with tip names and coverage

    #OPTIONAL PARAMETERS
    prefix = sim_ #optional prefix prepended to sequence names, default is using orginal sequence names

    #Optional evolutionary model parameters
    gamma_shape = 5 #dafault is no rate variation across sites

    #parameters for clustering of variable site locations (OPTIONAL)
    mutation_clustering = ON
    percent_clustered = 0.25 #The percentage of variable sites whose distance to another site is drawn from the clustering distribution
    exponential_mean = 125 #Minimum allowed value = 2

    #ART Optional parameters (for more fine grained control ART can be run seperately on the mutated genomes found in outdir/fasta_files)
    error_model1 = example/ErrprofR1.txt  # If you haven't generated have one of your own using ART, you can use one supplied by ART.
    error_model2 = example/ErrprofR2.txt  # Un-comment these lines (delete the first #) to set a non-default error profile
    read_length = 150 #maximum value with example error profile is 150, use a default or generate adifferent error profile for longer reads.
    fragment_size = 380
    stdev_frag_size = 120

