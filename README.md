# Synopsis

Creating a dataset for a "true" tree for SARS-CoV-2.
Where we know the actual tree and phylogeny and create the genome assemblies to go with it.
Tree-making pipelines can test their algorithm against the true tree
and benchmark how well they reconstructed it.

# True tree

The true tree is where we find a tree and say "this is the tree" and accept that as fact.
In other words, all ancestor nodes are no longer hypothetical.
Branch lengths are no longer estimated.
All properties of the tree are _fact_.

Next, find an anchor genome. I have used Wuhan-1 as the anchor genome here.

After the anchor genome and true tree are determined, other parameters are determined
such as the mutation rates for different nucleotides.
More details can be found in [METHODS.md](METHODS.md).

Finally, the TreeToReads pipeline is run and all leaves of the tree are simulated
into genome assembly files.

# Why?

If you want to benchmark a phylogenetic method, then you need to know what the ideal tree is.
Since phylogenetic methods _infer_ a phylogeny, there is no real way to know the ideal or ground truth.
Therefore to understand what the ideal tree would be from a phylogenetic method, you would want to compare
against a true tree.
For SARS-CoV-2, this repo offers a true tree.

# Comparisons

Some suggestions for comparing your tree from your phylogenetic or clustering program are below

## Robinson-Foulds

This implementation is in my script under [lskScripts](https://github.com/lskatz/lskScripts/tree/master/scripts).
There are many other implementations of Robinson-Foulds (RF).
In my implementation, I also make 100 random trees
and compare the observed RF between the tree and true tree
vs the set of RFs between the random trees and true tree.

    treedist_wrapper.pl --method rf simtree.tre mytree.dnd --numtrees 100 --numcpus 10 \
      > rf.tsv 2> rf.og &

This was my output for Mashtree v1.
The interpretation is that there is an observable distance between the true tree and the Mashtree tree (obs=6740)
but that it is much closer to the true tree than chance alone (avg=7542, Z=-137, p < 1e101).

    Ref          Query         num  obs        avg        stdev   Z          p
    simtree.tre  mashtree.dnd  100  6740.0000  7542.2600  5.8579  -136.9524  0.00e+00

