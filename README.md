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

