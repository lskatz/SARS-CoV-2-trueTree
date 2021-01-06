# Tips

Random thoughts go here 

## Recreate the alignment

If you want to have a fasta alignment to play around with,
the basic workflow is to split the vcf into samples
and then make a fasta file from each alignment,
and finally a multifasta file.

    # index the vcf file
    bcftools index sim.vcf.gz

    mkdir -v makeAlignment
    cd makeAlignment

    # 1. query -l to find all sample names
    # 2. send to xargs to multithread (-P)
    # 3. bcftools view to get the individual vcf files
    # 4. index the individual vcf file
    bcftools query -l ../sim.vcf.gz | \
      xargs -n 1 -P 12 bash -c '
        echo -n . >&2;
        bcftools view -l 9 -s $0 -Oz -o $0.vcf.gz ../sim.vcf.gz
        bcftools index $0.vcf.gz
      '
    
    # bcftools consensus step gives you a consensus from a vcf file
    for file in *.vcf.gz; do 
      name=$(basename $file .vcf.gz);
      bcftools consensus --fasta-ref ../../data/nextstrain-2020-01-04/wuhan-1.fasta --sample $name $file | \
        sed "s/>.*/>$name/" > $name.fasta;
    done;
    
    # cat the fasta files to make the alignment
    cat *.fasta > aln.fas

    # Save your work
    mv aln.fas ../sim.aln.fas

