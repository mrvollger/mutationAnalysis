



parallel --jobs 12 'freebayes -f s288c.fasta {} > {/.}.vcf' ::: ../sam/*bam






