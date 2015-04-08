#bowtie2-build s288c.fasta index/W303Index > index.results


parallel --jobs 12 'bowtie2 -t -q --very-sensitive -x index/W303Index -U {} -S ../sam/{.}.sam  &> freebayesAlign.results' ::: *fastq






