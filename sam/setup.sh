

parallel --jobs 12 'samtools view -bS {} -o {.}.bam' ::: *sam
parallel --jobs 12 'samtools sort {} {.}' ::: *bam
parallel --jobs 12 'samtools index {}' ::: *bam

