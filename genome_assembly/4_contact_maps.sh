## GENERATING HiC contact maps
BIN='/mnt/home5/zoology/rm786/beetle_HiC/noHiC/yahs_output/yahs.out.bin'
AGP='/mnt/home5/zoology/rm786/beetle_HiC/noHiC/yahs_output/yahs.out_scaffolds_final.agp'
FAI='/mnt/home5/zoology/rm786/beetle_HiC/noHiC/Ni2_pacbio.asm.fa.fai'

(~/software/yahs/juicer pre $BAM $AGP $FAI | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)

(java -jar -Xmx32G ~/software/juicer_tools_1.22.01.jar pre alignments_sorted.txt out.hic.part ../scaffolds_final.chrom.sizes) && (mv out.hic.part out.hic)

## files for JBAT (JUICEBOX MANUAL CURATION)
~/software/yahs/juicer pre -a -o out_JBAT $BIN $AGP $FAI > out_JBAT.log 2>&1

(java -jar -Xmx32G ~/software/juicer_tools_1.22.01.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)


## manually fix errors and regenerate assembly fasta

PREFIX='Nic2_outJBAT'
REVIEW='out_JBAT.rev.review.assembly'
AGP='/mnt/home5/zoology/rm786/beetle_HiC/noHiC/yahs_output/yahs.out_scaffolds_final.agp'
## PRE SCAFFOLD
CONTIGS='/mnt/home5/zoology/rm786/beetle_HiC/noHiC/Ni2_pacbio.asm.fa'

~/software/yahs/juicer post -o $PREFIX $REVIEW $AGP $CONTIGS
