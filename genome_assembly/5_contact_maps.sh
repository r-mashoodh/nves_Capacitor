## GENERATING HiC contact maps

## these are yahs output files
BIN='/mnt/home5/zoology/rm786/beetle_HiC/pilon/round2/map_hic/scaffold/Ni2.purged.polished.bin'
AGP='/mnt/home5/zoology/rm786/beetle_HiC/pilon/round2/map_hic/scaffold/Ni2.purged.polished_scaffolds_final.agp'

## this is meant to be the fai / sizes of the fasta used to scaffold (purged/polished)
#samtools faidx Ni2_purged.pilon2.fa.fasta
#cut -f1,2 Ni2_purged.pilon2.fa.fasta.fai > Ni2_purged.pilon2.fa.fasta.sizes
FAI='/mnt/home5/zoology/rm786/beetle_HiC/pilon/round2/Ni2_purged.pilon2.fa.fasta.fai'
SIZES='/mnt/home5/zoology/rm786/beetle_HiC/pilon/round2/Ni2_purged.pilon2.fa.fasta.sizes'

## getting files for Juicer
(~/software/yahs/juicer pre $BIN $AGP $FAI | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)

(java -jar -Xmx32G ~/software/juicer_tools_1.22.01.jar pre alignments_sorted.txt out.hic.part $SIZES) && (mv out.hic.part out.hic)

## files for JBAT (JUICEBOX MANUAL CURATION)
~/software/yahs/juicer pre -a -o out_JBAT $BIN $AGP $FAI > out_JBAT.log 2>&1

(java -jar -Xmx32G ~/software/juicer_tools_1.22.01.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)

# rename files to keep track of assembly version
for file in *out*; do mv "$file" "${file/out/Ni2.purged.filt.scaff}"; done

# get assembly stats
~/software/gfastats/build/bin/gfastats Ni2_purged.pilon2.fa.fasta --nstar-report --stats
