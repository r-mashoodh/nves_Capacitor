GENRICH=/home/rm786/rds/rds-rmash-sHE3M4BjkMM/software/Genrich/Genrich
SORTDIR=/home/rm786/rds/rds-rmash-sHE3M4BjkMM/CutTag/sort_n
OUTDIR=/home/rm786/rds/rds-rmash-sHE3M4BjkMM/CutTag/peaks

mkdir -p $OUTDIR

# H3K4me3 - no control
$GENRICH \
    -t "$SORTDIR/CT1_H3K4_FC.sorted.n.bam,$SORTDIR/CT1_H3K4_NC.sorted.n.bam,$SORTDIR/CT2_H3K4_FC.sorted.n.bam,$SORTDIR/CT2_H3K4_NC.sorted.n.bam,$SORTDIR/CT3_H3K4_FC.sorted.n.bam,$SORTDIR/CT3_H3K4_NC.sorted.n.bam,$SORTDIR/CT4_H3K4_FC.sorted.n.bam,$SORTDIR/CT4_H3K4_NC.sorted.n.bam" \
    -o $OUTDIR/H3K4.noControl.narrowPeak \
    -f $OUTDIR/H3K4.noControl.log \
    -j -r -q 0.05 -y -v

## H3K4me3 - with IgG control (CT1 and CT2 only)
#$GENRICH \
#    -t "$SORTDIR/CT1_H3K4_FC.sorted.n.bam,$SORTDIR/CT1_H3K4_NC.sorted.n.bam,$SORTDIR/CT2_H3K4_FC.sorted.n.bam,$SORTDIR/CT2_H3K4_NC.sorted.n.bam,$SORTDIR/CT3_H3K4_FC.sorted.n.bam,$SORTDIR/CT3_H3K4_NC.sorted.n.bam,$SORTDIR/CT4_H3K4_FC.sorted.n.bam,$SORTDIR/CT4_H3K4_NC.sorted.n.bam" \
#    -c "$SORTDIR/CT1_IgG_FC.sorted.n.bam,$SORTDIR/CT1_IgG_NC.sorted.n.bam,$SORTDIR/CT2_IgG_FC.sorted.n.bam,$SORTDIR/CT2_IgG_NC.sorted.n.bam" \
#    -o $OUTDIR/H3K4.withControl.narrowPeak \
#    -f $OUTDIR/H3K4.withControl.log \
#    -j -r -q 0.05 -y -v

# H3K9 - no control
$GENRICH \
    -t "$SORTDIR/CT1_H3K9_FC.sorted.n.bam,$SORTDIR/CT1_H3K9_NC.sorted.n.bam,$SORTDIR/CT2_H3K9_FC.sorted.n.bam,$SORTDIR/CT2_H3K9_NC.sorted.n.bam,$SORTDIR/CT3_H3K9_FC.sorted.n.bam,$SORTDIR/CT3_H3K9_NC.sorted.n.bam,$SORTDIR/CT4_H3K9_FC.sorted.n.bam,$SORTDIR/CT4_H3K9_NC.sorted.n.bam" \
    -o $OUTDIR/H3K9.noControl.narrowPeak \
    -f $OUTDIR/H3K9.noControl.log \
    -j -r -q 0.05 -y -v

# H3K9 - with IgG control (CT1 and CT2 only)
#$GENRICH \
#    -t "$SORTDIR/CT1_H3K9_FC.sorted.n.bam,$SORTDIR/CT1_H3K9_NC.sorted.n.bam,$SORTDIR/CT2_H3K9_FC.sorted.n.bam,$SORTDIR/CT2_H3K9_NC.sorted.n.bam,$SORTDIR/CT3_H3K9_FC.sorted.n.bam,$SORTDIR/CT3_H3K9_NC.sorted.n.bam,$SORTDIR/CT4_H3K9_FC.sorted.n.bam,$SORTDIR/CT4_H3K9_NC.sorted.n.bam" \
#    -c "$SORTDIR/CT1_IgG_FC.sorted.n.bam,$SORTDIR/CT1_IgG_NC.sorted.n.bam,$SORTDIR/CT2_IgG_FC.sorted.n.bam,$SORTDIR/CT2_IgG_NC.sorted.n.bam" \
#    -o $OUTDIR/H3K9.withControl.narrowPeak \
#    -f $OUTDIR/H3K9.withControl.log \
#    -j -r -q 0.05 -y -v
