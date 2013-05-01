#!/bin/bash
# Runar Furenes
# 2013-05-01

# Input variables
DIR=$1
NGUIDES=$2
WINDOWSIZE=$3
THRESHOLD=$4
NCUT=$5

GUIDE_DIR=$DIR/guidingGenomes

command_exists () {
    type "$1" &> /dev/null ;
}

################################################################################
##                           MAIN PIPELINE                                    ##
################################################################################

# Make paired ends of contigs
extractContigEnds.py --inputFile  $DIR/contigs.fasta \
                     --outputFile $DIR/contigEnds.fasta \
					 --nCut $NCUT
echo "> Made contig pairs with nCut=$NCUT" >&2

# Map the contig-pairs to each of the guiding genomes in parallel
# using GNU Parallel if available
if command_exists parallel
then
	parallel "runxmer.sh {} $DIR/contigEnds.fasta $NCUT" ::: $GUIDE_DIR/*.fasta
else
    for f in $(ls $GUIDE_DIR/*.fasta)
    do
        runxmer.sh $f $DIR/contigEnds.fasta $NCUT
    done
fi
echo "> Aligned contig ends to all guiding genomes" >&2

# Make contig links from tiling files produced by nucmer or promer
# using only the ones produced with NCUT ends, if several is present
makeContigLinks.py --inputFiles $GUIDE_DIR/*.fasta.$NCUT.tiling \
                   --output     $DIR/contigLinks \
				   --nGuides    $NGUIDES \
				   --windowSize $WINDOWSIZE \
				   --threshold  $THRESHOLD
echo "> Created contig-links from tiling-files" >&2

# Use contig links to build scaffolds in FASTA format
makeScaffolds.py --inputFile   $DIR/contigLinks \
                 --outputFile  $DIR/scaffolds.fasta \
				 --contigsFile $DIR/contigs.fasta
echo "> Created scaffolds in scaffolds.fasta" >&2

################################################################################
##                           EVALUATION                                       ##
################################################################################

# If not genome.tiling is already present,
# align contigs to true genome
if [ -e "$DIR/genome.tiling" ]
then
	:
else
	# Map the original contigs to the target genome for later evaluation
	nucmer -o -p $DIR/genome $DIR/genome.fasta $DIR/contigs.fasta
	show-tiling -c $DIR/genome.delta > $DIR/genome.tiling
	echo "> Mapped contigs to target genome and created tiling" >&2
fi

# Calculate values
TARGET_GENOME=$(head -n1 $DIR/genome.fasta | awk '{print $2 " " $3}')
GENOME_SIZE=$(genomeSize.py $DIR/genome.fasta | awk -F"\t" '{sum+=$1} END {print sum}')
CONTIG_SET=$(basename $(readlink -m $DIR/contigs.fasta))
CONTIGS=$(grep -c '^>' $DIR/contigs.fasta)
MIN_CONTIG_LENGTH=$(genomeSize.py $DIR/contigs.fasta | sort -n | head -n1 | awk -F"\t" '{print $1}')
MAX_CONTIG_LENGTH=$(genomeSize.py $DIR/contigs.fasta | sort -n | tail -n1 | awk -F"\t" '{print $1}')
SCAFFOLDS=$(grep -c '^>' $DIR/contigLinks)
MIN_SCAFFOLD_LENGTH=$(genomeSize.py $DIR/scaffolds.fasta | grep 'Scaffold' | sort -n | head -n1 | awk -F"\t" '{print $1}')
MAX_SCAFFOLD_LENGTH=$(genomeSize.py $DIR/scaffolds.fasta | grep 'Scaffold' | sort -n | tail -n1 | awk -F"\t" '{print $1}')
SCAFFOLDS_SIZE=$(genomeSize.py $DIR/scaffolds.fasta | awk -F"\t" '{sum+=$1} END {print sum}')
COVERAGE=$(echo -e "scale=2;$SCAFFOLDS_SIZE / ($GENOME_SIZE * 1.0)" | bc -l)
perl assemblathon_stats.pl $DIR/contigs.fasta > $DIR/assemblathonStatsContigs.txt
perl assemblathon_stats.pl $DIR/scaffolds.fasta > $DIR/assemblathonStatsScaffolds.txt
N50_CONTIGS=$(grep "N50 contig length" $DIR/assemblathonStatsContigs.txt | awk '{print $4}')
N50_SCAFFOLDS=$(grep "N50 scaffold length" $DIR/assemblathonStatsScaffolds.txt | awk '{print $4}')
N_CONTIGS_USED=$(grep -c '[+-]' $DIR/contigLinks)
REL_CONTIGS_USED=$(echo -e "scale=3;$N_CONTIGS_USED / ($CONTIGS * 1.0)" | bc -l)
countBreakPoints.py --correctContigOrderFile $DIR/genome.tiling --suggestedContigOrderFile $DIR/contigLinks > $DIR/breakpoints

# Report
echo -e "NCUT\t$NCUT"
echo -e "NGUIDES\t$NGUIDES"
echo -e "WINDOWSIZE\t$WINDOWSIZE"
echo -e "THRESHOLD\t$THRESHOLD"
echo -e "TARGET_GENOME\t$TARGET_GENOME"
echo -e "GENOME_SIZE\t$GENOME_SIZE"
echo -e "CONTIG_SET\t$CONTIG_SET"
echo -e "NCONTIGS\t$CONTIGS"
echo -e "MIN_CONTIG_LENGTH\t$MIN_CONTIG_LENGTH"
echo -e "MAX_CONTIG_LENGTH\t$MAX_CONTIG_LENGTH"
echo -e "NSCAFFOLDS\t$SCAFFOLDS"
echo -e "MIN_SCAFFOLD_LENGTH\t$MIN_SCAFFOLD_LENGTH"
echo -e "MAX_SCAFFOLD_LENGTH\t$MAX_SCAFFOLD_LENGTH"
echo -e "COVERAGE\t$COVERAGE"
echo -e "N50_CONTIGS\t$N50_CONTIGS"
echo -e "N50_SCAFFOLDS\t$N50_SCAFFOLDS"
echo -e "N_CONTIGS_USED\t$N_CONTIGS_USED"
echo -e "REL_CONTIGS_USED\t$REL_CONTIGS_USED"
cat $DIR/breakpoints
echo -e "\n"

