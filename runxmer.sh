#!/bin/bash

GENOME=$1  # Filename of FASTA-formatted guiding genome
QUERY=$2   # Filename of FASTA-formatted contig-pairs
NCUT=$3
LIMIT="90" # Threshold determining whether to use nucmer or promer
PREFIX=$GENOME.$NCUT
PROMERFLAGS="-o -p $PREFIX"

# Test if mappings are already present
if test -n "$(shopt -s nullglob; echo $PREFIX.*)"
then
    :
else
    # Run nucmer first to calculate average identity value
    nucmer -o -p $PREFIX $GENOME $QUERY

    # Calculate the average identity score for all contig mappings
    AVGID=$(cat $PREFIX.coords | awk -F " " 'NF==13{sum+=$10; count++} END {print sum/count}')
    echo "Average identity value is $AVGID" >&2
    belowThreshold=$(echo "$AVGID < $LIMIT" | bc) # Use bc to check this

    # If average identity score is below the preset threshold, run promer
    if [ $belowThreshold -eq "1" ]
    then
        echo "Using promer" >&2
        promer $PROMERFLAGS $GENOME $QUERY
    else
        echo "Using results from nucmer" >&2
    fi

fi

show-tiling $PREFIX.delta > $PREFIX.tiling

