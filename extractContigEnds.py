#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Runar Furenes'
__email__  = 'runarfu@ifi.uio.no'

import argparse
import sys

FASTA_TEXTWIDTH = 80 # Use this width when breaking sequences into lines

def readMultiFASTA(filename):
    """
    Read a multi-FASTA file and return each sequence in a dictionary
    with FASTA headers as keys and (flattened) sequences as values.
    """
    f = open(filename)
    d = {}
    try:
        c = f.read()
        f.close()
        entries = c.split('>')
        for entry in entries[1:]:
            lines = entry.split('\n')
            header = lines[0]
            d[header] = ''.join(lines[1:])
    except:
        sys.stderr.write('Failed to read/parse FASTA-file %s\n' % filename)
    return d

def makePairedEnds(seqDict, nCut, makePairs=True):
    """
    Look at each sequence and make a left- and a right-end of it.
    If the length of a contig is < nCut, two parts is still created,
    but they are identical.
    """
    pairedContigs = {}
    for fastaKey in seqDict.keys():
        l = len(seqDict[fastaKey])
        if not makePairs or l < 2*nCut:
            allName = 'ALL_%s' % fastaKey
            pairedContigs[allName] = seqDict[fastaKey]
            continue
        minCut = min(nCut, l) # In case length is less than nCut, use whole contig
        left  = seqDict[fastaKey][:minCut]
        right = seqDict[fastaKey][l-minCut:]
        leftName  = 'LFT_%s' % fastaKey
        rightName = 'RGT_%s' % fastaKey
        pairedContigs[leftName]  = left
        pairedContigs[rightName] = right
    return pairedContigs

def writeContigsToFile(a, pairedContigs):
    """
    Write contig-end pairs (in a dict) to file in multi-FASTA format.
    """
    f = open(a.outputFile, 'w')
    try:
        for key in pairedContigs.keys():
            seq = pairedContigs[key]
            seq = '\n'.join(seq[i:i+FASTA_TEXTWIDTH] \
                for i in xrange(0, len(seq), FASTA_TEXTWIDTH))
            f.write('>%s\n%s\n' % (key, seq)) 
    except:
        sys.stderr.write('Failed to write contigs to %s\n' % a.outputFile)

if __name__ == '__main__':
    # Parse command line arguments for input and output files
    p = argparse.ArgumentParser(description='Extract ends from each side of\
    the sequences in a FASTA file.')
    p.add_argument('-i', '--inputFile', dest='inputFile', required=True, help='Input filename')
    p.add_argument('-o', '--outputFile', dest='outputFile', required=True, help='Output filename')
    p.add_argument('-n', '--nCut', dest='nCut', type=int, 
        help='Number of bases to cut from each contig-end')

    a = p.parse_args()

    seqDict = readMultiFASTA(a.inputFile)
    if a.nCut == 0:
        pairedContigs = makePairedEnds(seqDict, a.nCut, makePairs=False)
    else:
        pairedContigs = makePairedEnds(seqDict, a.nCut)
        
    writeContigsToFile(a, pairedContigs)

