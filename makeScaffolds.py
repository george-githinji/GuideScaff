#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Runar Furenes'
__email__ = 'runarfu@ifi.uio.no'

import argparse
import sys
from collections import namedtuple
from itertools import izip, tee

FASTA_TEXTWIDTH = 80
ContigLinksEntry = namedtuple('ContigLinksEntry', 'GAP,ORIENTATION,CONTIGID')


def pairwise(it):
    """
    Creates an iterator with pairwise elements from it.
    """
    a, b = tee(it)
    next(b, None)
    return izip(a, b)


def reverseComplement(seq):
    """
    Create the reverse complement of a nucleotide sequence,
    i.e. replacing each nucleotide (all legal FASTA-symbols)
    with its complementary symbol, and reversing this sequence.
    """
    trans = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
        'M': 'K',
        'R': 'Y',
        'W': 'W',
        'S': 'S',
        'Y': 'R',
        'K': 'M',
        'V': 'B',
        'H': 'D',
        'D': 'H',
        'B': 'V',
        'X': 'X',
        'N': 'N'
    }
    #                            Replace   Reverse
    return ''.join(map(lambda x: trans[x], seq[::-1]))


def readMultiFASTA(filename):
    """
    Read a multi-FASTA file and return each sequence in a dictionary
    with FASTA headers as keys and (flattened) sequences as values.
    """
    d = {}
    try:
        with open(filename) as f:
            c = f.read()
    except IOError:
        sys.stderr.write('Failed to read/parse FASTA-file %s\n' % filename)
    entries = c.split('>')
    for entry in entries[1:]:
        lines = entry.split('\n')
        header = lines[0].split()[0]
        d[header] = ''.join(lines[1:])
    return d


def parseFile(filename, tupleType):
    """
    Parse file with contig-links with columns as specified in tupleType.
    Return dict with clusters as keys and a list of entries for each key.
    """
    try:
        with open(filename) as f:
            content = f.read()
    except IOError:
        sys.stderr.write('Failed to read %s\n' % filename)
        sys.exit(1)
    clusters = content.split('\n>')
    d = {}
    for c in clusters:
        lines = c.split('\n')
        header = lines[0].replace('>', '')
        entries = []
        for line in lines[1:]:
            elements = line.split('\t')
            # If something is odd about the #columns
            if len(elements) != len(tupleType._fields):
                continue
            entry = tupleType(*elements)
            entries.append(entry)
        d[header] = entries
    return d


def findUnusedContigs(d, contigs):
    """
    Figure out which contigs are not present in contig-links and
    return a dictionary with these.
    """
    allContigs = set(contigs.keys())
    used = set()
    for cluster in d.values():
        for entry in cluster:
            used.add(entry.CONTIGID)
    unusedContigs = allContigs.difference(used)
    d = {}
    for contig in unusedContigs:
        d[contig] = contigs[contig]
    return d


def merge(s1, s2, o):
    """
    Merge two string (sequences) together according to an overlap
    measure already calculated.
    """
    return s1 + s2[o:]


def nOverlap(s1, s2):
    """
    Calculate number of overlapping symbols between end of s1 and start of s2.
    """
    x = min(len(s1), len(s2))
    while x > 0:
        if s1[-x:] == s2[:x]:
            break
        x -= 1
    return x


def makeScaffolds(d, contigs):
    """
    Build actual scaffolds from the contig-links in dictionary d,
    with contig-sequences stored in dictionary contigs.
    """
    scaffolds = {} # Resulting scaffolds-dictionary
    for cluster in d.keys():
        scaffold = []
        i = 0
        while i < len(d[cluster]):
            c1 = d[cluster][i]
            gap = int(c1.GAP)
            if c1.ORIENTATION == '-':
                seq1 = reverseComplement(contigs[c1.CONTIGID])
            else:
                seq1 = contigs[c1.CONTIGID]
            # If gap-estimate is positive, append a corresponding number of
            # N's to it.
            if gap >= 0:
                scaffold.append(seq1)
                ns = 'N' * gap
                scaffold.append(ns)
                i += 1  # Go to next contig in contig-links
            # If gap-estimate is negative, try to merge the current contig
            # with the next.
            else:
                c2 = d[cluster][i + 1]
                if c2.ORIENTATION == '-':
                    seq2 = reverseComplement(contigs[c2.CONTIGID])
                else:
                    seq2 = contigs[c2.CONTIGID]
                    # Figure out if there is an ACTUAL overlap between the two
                # sequences.
                overlap = nOverlap(seq1, seq2)
                # If there is, merge them and append the merged contigs to
                # the scaffold.
                if overlap > 0:
                    merged = merge(seq1, seq2, overlap)
                    scaffold.append(merged)
                    i += 2  # Now both contigs are processed, skip to next one
                else:
                    # In this case, there were no actual overlap.
                    # Append sequence to scaffold and go to next contig in
                    # contig-links.
                    scaffold.append(seq2)
                    i += 1
        # Merge the entire scaffold created and add to scaffolds dictionary
        scaffolds[cluster] = ''.join(scaffold)
    return scaffolds


def writeScaffoldsToFile(a, scaffolds):
    """
    Write scaffolds to file in a multi-FASTA format.
    """
    try:
        f = open(a.outputFile, 'w')
        for cluster in sorted(scaffolds.keys()):
            header = cluster
            # Split sequence into lines according to pre-defined
            # FASTA width.
            seq = '\n'.join(scaffolds[cluster][i:i + FASTA_TEXTWIDTH] \
                            for i in xrange(0, len(scaffolds[cluster]), FASTA_TEXTWIDTH))
            f.write('>%s\n%s\n' % (header, seq))
    except IOError:
        sys.stderr.write('Failed to write scaffolds to %s\n' % a.outputFile)


if __name__ == '__main__':
    # Parse command line arguments for input and output files.
    p = argparse.ArgumentParser(description='Create scaffolds from contig-links\
    file')

    p.add_argument('-i', '--inputFile', dest='inputFile',
                   help='Input file with contig-links, orientations and gap estimates', required=True)
    p.add_argument('-o', '--outputFile', dest='outputFile',
                   help='Output file for scaffolds and unused contigs', required=True)
    p.add_argument('-c', '--contigsFile', dest='contigsFile',
                   help='File in FASTA-format containing all the contigs used in contig-links\
    file', required=True)
    a = p.parse_args()

    # Create dictionaries of contig-links and contig sequences
    contigLinksDict = parseFile(a.inputFile, ContigLinksEntry)
    contigsDict = readMultiFASTA(a.contigsFile)
    # Figure out which contigs are not used in any of the contig-links
    unusedContigs = findUnusedContigs(contigLinksDict, contigsDict)
    # Create scaffolds from the contig-links
    scf = makeScaffolds(contigLinksDict, contigsDict)
    # Merge created scaffolds with the unused contigs to a total dictionary
    scaffoldsAndContigs = dict(scf, **unusedContigs)
    # Write all these to a multi-FASTA file
    writeScaffoldsToFile(a, scaffoldsAndContigs)

