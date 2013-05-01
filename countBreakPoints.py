#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
from itertools import tee, izip, ifilter
from collections import namedtuple


def pairwise(it):
    """
    Creates an iterator with pairwise elements from it.
    """
    a, b = tee(it)
    next(b, None)
    return izip(a, b)


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


def getOverlap(a, b, c, d):
    """
    Calculate overlap in ranges [a,b] and [c,d].
    """
    r = 0 if a == c and b == d else min(b, d) - max(a, c)
    if r >= 0:
        return r


def getDistance(a, b, c, d):
    """
    Calculate distance between two ranges [a,b] and [c,d],
    independent on their order.
    If there is overlap, return the overlap negated as distance.
    """
    overlap = getOverlap(a, b, c, d)
    if overlap:
        return -overlap
    return min([abs(max(a, b) - min(c, d)),
                abs(max(c, d) - min(a, b))])


def makeContigInfo(correctDict):
    """
    Create a new dictionary with original entry and a new field
    for which chromosome a contig resides in.
    """
    d = {}
    for chromosome, entries in correctDict.items():
        for entry in entries:
            contig = entry.CONTIGID
            d[contig] = ContigInfo(chromosome, entry)
    return d


def differentChromosomes(contigInfo, e1, e2):
    """
    Check if two entries maps to different chromosomes.
    """
    c1, c2 = e1.CONTIGID, e2.CONTIGID
    try:
        chr1, chr2 = contigInfo[c1].CHROMOSOME, contigInfo[c2].CHROMOSOME
        return chr1 != chr2
    except KeyError:
        return False


def differentRelativeOrientations(contigInfo, cluster):
    """
    Check every entry in a cluster and its orientation.
    Count as "alike" when orientation is equal to the correct
    orientation, and "unlike" when orientation differs.
    Report minimum value of these as total error.
    """
    alike = 0
    unlike = 0
    for entry in cluster:
        try:
            ci = contigInfo[entry.CONTIGID]
        except KeyError:
            sys.stderr.write("%s is not in genome.tiling\n" % entry.CONTIGID)
            continue
        if ci.TILINGENTRY.ORIENTATION != entry.ORIENTATION:
            unlike += 1
        else:
            alike += 1
    return min(unlike, alike)


def differentRelativeOrders(correctLists, cluster):
    """
    For each pair of contigs in a cluster, count how many of these
    maps to positions in the target genome which differs from the
    order implicitly given by the pair.
    """
    tot = 0
    for e1, e2 in pairwise(cluster):
        i1, i2 = None, None
        for correctList in correctLists:
            if e1.CONTIGID in correctList and e2.CONTIGID in correctList:
                i1 = correctList.index(e1.CONTIGID)
                i2 = correctList.index(e2.CONTIGID)
        if i1 is not None and i2 is not None:
            if i2 < i1:
                tot += 1
    return tot


def gapEstimatesExceedsDelta(contigInfo, e1, e2, delta=500):
    """
    Compare gap-estimate with true distance between two contigs.
    """
    gapEstimate = int(e1.GAP)
    correctGap = trueDistance(contigInfo, e1, e2)
    if abs(gapEstimate - correctGap) > delta:
        return True
    return False


def trueDistance(contigInfo, e1, e2):
    """
    Calculate distance between two contigs in the target genome,
    if they are present in the mapping.
    """
    try:
        c1, c2 = e1.CONTIGID, e2.CONTIGID
        start1 = int(contigInfo[c1].TILINGENTRY.START)
        end1 = int(contigInfo[c1].TILINGENTRY.END)
        start2 = int(contigInfo[c2].TILINGENTRY.START)
        end2 = int(contigInfo[c2].TILINGENTRY.END)
        dist = getDistance(start1, end1, start2, end2)
    except KeyError:
        return 0
    return dist


def countBreakPoints(correctDict, suggested, contigInfo, deltaValues):
    """
    Count breakpoints for each pair of contigs in suggested contig-links list.
    Print number of breakpoints and number of different errors to stdout.
    """
    dc = 0  # Different chromosomes
    dori = 0  # Different orientation
    dord = 0  # Different order
    gaps = [0 for _ in range(len(deltaValues))]  # Gap estimate exceeded delta value
    nPairs = 0.0
    nContigs = 0.0

    # Make a correct-list for easy lookup
    correctLists = []
    for cluster in correctDict.keys():
        correctLists.append(map(lambda entry: entry.CONTIGID, correctDict[cluster]))

    # Go through each cluster in the contig-links file
    # and count breakpoints.
    for cluster in suggested.keys():
        nContigs += len(suggested[cluster])
        for e1, e2 in pairwise(suggested[cluster]):
            nPairs += 1
            if differentChromosomes(contigInfo, e1, e2):
                dc += 1
            for i in range(len(deltaValues)):
                if gapEstimatesExceedsDelta(contigInfo, e1, e2,
                                            delta=deltaValues[i]):
                    gaps[i] += 1

        # Count number of wrong relative orientations in the entire cluster
        dori += differentRelativeOrientations(contigInfo, suggested[cluster])

        # Count number of relative wrong orders in the cluster both ways,
        # use least of these two numbers.
        dord1 = differentRelativeOrders(correctLists, suggested[cluster])
        dord2 = differentRelativeOrders(correctLists, reversed(suggested[cluster]))
        dord += min(dord1, dord2)

    print('N_PAIRS\t%d' % nPairs)
    print('DIFF_CHROMOSOMES\t%d' % dc)
    print('DIFF_ORIENTATION\t%d' % dori)
    print('DIFF_ORDER\t%d' % dord)
    for i in range(len(deltaValues)):
        print('GAP_ERROR_%d\t%d' % (deltaValues[i], gaps[i]))

    if nContigs > 0:
        rel_dc = dc / nContigs
        rel_dori = dori / nContigs
    else:
        rel_dc = 0
        rel_dori = 0
    if nPairs > 0:
        rel_ord = dord = dord / nPairs
    else:
        rel_ord = 0

    print('REL_DIFF_CHROMOSOMES\t%f' % rel_dc)
    print('REL_DIFF_ORIENTATION\t%f' % rel_dori)
    print('REL_DIFF_ORDER\t%f' % rel_ord)
    for i, deltaValue in enumerate(deltaValues):
        if nPairs > 0:
            rel_gap = gaps[i] / nPairs
        else:
            rel_gap = 0
        print('REL_GAP_ERROR_%d\t%f' % (deltaValue, rel_gap))


if __name__ == '__main__':
    # Parse command line arguments for input and output files
    p = argparse.ArgumentParser(description='Validate contig ordering')
    p.add_argument('--suggestedContigOrderFile', '-s', dest='suggestedFile',
                   help='Suggested contig order file', required=True)
    p.add_argument('--correctContigOrderFile', '-c', dest='correctFile',
                   help='Correct contig order file', required=True)
    a = p.parse_args()

    TilingEntry = namedtuple('TilingEntry',
                             'START,END,GAP,LENGTH,COV,AVGID,ORIENTATION,CONTIGID')
    ContigLinksEntry = namedtuple('ContigLinksEntry',
                                  'GAP,ORIENTATION,CONTIGID')
    ContigInfo = namedtuple('ContigInfo',
                            'CHROMOSOME,TILINGENTRY')

    correctDict = parseFile(a.correctFile, TilingEntry)
    suggested = parseFile(a.suggestedFile, ContigLinksEntry)
    contigInfo = makeContigInfo(correctDict)

    # Delta-values determines the resolution used when counting incorrect
    # gap-estimates.
    deltaValues = [100, 500, 1000, 10000]

    # Count breakpoints.
    countBreakPoints(correctDict, suggested, contigInfo, deltaValues)