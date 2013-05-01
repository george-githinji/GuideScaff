#!/usr/bin/env python
# -*- coding: utf-8 -*-
__AUTHOR__ = 'Runar Furenes'
__email__  = 'runarfu@ifi.uio.no'

from collections import namedtuple
from itertools import *
from numpy import median, average
from operator import itemgetter, attrgetter
from os import path
import argparse
import glob
import sys

TilingEntry = namedtuple('TilingEntry',\
    'START,END,GAP,LENGTH,COV,AVGID,ORIENTATION,CONTIGID')
ContigEntry = namedtuple('ContigEntry',\
    'ORIENTATION,CONTIGID,END')
ScaffoldLine = namedtuple('ScaffoldLine',\
    'GAP,ORIENTATION,CONTIG')

def harmonicMean(l):
    fil = filter(lambda x: x!= 0, l)
    if len(fil) > 0:
        return len(fil) / sum(1.0 / x for x in fil)
    return average(l)

def oppositeContigEntry(contigEntry):
    """
    Pick the contig end at the opposite end, but in the same
    orientation.
    """
    if contigEntry.END == 'LFT':
        return ContigEntry(contigEntry.ORIENTATION, contigEntry.CONTIGID, 'RGT')
    if contigEntry.END == 'RGT':
        return ContigEntry(contigEntry.ORIENTATION, contigEntry.CONTIGID, 'LFT')
    return ContigEntry(contigEntry.ORIENTATION, contigEntry.CONTIGID, 'ALL')

def oppositeOrientation(contigEntry):
    """
    Pick the contig end at the opposite end, but in the same
    orientation.
    """
    if contigEntry.ORIENTATION == '+':
        return ContigEntry('-', contigEntry.CONTIGID, contigEntry.END)
    return ContigEntry('+', contigEntry.CONTIGID, contigEntry.END)

def pairwise(it):
    """
    Creates an iterator with pairwise elements from it.
    """
    a, b = tee(it)
    next(b, None)
    return izip(a, b)

def parseTilingFile(filename):
    """
    Parse tiling file (from nucmer or promer -> show-tiling).
    Return list of clusters for each chromosome in the tiling-file.
    """
    clusters = []
    try:
        f = open(filename)
        _ = f.readline()
        cluster = []
        for line in f:
            if line == '\n':
                break
            if line.startswith('>'):
                clusters.append(cluster)
                cluster = []
                continue
            entry = TilingEntry(*line.split())
            cluster.append(entry)
        if len(cluster) > 0:
            clusters.append(cluster)
        f.close()
    except:
        sys.stderr.write('Failed to read/parse file %s\n' % filename)
    return clusters

def overlap(a,b,c,d):
    """
    If ranges [a,b] and [c,d] overlap, report it negated.
    overlap(0,10,8,20) -> 2
    """
    r = 0 if a==c and b==d else min(b,d)-max(a,c)
    if r>=0: return r

def distance(a,b,c,d):
    """
    Calculate distance between two contigs A and B
    where A starts at index a and ends at index b and
          B starts at index c and ends at index c.
    Return overlap as negated number if an overlap exist.
    """
    o = overlap(a,b,c,d)
    if o:
        return -o
    return min([abs(max(a,b) - min(c,d)),\
                abs(max(c,d) - min(a,b))])

def buildMatrix(tilingsList, w):
    """
    Build a distance matrix from the tilings list, where each entry within
    a sliding window of size w is used to make entries.
    """
    M = {}
    for tilings in tilingsList:
        for chromosome in tilings:
            for i in xrange(len(chromosome)-w):
                e1 = chromosome[i]
                # Look at a window from e1 and ahead in the list
                for e2 in islice(chromosome, i+1, i+w):
                    start1, end1 = int(e1.START), int(e1.END)
                    start2, end2 = int(e2.START), int(e2.END)

                    # Remove labels indicating which end is used
                    c1 = e1.CONTIGID[4:]
                    c2 = e2.CONTIGID[4:]

                    # Extract labels indicating which end is used
                    paired1 = e1.CONTIGID[0:3]
                    paired2 = e2.CONTIGID[0:3]
                    ce1 = ContigEntry(e1.ORIENTATION, c1, paired1)
                    ce2 = ContigEntry(e2.ORIENTATION, c2, paired2)

                    # Make an entry for the opposite relative order
                    # and orientation
                    ce1opp = oppositeOrientation(ce1)
                    ce2opp = oppositeOrientation(ce2)

                    # Find distance between contigs
                    dist = distance(start1, end1, start2, end2)

                    # Append distance to matrix entries
                    M.setdefault(ce1, {})
                    M[ce1].setdefault(ce2, [])
                    M[ce1][ce2].append(dist)

                    M.setdefault(ce2opp, {})
                    M[ce2opp].setdefault(ce1opp, [])
                    M[ce2opp][ce1opp].append(dist)

                    #M.setdefault(ce2, {})
                    #M[ce2].setdefault(ce1, [])
                    #M[ce2][ce1].append(dist)

    return M

def makeConsensusMatrix(M, threshold, distancesFunc):
    """
    Process distance matrix M with multiple distances in lists.
    Apply distancesFunc on each of these lists to get an estimate
    of the distances.
    """
    MC = {}
    for m1 in M.keys():
        for m2 in M[m1].keys():
            # If less than threshold of the guiding genomes had distance
            # information about this (within the window when distance matrix
            # were created) - discard these values.
            if len(M[m1][m2]) < threshold:
                continue
            consensusDistance = distancesFunc(M[m1][m2])
            MC.setdefault(m1, {}) # If it doesn't already exist
            MC[m1][m2] = int(consensusDistance)
    return MC

def buildPaths(M):
    """
    Based on a consensus distance matrix M, build paths of contig ends
    starting with a random unused contig end not yet in one of the created
    paths.
    Return a list of paths.
    """
    paths = []
    seenContigs = set()
    unseenContigs = set(M.keys())
    # As long as there are more contigs not yet included in any of the paths
    while len(unseenContigs) > 0:
        startNode = unseenContigs.pop()
        if startNode.CONTIGID in seenContigs:
            continue
        seenContigs.add(startNode.CONTIGID)
        # Grow a new path starting with the chosen contig end
        path = expandPath(M, seenContigs, startNode)
        # Keep only paths with at least two entries
        if len(path) > 1:
            paths.append(path)
            # Update the sets keeping track of seen and unseen contig ends
            for entry in path:
                seenContigs.add(entry.CONTIGID)
                unseenContigs.discard(entry)
    return paths

def expandPath(M, dropList, initial):
    """
    Grow a path from initial as long as there is a corresponding entry
    in the distance matrix not yet processed (in the dropList).
    """
    path = [initial]
    while True:
        # Start by skipping from one contig end to the corresponding
        # opposite end (these belongs to the same contig, and should
        # therefore always follow each other in paths.
        oppositeEnd = oppositeContigEntry(path[-1])
        if not M.has_key(oppositeEnd) or oppositeEnd in dropList:
            break
        path.append(oppositeEnd)
        dropList.add(oppositeEnd.CONTIGID)
        # Try to find the closest match to the latest entry in the path
        n = bestNextMatch(M, dropList, path[-1])
        if n:
            path.append(n)
            continue
        break
    return path

def bestNextMatch(M, dropList, n):
    """
    Find the closest match for n in distance matrix M not in the dropList.
    """
    candidates = filter(lambda (entry,dist): entry.CONTIGID not in dropList,
                                              M[n].items())
    if len(candidates) > 0:
        return sorted(candidates, key=itemgetter(1))[0][0]
    return None

def makeScaffolds(M, paths):
    """
    Make scaffold-lines for each path in paths.
    """
    scaffolds = []
    for path in paths:
        scaffold = []
        for ce1, ce2 in pairwise(path):
            if ce1.CONTIGID == ce2.CONTIGID:
                continue
            orientation = ce1.ORIENTATION
            gap         = M[ce1][ce2]
            scaffoldLine = ScaffoldLine(gap, orientation, ce1.CONTIGID)
            scaffold.append(scaffoldLine)
        ceLast = path[-1]
        scaffoldLine = ScaffoldLine(0, ceLast.ORIENTATION, ceLast.CONTIGID)
        scaffold.append(scaffoldLine)
        scaffolds.append(scaffold)
    return scaffolds

def writeScaffoldLinesToFile(filename, scaffolds):
    """
    Print scaffolds to line in format similar to output of MUMmer's
    show-tiling.
    """
    f = open(filename, 'w')
    i = 1
    for scaffold in scaffolds:
        if len(scaffold) < 2:
            continue
        f.write('>Scaffold%d\n' % i)
        for scaffoldLine in scaffold:
            f.write('%d\t%s\t%s\n' % (scaffoldLine.GAP,\
                                      scaffoldLine.ORIENTATION,\
                                      scaffoldLine.CONTIG))
        i += 1
    f.write('\n')
    f.close()

if __name__ == '__main__':
    # Parse command line arguments for input and output files
    p = argparse.ArgumentParser(description='Make contig links')
    p.add_argument('-i', '--inputFiles', dest='inputFiles', 
        help='.tiling-files to use', required=True, nargs='+', type=str)
    p.add_argument('-o', '--output', dest='output', 
        help='Filename with output contig-links', required=True)
    p.add_argument('-n', '--nGuides', dest='nGuides', 
        help='Only use the first n guiding genomes', type=int, default=sys.maxint)
    p.add_argument('-w', '--windowSize', dest='windowSize', 
        help='Size of sliding window to use when building consensus matrix.\
        Must be 2 or greater', type=int, default=2)
    p.add_argument('-t', '--threshold', dest='threshold', 
        help='How many guiding genomes must have contigs within the window\
        size in their tiling files in order to create a contig-link', type=int,
        default=1)
    a = p.parse_args()

    # Process filenames and pruning of these
    tilingFilenames = a.inputFiles
    n = min(a.nGuides, len(tilingFilenames))
    tilingFilenamesPruned = sorted(tilingFilenames)[0:n]
    tilingsList = [parseTilingFile(t) for t in tilingFilenamesPruned]

    # Build a distance matrix from the tilings list using specified window size
    M = buildMatrix(tilingsList, a.windowSize)

    # Make a consensus distance matrix with some function reducing a list
    # of distances to one "consensus"
    MC = makeConsensusMatrix(M, a.threshold, median)

    # Grow paths based on the consensus matrix
    paths = buildPaths(MC)

    # Create scaffolds from the paths created
    scaffolds = makeScaffolds(MC, paths)

    # Write scaffolds to file
    writeScaffoldLinesToFile(a.output, scaffolds)

