#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Runar Furenes'
__email__  = 'runarfu@ifi.uio.no'

from sys import argv

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

if __name__ == '__main__':
    fastas = readMultiFASTA(argv[1])
    for fasta in fastas.keys():
        l = len(fastas[fasta])
        print('%d\t%s' % (l, fasta))

