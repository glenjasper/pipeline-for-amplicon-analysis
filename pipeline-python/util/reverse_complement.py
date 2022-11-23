#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

def get_reverse(sequence):
    return sequence.upper()[::-1]

def get_complement(sequence):
    # https://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
    complement = {'A': 'T',
                  'T': 'A',
                  'U': 'A',
                  'G': 'C',
                  'C': 'G',
                  'Y': 'R',
                  'R': 'Y',
                  'S': 'S',
                  'W': 'W',
                  'K': 'M',
                  'M': 'K',
                  'B': 'V',
                  'D': 'H',
                  'H': 'D',
                  'V': 'B',
                  'N': 'N'}
    arr_complement = []
    for base in list(sequence):
        if base in complement:
            arr_complement.append(complement[base])
        else:
            arr_complement.append(base)

    return ''.join(arr_complement)

def main(args):
    if len(args) == 1:
        message = 'Use:\n  python3 reverse_complement.py <SEQUENCE>\n'
    else:
        message = 'Reverse-complement: %s\n' % get_complement(get_reverse(args[1]))
    print(message)

if __name__ == '__main__':
    main(sys.argv)
