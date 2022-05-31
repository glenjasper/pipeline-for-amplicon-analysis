#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

'''
python3 rename_silva.py SILVA_138.1_SSURef_NR99_tax_silva.fasta
'''

def edit_silva_fasta_file(fasta_file):
    records = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        _id = record.id
        _desc = record.description
        _desc = _desc.replace(_id, '').strip()

        if _desc.startswith('Archaea;') or _desc.startswith('Bacteria;'):
            lineage = _desc.split(';')

            lineage_format = 'tax='
            for index, tax in enumerate(lineage):
                index = 7 - index
                tax = tax.replace(' ', '_')
                if index == 7:
                    tax = 'd:%s,' % tax
                elif index == 6:
                    tax = 'p:%s,' % tax
                elif index == 5:
                    tax = 'c:%s,' % tax
                elif index == 4:
                    tax = 'o:%s,' % tax
                elif index == 3:
                    tax = 'f:%s,' % tax
                elif index == 2:
                    tax = 'g:%s,' % tax
                elif index == 1:
                    tax = tax.replace(',', '.')
                    tax = 's:%s,' % tax
                lineage_format += tax

            if lineage_format[-1] == ',':
                lineage_format = lineage_format[0:len(lineage_format) - 1]

            _id = '%s;%s;' % (_id, lineage_format)

            rec = SeqRecord(Seq(str(record.seq)), id = _id, description = '')
            records.append(rec)

    _path = os.path.dirname(fasta_file)
    _prefix, _ = os.path.splitext(os.path.basename(fasta_file))
    output_file = '%s_for_asv.fasta' % _prefix
    output_file = os.path.join(_path, output_file)

    SeqIO.write(records, output_file, 'fasta')

def main(args):
    if len(args) <= 1:
        message = 'One argument needed: SILVA_138.1_SSURef_NR99_tax_silva.fasta\n'
        print(message)
    else:
        silva_fasta = args[1]

        asvs = edit_silva_fasta_file(silva_fasta)

if __name__ == '__main__':
    main(sys.argv)
