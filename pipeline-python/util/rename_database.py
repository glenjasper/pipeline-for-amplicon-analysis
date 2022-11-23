#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

'''
python3 rename_database.py silva SILVA_138.1_SSURef_NR99_tax_silva.fasta
python3 rename_database.py unite uchime_reference_dataset_untrimmed_28.06.2017.fasta
'''

def format_lineage(lineage):
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

    return lineage_format

def edit_database_fasta(fasta_file, db_type = 'silva'):
    records = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        _id = record.id
        _desc = record.description
        _desc = _desc.replace(_id, '').strip()

        if db_type == 'silva':
            if _desc.startswith('Archaea;') or _desc.startswith('Bacteria;'):
                lineage = _desc.split(';')
                _lineage_format = format_lineage(lineage)
                _id = '%s;%s;' % (_id, _lineage_format)

                rec = SeqRecord(Seq(str(record.seq)), id = _id, description = '')
                records.append(rec)
        elif db_type == 'unite':
            id_seq = _id.split('|k__')[0].strip()
            id_seq = id_seq.split('|')
            _id_seq = '|'.join(id_seq[0:3])
            _id_alpha = id_seq[1]

            lineage = _id.split('|k__')[1].strip()
            lineage = 'k__%s' % lineage
            lineage = lineage.replace('k__', '').replace('p__', '').replace('c__', '').replace('o__', '').replace('f__', '').replace('g__', '').replace('s__', '')
            lineage = lineage.split(';')
            _lineage_format = format_lineage(lineage)
            # _lineage_format = _lineage_format.replace(':unidentified,', ':unidentified_%s,' % _id_alpha).replace(':unidentified;', ':unidentified_%s;' % _id_alpha)
            _id = '%s;%s;' % (_id_seq, _lineage_format)

            rec = SeqRecord(Seq(str(record.seq)), id = _id, description = '')
            records.append(rec)

    _path = os.path.dirname(fasta_file)
    _prefix, _ = os.path.splitext(os.path.basename(fasta_file))
    output_file = '%s_for_asv.fasta' % _prefix
    output_file = os.path.join(_path, output_file)

    SeqIO.write(records, output_file, 'fasta')

def main(args):
    if len(args) <= 2:
        message = 'Two arguments needed: <silva|unite> database.fasta\n'
        print(message)
    else:
        database_type = args[1]
        database_fasta = args[2]

        asvs = edit_database_fasta(database_fasta, database_type)

if __name__ == '__main__':
    main(sys.argv)
