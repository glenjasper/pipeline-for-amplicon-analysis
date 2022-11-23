#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

def read_fasta2_file(fasta_file):
    heads = {}
    with open(fasta_file, 'r') as fr:
        for line in fr:
            line = line.strip()
            if line.startswith('>'):
                _id = str(line).replace('>', '').split(';')[0]
                heads.update({_id: 1})
    fr.close()

    return heads

def read_uc_file(uc_file, dictionary):
    with open(uc_file, 'r') as fr:
        for line in fr:
            line = line.strip().split('\t')

            hit_type = line[0].split(';')[0]
            query_seq = line[8].split(';')[0]
            target_seq = line[9].split(';')[0]

            if hit_type == 'H' and target_seq in dictionary:
                dictionary.update({query_seq: 1})
    fr.close()

    return dictionary

def read_fasta1_file(fasta_file, dictionary, output_fasta):
    with open(fasta_file, 'r') as fr, open(output_fasta, 'w') as fw:
        save_seq = False
        for line in fr:
            line = line.strip()
            if line.startswith('>'):
                _id = line.replace('>', '').split(';')[0]
                if _id in dictionary:
                    sequence_id = '%s\n' % line
                    fw.write(sequence_id)
                    save_seq = True
                else:
                    dictionary.update({_id: None})
                    save_seq = False
            else:
                if save_seq:
                    sequence = '%s\n' % str(line)
                    fw.write(sequence)
    fr.close()
    fw.close()

def main(args):
    if len(args) <= 4:
        message = 'Four arguments needed: fasta1, uc, fasta2 and outfasta\n'
        print(message)
    else:
        fasta_file1 = args[1]
        uc_file = args[2]
        fasta_file2 = args[3]
        output_file = args[4]

        dict_heads = read_fasta2_file(fasta_file2)
        read_uc_file(uc_file, dict_heads)
        read_fasta1_file(fasta_file1, dict_heads, output_file)

if __name__ == '__main__':
    main(sys.argv)
