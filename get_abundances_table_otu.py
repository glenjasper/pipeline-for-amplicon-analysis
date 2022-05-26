#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import pandas as pd

def read_taxonomy_file(file, input_type = 'silva'):
    df = pd.read_csv(filepath_or_buffer = file, sep = '\t', header = None, usecols = [0, 1, 2])
    df = df.where(pd.notnull(df), '')
    # print(df)

    OTUs = {}
    for idx, row in df.iterrows():
        query_id = row[0]
        otu_id = query_id.split(';')[0].strip()

        if input_type == 'silva':
            subject_id = row[1]
            subject_title = row[2]
            lineage = subject_title.replace(subject_id, '').strip()
        elif input_type == 'rdp':
            subject_id = row[1]
            lineage = subject_id.split('tax=')[1].strip()

        if otu_id not in OTUs:
            OTUs.update({otu_id: lineage})

    return OTUs

def read_otu_file(file, dict_otus, output_file, input_type = 'silva'):
    df = pd.read_csv(filepath_or_buffer = file, sep = '\t', header = 0)
    df = df.where(pd.notnull(df), '')
    # print(df)

    arr_lineage = []
    arr_domain = []
    arr_phylum = []
    arr_class = []
    arr_order = []
    arr_family = []
    arr_genus = []
    arr_species = []
    arr_level = []
    for idx, row in df.iterrows():
        otu_id = row['#OTU ID']

        if otu_id in dict_otus:
            lineage = dict_otus[otu_id]

            if input_type == 'silva':
                lineage_split = lineage.split(';')
            elif input_type == 'rdp':
                lineage_split = lineage.replace('d:', '').replace('p:', '').replace('c:', '').replace('o:', '').replace('f:', '').replace('g:', '').replace('s:', '').replace('_', ' ').replace("\"", '')
                lineage_split = lineage_split.split(',')

            n_levels = len(lineage_split)

            level = ''
            if n_levels == 7:
                level = 'Species'
            elif n_levels == 6:
                level = 'Genus'
            elif n_levels == 5:
                level = 'Family'
            elif n_levels == 4:
                level = 'Order'
            elif n_levels == 3:
                level = 'Class'
            elif n_levels == 2:
                level = 'Phylum'
            elif n_levels == 1:
                level = 'Domain'
            arr_level.append(level)

            missing_levels = 7 - n_levels
            for n in range(missing_levels):
                lineage_split.append('')
        else:
            lineage = ''
            lineage_split = [''] * 7
            arr_level.append('Unknown')

        arr_lineage.append(lineage)
        arr_domain.append(lineage_split[0])
        arr_phylum.append(lineage_split[1])
        arr_class.append(lineage_split[2])
        arr_order.append(lineage_split[3])
        arr_family.append(lineage_split[4])
        arr_genus.append(lineage_split[5])
        arr_species.append(lineage_split[6])

    df.insert(loc = 1, column = 'Lineage', value = arr_lineage)
    df.insert(loc = 2, column = 'Domain', value = arr_domain)
    df.insert(loc = 3, column = 'Phylum', value = arr_phylum)
    df.insert(loc = 4, column = 'Class', value = arr_class)
    df.insert(loc = 5, column = 'Order', value = arr_order)
    df.insert(loc = 6, column = 'Family', value = arr_family)
    df.insert(loc = 7, column = 'Genus', value = arr_genus)
    df.insert(loc = 8, column = 'Species', value = arr_species)
    df.insert(loc = 9, column = 'Level', value = arr_level)
    # print(df)

    df.to_csv(output_file, sep = '\t', encoding = 'utf-8', index = False)

def main(args):
    if len(args) <= 4:
        message = 'Four arguments needed: <silva|rdp> file.blast, file.otu and abundance.txt\n'
        print(message)
    else:
        input_type = args[1]
        blast_file = args[2]
        otu_file = args[3]
        abundance_file = args[4]

        otus = read_taxonomy_file(blast_file, input_type)
        read_otu_file(otu_file, otus, abundance_file, input_type)

if __name__ == '__main__':
    main(sys.argv)
