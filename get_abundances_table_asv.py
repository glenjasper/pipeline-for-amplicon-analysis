#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import pandas as pd

def read_taxonomy_file(file):
    df = pd.read_csv(filepath_or_buffer = file, sep = '\t', header = None, usecols = [0, 3])
    df = df.where(pd.notnull(df), '')
    # print(df)

    ASVs = {}
    for idx, row in df.iterrows():
        asv_id = row[0]
        lineage = row[3].strip()

        if asv_id not in ASVs:
            ASVs.update({asv_id: lineage})

    return ASVs

def read_asv_file(file, dict_asvs, output_file):
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
        asv_id = row['#OTU ID']

        if asv_id in dict_asvs:
            lineage = dict_asvs[asv_id]
            _lineage = lineage.replace('d:', '').replace('p:', '').replace('c:', '').replace('o:', '').replace('f:', '').replace('g:', '').replace('s:', '').replace('_', ' ').replace("\"", '')
            lineage_split = _lineage.split(',')
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
    df.rename(columns = {'#OTU ID': '#ASV ID'}, inplace = True)
    # print(df)

    df.to_csv(output_file, sep = '\t', encoding = 'utf-8', index = False)

def main(args):
    if len(args) <= 3:
        message = 'Three arguments needed: taxonomy.txt, asv.txt and abundance.txt\n'
        print(message)
    else:
        taxonomy_file = args[1]
        asv_file = args[2]
        abundance_file = args[3]

        asvs = read_taxonomy_file(taxonomy_file)
        read_asv_file(asv_file, asvs, abundance_file)

if __name__ == '__main__':
    main(sys.argv)
