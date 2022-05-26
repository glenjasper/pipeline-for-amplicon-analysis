#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import time
import traceback
import pandas as pd
from colorama import init
init()

class Parse:

    def __init__(self):
        self.VERSION = 1.0
        self.ROOT = os.path.dirname(os.path.realpath(__file__))

        # Log
        self.LOG_NAME = "log_%s_%s.log" % (os.path.splitext(os.path.basename(__file__))[0], time.strftime('%Y%m%d'))
        self.LOG_FILE = os.path.join(self.ROOT, self.LOG_NAME)

        self.OUTPUT_PATH = os.path.join(self.ROOT, 'abundance_tables')

        self.RANK_SUPERKINGDOM = 'Domain' # [1] superkingdom
        self.RANK_KINGDOM = 'Kingdom'     # [2]
        self.RANK_PHYLUM = 'Phylum'       # [3]
        self.RANK_CLASS = 'Class'         # [4]
        self.RANK_ORDER = 'Order'         # [5]
        self.RANK_FAMILY = 'Family'       # [6]
        self.RANK_GENUS = 'Genus'         # [7]
        self.RANK_SPECIES = 'Species'     # [8]

        # Fonts
        self.RED = '\033[31m'
        self.BIRED = '\033[1;91m'
        self.BIGREEN = '\033[1;92m'
        self.END = '\033[0m'

    def show_print(self, message, logs = None, showdate = True, font = None):
        msg_print = message
        msg_write = message

        if font is not None:
            msg_print = "%s%s%s" % (font, msg_print, self.END)

        if showdate is True:
            _time = time.strftime('%Y-%m-%d %H:%M:%S')
            msg_print = "%s %s" % (_time, msg_print)
            msg_write = "%s %s" % (_time, message)

        print(msg_print)
        if logs is not None:
            for log in logs:
                if log is not None:
                    with open(log, 'a', encoding = 'utf-8') as f:
                        f.write("%s\n" % msg_write)
                        f.close()

    def start_time(self):
        return time.time()

    def finish_time(self, start, message = None):
        finish = time.time()
        runtime = time.strftime("%H:%M:%S", time.gmtime(finish - start))
        if message is None:
            return runtime
        else:
            return "%s: %s" % (message, runtime)

    def check_path(self, path):
        if len(path) > 0 and os.path.exists(path):
            return True
        else:
            return False

    def create_directory(self, path):
        output = True
        try:
            if len(path) > 0 and not os.path.exists(path):
                os.makedirs(path)
        except Exception as e:
            output = False
        return output

    def read_abundance_table(self, file, prefix):
        self.show_print("Abundance file: %s" % file, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

        df = pd.read_csv(filepath_or_buffer = file, sep = '\t', header = 0)
        df = df.where(pd.notnull(df), '')
        # print(df)

        arr_taxa = [self.RANK_SUPERKINGDOM,
                    self.RANK_PHYLUM,
                    self.RANK_CLASS,
                    self.RANK_ORDER,
                    self.RANK_FAMILY,
                    self.RANK_GENUS,
                    # self.RANK_SPECIES
                    ]

        for level_taxa in arr_taxa:
            df_cond = df[df[level_taxa] != '']
            df_level = df_cond.groupby([level_taxa]).sum()
            # print(df_level)

            # Converter a valores relativos as colunas (100%)
            df_level = (100. * df_level / df_level.sum()) #.round(3).astype(str)
            # print(df_level)
            # print(df_level.dtypes)

            # -----------------------------
            # Data for the heatmap
            # -----------------------------
            # Converter o Index (gerada no groupby) em coluna
            df_level_heatmap = df_level.reset_index(level = 0)
            # print(df_level_heatmap)

            # Hellinger transform
            for colname, coltype in df_level_heatmap.dtypes.iteritems():
                if colname != level_taxa:
                    df_level_heatmap[colname] = df_level_heatmap[[colname]] ** .5
            # print(df_level_heatmap)

            # Renomeando a coluna Taxon
            df_level_heatmap = df_level_heatmap.rename({level_taxa: ''}, axis = 'columns')
            # print(df_level_heatmap)

            output_file = '%s_data-heatmap-abundances_%s.txt' % (prefix, level_taxa)
            output_file = os.path.join(self.OUTPUT_PATH, output_file.lower())

            df_level_heatmap.to_csv(output_file, sep = '\t', encoding = 'utf-8', index = False)

            self.show_print("  [Heatmap] Abundance table: %s" % output_file, [self.LOG_FILE], font = self.BIGREEN)
            self.show_print("", [self.LOG_FILE])

            # -----------------------------
            # Data for venn diagram
            # -----------------------------
            # Converter o Index (gerada no groupby) em coluna
            df_level_venn = df_level.reset_index(level = 0)

            numeric_columns = df_level_venn._get_numeric_data()
            numeric_columns[numeric_columns > 0] = 1 # presence
            # print(df_level_venn)

            output_file = '%s_data-venn-abundances_%s.txt' % (prefix, level_taxa)
            output_file = os.path.join(self.OUTPUT_PATH, output_file.lower())

            df_level_venn.to_csv(output_file, sep = '\t', encoding = 'utf-8', index = False)

            self.show_print("  [Venn] Abundance table: %s" % output_file, [self.LOG_FILE], font = self.BIGREEN)
            self.show_print("", [self.LOG_FILE])

            # -----------------------------
            # Data for the Bar-plot
            # -----------------------------
            # Converter as várias colunas em colunas verticais (precisa da coluna Index)
            df_level = df_level.stack().reset_index()
            # print(df_level)

            # Renomeando as colunas
            df_level.columns = [level_taxa, 'sample', 'value']
            # print(df_level)

            output_file = '%s_data-bar-abundances_%s.txt' % (prefix, level_taxa)
            output_file = os.path.join(self.OUTPUT_PATH, output_file.lower())

            df_level.to_csv(output_file, sep = '\t', encoding = 'utf-8', index = False)

            self.show_print("  [Bar] Abundance table: %s" % output_file, [self.LOG_FILE], font = self.BIGREEN)
            self.show_print("", [self.LOG_FILE])

def main(args):
    try:
        if len(args) <= 2:
            message = 'Two arguments needed: raw-data.txt, output_prefix\n'
            print(message)
            exit()
        else:
            input_file = args[1]
            output_prefix = args[2]

        start = oparse.start_time()
        oparse.create_directory(oparse.OUTPUT_PATH)
        oparse.LOG_FILE = os.path.join(oparse.OUTPUT_PATH, oparse.LOG_NAME)

        oparse.show_print("###########################################################", [oparse.LOG_FILE], font = oparse.BIGREEN)
        oparse.show_print("########################### RUN ###########################", [oparse.LOG_FILE], font = oparse.BIGREEN)
        oparse.show_print("###########################################################", [oparse.LOG_FILE], font = oparse.BIGREEN)

        oparse.read_abundance_table(input_file, output_prefix)

        oparse.show_print(oparse.finish_time(start, "Elapsed time"), [oparse.LOG_FILE])
        oparse.show_print("Done!", [oparse.LOG_FILE])
    except Exception as e:
        oparse.show_print("\n%s" % traceback.format_exc(), [oparse.LOG_FILE], font = oparse.RED)
        oparse.show_print(oparse.finish_time(start, "Elapsed time"), [oparse.LOG_FILE])
        oparse.show_print("Done!", [oparse.LOG_FILE])

if __name__ == '__main__':
    oparse = Parse()
    main(sys.argv)
