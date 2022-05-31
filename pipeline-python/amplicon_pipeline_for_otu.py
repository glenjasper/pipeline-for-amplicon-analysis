#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import re
import sys
import time
import argparse
import traceback
import subprocess
import configparser
from Bio import SeqIO
from colorama import init
init()

def menu(args):
    parser = argparse.ArgumentParser(description = opipe.PIPELINE, epilog = "Thank you!")
    parser.add_argument("-c", "--config_file", metavar = "FILE", required = True, help = "Configuration file")
    parser.add_argument("--version", action = "version", version = "%s %s" % ('%(prog)s', opipe.VERSION))
    args = parser.parse_args()

    if opipe.check_path(args.config_file):
        _path = os.path.dirname(args.config_file)
        if _path is None or _path == "":
            _path = opipe.ROOT

        _file = os.path.join(_path, os.path.basename(args.config_file))
        opipe.SETTINGS_FILE = _file
        opipe.SETTINGS_FILE_NAME = os.path.basename(_file)

        opipe.read_keys()
    else:
        opipe.show_print("File '%s' doesn't exist" % args.config_file, showdate = False, font = opipe.YELLOW)
        exit()

class Pipeline:

    def __init__(self):
        self.VERSION = 1.0
        self.ROOT = os.path.dirname(os.path.realpath(__file__))
        self.PIPELINE = "Pipeline for the analysis of 16s rRNA amplicons, through OTUs (Operational Taxonomic Unit)"

        # Config file
        self.SETTINGS_FILE = None
        self.SETTINGS_FILE_NAME = None

        # Log
        self.LOG_NAME = "log_%s_%s.log" % (os.path.splitext(os.path.basename(__file__))[0], time.strftime('%Y%m%d'))
        self.LOG_FILE = None

        # Utilities folder
        self.UTILITIES_PATH = 'utilities'

        # Programs | Scripts
        self.PROGRAM_VSEARCH = None
        self.PROGRAM_USEARCH = None
        self.PROGRAM_CUTADAPT = None
        self.PROGRAM_BLASTN = None
        self.PROGRAM_FASTQC = 'fastqc'
        self.PROGRAM_RC = 'reverse_complement.py'
        self.PROGRAM_MAP = 'map.py'
        self.PROGRAM_ABUNDANCE_TABLE = 'get_abundances_table_otu.py'

        # Key parameters
        self.KEY_SAMPLES_PATH = None
        self.KEY_DATABASE_PATH = None
        self.KEY_UTIL_PATH = None
        self.KEY_OUTPUT_PATH = None
        self.KEY_PRIMERS_FILE = None
        self.KEY_DATABASE_FASTA = None
        self.KEY_DATABASE_BIN = None
        self.KEY_THREADS = None
        self.KEY_PYTHON_VERSION = None
        self.KEY_OS_TYPE = None

        # Sections
        self.SECTION_PARAMETERS = "PARAMETERS"

        # Keys
        self.PARAMETER_SAMPLES_PATH = "SAMPLES_PATH"
        self.PARAMETER_DATABASE_PATH = "DATABASE_PATH"
        self.PARAMETER_UTIL_PATH = "UTIL_PATH"
        self.PARAMETER_OUTPUT_PATH = "OUTPUT_PATH"
        self.PARAMETER_PRIMERS_FILE = "PRIMERS_FILE"
        self.PARAMETER_DATABASE_FASTA = "DATABASE_FASTA"
        self.PARAMETER_DATABASE_BIN = "DATABASE_BIN"
        self.PARAMETER_THREADS = "THREADS"
        self.PARAMETER_PYTHON_VERSION = "PYTHON_VERSION"
        self.PARAMETER_OS_TYPE = "OS_TYPE"

        self.OS_TYPE_GNULINUX = "gnulinux"
        self.OS_TYPE_WINDOWS = "win"

        # Fonts
        self.RED = '\033[31m'
        self.YELLOW = '\033[33m'
        self.ICYAN = '\033[96m'
        self.IGREEN = '\033[92m'
        self.BIGREEN = '\033[1;92m'
        self.END = '\033[0m'

    def read_settings(self, file, section, key):
        config = configparser.ConfigParser()
        config.read(file)
        try:
            value = config.get(section, key)
        except Exception as e:
            value = ""
        return value

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

    def count_sequences(self, fasta_file):
        n = len([1 for line in open(fasta_file, 'r') if line.startswith('>')])
        return n

    def read_keys(self):
        self.show_print("[Checking parameters and programs to use]", showdate = False, font = self.ICYAN)

        self.KEY_SAMPLES_PATH = self.read_settings(self.SETTINGS_FILE, self.SECTION_PARAMETERS, self.PARAMETER_SAMPLES_PATH)
        self.KEY_DATABASE_PATH = self.read_settings(self.SETTINGS_FILE, self.SECTION_PARAMETERS, self.PARAMETER_DATABASE_PATH)
        self.KEY_UTIL_PATH = self.read_settings(self.SETTINGS_FILE, self.SECTION_PARAMETERS, self.PARAMETER_UTIL_PATH)
        self.KEY_OUTPUT_PATH = self.read_settings(self.SETTINGS_FILE, self.SECTION_PARAMETERS, self.PARAMETER_OUTPUT_PATH)
        self.KEY_DATABASE_FASTA = self.read_settings(self.SETTINGS_FILE, self.SECTION_PARAMETERS, self.PARAMETER_DATABASE_FASTA)
        self.KEY_DATABASE_BIN = self.read_settings(self.SETTINGS_FILE, self.SECTION_PARAMETERS, self.PARAMETER_DATABASE_BIN)
        self.KEY_PRIMERS_FILE = self.read_settings(self.SETTINGS_FILE, self.SECTION_PARAMETERS, self.PARAMETER_PRIMERS_FILE)
        self.KEY_THREADS = self.read_settings(self.SETTINGS_FILE, self.SECTION_PARAMETERS, self.PARAMETER_THREADS)
        self.KEY_PYTHON_VERSION = self.read_settings(self.SETTINGS_FILE, self.SECTION_PARAMETERS, self.PARAMETER_PYTHON_VERSION)
        self.KEY_OS_TYPE = self.read_settings(self.SETTINGS_FILE, self.SECTION_PARAMETERS, self.PARAMETER_OS_TYPE)

        # Output path
        if not self.KEY_OUTPUT_PATH:
            self.show_print("[WARNING] Value of parameter '%s' not specified" % (self.PARAMETER_OUTPUT_PATH.lower()), showdate = False, font = self.YELLOW)
            exit()
        else:
            if not self.create_directory(self.KEY_OUTPUT_PATH):
                self.show_print("[WARNING] Could not create '%s' directory (parameter '%s')\n" % (self.KEY_OUTPUT_PATH, self.PARAMETER_OUTPUT_PATH.lower()), showdate = False, font = self.YELLOW)
                exit()
            else:
                self.create_directory(self.KEY_OUTPUT_PATH)
                self.LOG_FILE = os.path.join(self.KEY_OUTPUT_PATH, self.LOG_NAME)

        # Samples path
        if not self.KEY_SAMPLES_PATH:
            self.show_print("[WARNING] Value of parameter '%s' not specified" % (self.PARAMETER_SAMPLES_PATH.lower()), showdate = False, font = self.YELLOW)
            exit()
        else:
            if not self.check_path(self.KEY_SAMPLES_PATH):
                self.show_print("[WARNING] Path '%s' of parameter '%s' doesn't exist" % (self.KEY_SAMPLES_PATH, self.PARAMETER_SAMPLES_PATH.lower()), showdate = False, font = self.YELLOW)
                exit()
            else:
                are_there_files = False
                for subdir, dirs, files in os.walk(self.KEY_SAMPLES_PATH):
                    for file in files:
                        if re.search('[_][Rr][1][_](\w|[-])+\.([Ff][Aa][Ss][Tt][Qq]|[Ff][Qq])$', file):
                            are_there_files = True
                            break
                if not are_there_files:
                    self.show_print("[WARNING] Path '%s' doesn't contain any FASTQ file with the format <part1>_R1_<part2>.fastq" % (self.KEY_SAMPLES_PATH), showdate = False, font = self.YELLOW)
                    exit()

        # Database path
        if not self.KEY_DATABASE_PATH:
            self.show_print("[WARNING] Value of parameter '%s' not specified" % (self.PARAMETER_DATABASE_PATH.lower()), showdate = False, font = self.YELLOW)
            exit()
        else:
            if not self.check_path(self.KEY_DATABASE_PATH):
                self.show_print("[WARNING] Path '%s' of parameter '%s' doesn't exist" % (self.KEY_DATABASE_PATH, self.PARAMETER_DATABASE_PATH.lower()), showdate = False, font = self.YELLOW)
                exit()

        # Util path
        if not self.KEY_UTIL_PATH:
            self.show_print("[WARNING] Value of parameter '%s' not specified" % (self.PARAMETER_UTIL_PATH.lower()), showdate = False, font = self.YELLOW)
            exit()
        else:
            if not self.check_path(self.KEY_UTIL_PATH):
                self.show_print("[WARNING] Path '%s' of parameter '%s' doesn't exist" % (self.KEY_UTIL_PATH, self.PARAMETER_UTIL_PATH.lower()), showdate = False, font = self.YELLOW)
                exit()

        # Database fasta file
        if not self.KEY_DATABASE_FASTA:
            self.show_print("[WARNING] Value of parameter '%s' not specified" % (self.PARAMETER_DATABASE_FASTA.lower()), showdate = False, font = self.YELLOW)
            exit()
        else:
            _file = os.path.join(self.KEY_DATABASE_PATH, self.KEY_DATABASE_FASTA)
            if not os.path.isfile(_file):
                self.show_print("[WARNING] File '%s' of parameter '%s' doesn't exist" % (_file, self.PARAMETER_DATABASE_FASTA.lower()), showdate = False, font = self.YELLOW)
                exit()

        # Database binary file
        if not self.KEY_DATABASE_BIN:
            self.show_print("[WARNING] Value of parameter '%s' not specified" % (self.PARAMETER_DATABASE_BIN.lower()), showdate = False, font = self.YELLOW)
            exit()
        else:
            _file = os.path.join(self.KEY_DATABASE_PATH, self.KEY_DATABASE_BIN)
            if not os.path.isfile('%s.nhr' % _file):
                self.show_print("[WARNING] File '%s.*' of parameter '%s' doesn't exist" % (_file, self.PARAMETER_DATABASE_BIN.lower()), showdate = False, font = self.YELLOW)
                exit()

        # Primers file
        if not self.KEY_PRIMERS_FILE:
            self.show_print("[WARNING] Value of parameter '%s' not specified" % (self.PARAMETER_PRIMERS_FILE.lower()), showdate = False, font = self.YELLOW)
            exit()
        else:
            _file = os.path.join(self.KEY_DATABASE_PATH, self.KEY_PRIMERS_FILE)
            if not os.path.isfile(_file):
                self.show_print("[WARNING] File '%s' of parameter '%s' doesn't exist" % (_file, self.PARAMETER_PRIMERS_FILE.lower()), showdate = False, font = self.YELLOW)
                exit()

        # Threads
        if not self.KEY_THREADS:
            self.show_print("[WARNING] Value of parameter '%s' not specified" % (self.PARAMETER_THREADS.lower()), showdate = False, font = self.YELLOW)
            exit()
        else:
            if (not self.KEY_THREADS.isdigit()) or (int(self.KEY_THREADS) == 0):
                self.show_print("[WARNING] Value '%s' of parameter '%s' is not a positive integer" % (self.KEY_THREADS, self.PARAMETER_THREADS.lower()), showdate = False, font = self.YELLOW)
                exit()

        # OS type
        if not self.KEY_OS_TYPE:
            self.show_print("[WARNING] Value of parameter '%s' not specified" % (self.PARAMETER_OS_TYPE.lower()), showdate = False, font = self.YELLOW)
            exit()
        else:
            if not self.KEY_OS_TYPE.lower() in [self.OS_TYPE_GNULINUX, self.OS_TYPE_WINDOWS]:
                self.show_print("[WARNING] You must specify some value for the '%s' parameter: %s (for GNU/Linux) or %s (for Windows)" % (self.PARAMETER_PYTHON_VERSION.lower(), self.OS_TYPE_GNULINUX, self.OS_TYPE_WINDOWS), showdate = False, font = self.YELLOW)
                exit()

        # Python version
        if not self.KEY_PYTHON_VERSION:
            self.show_print("[WARNING] Value of parameter '%s' not specified" % (self.PARAMETER_PYTHON_VERSION.lower()), showdate = False, font = self.YELLOW)
            exit()
        else:
            if not self.KEY_PYTHON_VERSION.lower() in ['python', 'python3']:
                self.show_print("[WARNING] You must specify some value for the '%s' parameter: python3 (for Python 3.X in GNU/Linux) or python (for Python3.X in Windows)" % (self.PARAMETER_PYTHON_VERSION.lower()), showdate = False, font = self.YELLOW)
                exit()
            else:
                if self.KEY_OS_TYPE == self.OS_TYPE_GNULINUX and self.KEY_PYTHON_VERSION == 'python':
                    self.show_print("[WARNING] You are using a GNU/Linux platform, and you specified the parameter '%s' as '%s', you may have problems, we recommend using '%s'" % (self.PARAMETER_PYTHON_VERSION.lower(), self.KEY_PYTHON_VERSION, 'python3'), showdate = False, font = self.YELLOW)

        self.UTILITIES_PATH = os.path.join(self.ROOT, self.UTILITIES_PATH)
        self.UTILITIES_PATH = os.path.join(self.UTILITIES_PATH, self.KEY_OS_TYPE.lower())

        if self.KEY_OS_TYPE.lower() == self.OS_TYPE_GNULINUX:
            self.PROGRAM_VSEARCH = 'vsearch'
            self.PROGRAM_USEARCH = 'usearch'
            self.PROGRAM_CUTADAPT = 'cutadapt'
            self.PROGRAM_BLASTN = 'blastn'
        elif self.KEY_OS_TYPE.lower() == self.OS_TYPE_WINDOWS:
            self.PROGRAM_VSEARCH = 'vsearch.exe'
            self.PROGRAM_USEARCH = 'usearch.exe'
            self.PROGRAM_CUTADAPT = 'cutadapt.exe'
            self.PROGRAM_BLASTN = 'blastn.exe'

        prog_vsearch = os.path.join(self.UTILITIES_PATH, self.PROGRAM_VSEARCH)
        prog_usearch = os.path.join(self.UTILITIES_PATH, self.PROGRAM_USEARCH)
        prog_cutadapt = os.path.join(self.UTILITIES_PATH, self.PROGRAM_CUTADAPT)
        prog_blastn = os.path.join(self.UTILITIES_PATH, self.PROGRAM_BLASTN)

        fastqc_path = os.path.join(os.path.dirname(self.UTILITIES_PATH), 'common', 'FastQC')
        prog_fastqc = os.path.join(fastqc_path, self.PROGRAM_FASTQC)

        if self.KEY_OS_TYPE.lower() == self.OS_TYPE_GNULINUX:
            os.chmod(prog_vsearch, int('755', base = 8))
            os.chmod(prog_usearch, int('755', base = 8))
            os.chmod(prog_cutadapt, int('755', base = 8))
            os.chmod(prog_blastn, int('755', base = 8))
            os.chmod(prog_fastqc, int('755', base = 8))

        self.check_version('%s --version' % prog_vsearch, self.PROGRAM_VSEARCH)
        self.check_version('%s --version' % prog_usearch, self.PROGRAM_USEARCH)
        self.check_version('%s --version' % prog_cutadapt, self.PROGRAM_CUTADAPT)
        self.check_version('%s -version' % prog_blastn, self.PROGRAM_BLASTN)

        if self.KEY_OS_TYPE.lower() == self.OS_TYPE_GNULINUX:
            self.check_version('%s --version' % prog_fastqc, self.PROGRAM_FASTQC)
        elif self.KEY_OS_TYPE.lower() == self.OS_TYPE_WINDOWS:
            jar1_path = os.path.join(fastqc_path, 'sam-1.103.jar')
            jar2_path = os.path.join(fastqc_path, 'jbzip2-0.9.jar')
            program_path_fqc = 'java -Xmx250m -Dfastqc.show_version=true -Djava.awt.headless=true -classpath %s;%s;%s uk.ac.babraham.FastQC.FastQCApplication' % (fastqc_path, jar1_path, jar2_path)
            self.check_version(program_path_fqc, self.PROGRAM_FASTQC)

        self.show_print("Ok", showdate = False, font = self.IGREEN)

    def check_version(self, cmd, program):
        p = subprocess.Popen(cmd, shell = True, stdin = None, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (checkStdout, checkStderr) = p.communicate()
        checkStdout = checkStdout.decode('utf-8').strip()
        checkStderr = checkStderr.decode('utf-8').strip()

        if not checkStdout:
            self.show_print("[Check %s version]" % (program), showdate = False, font = self.ICYAN)
            self.show_print("There are problems with the '%s' program, check your installation" % program, showdate = False, font = self.YELLOW)
            exit()
        # else:
        #     self.show_print("[Check %s version]" % (program), showdate = False, font = self.ICYAN)
        #     self.show_print("%s" % checkStdout, showdate = False, font = self.IGREEN)

    def get_cmd_information(self, command):
        self.show_print("Command information:", [self.LOG_FILE])
        for item in command:
            self.show_print("  %s" % item, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

    def run_program(self, program, command, success_words = None, **kwargs):
        _extra_info = kwargs.get('extra_info') if kwargs.get('extra_info') else ''

        self.show_print("---------------------------------------------------------------------------------", [self.LOG_FILE], font = self.IGREEN)
        self.show_print("[Run %s] %s" % (program, _extra_info), [self.LOG_FILE], font = self.BIGREEN)
        self.show_print("---------------------------------------------------------------------------------", [self.LOG_FILE], font = self.IGREEN)
        self.get_cmd_information(command)
        start = self.start_time()

        _command = " ".join(command)
        try:
            self.show_print("Running...", [self.LOG_FILE])
            p = subprocess.Popen(_command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
        except Exception as e:
            self.show_print("Error %s while executing command %s" % (e, _command), [self.LOG_FILE], font = self.YELLOW)

        successful = False
        output_run = ''
        for line in iter(p.stdout.readline, b''):
            _line = line.decode('ISO-8859-1').rstrip()
            # _line = line.decode('utf-8').rstrip()

            if _line.startswith(success_words) or success_words in _line:
                successful = True

                if program == self.PROGRAM_RC:
                    primer_rev_rc = _line.split()
                    primer_rev_rc = primer_rev_rc[1]
                    output_run = primer_rev_rc

            self.show_print(_line, [self.LOG_FILE])

        if program in [self.PROGRAM_MAP, self.PROGRAM_BLASTN, self.PROGRAM_ABUNDANCE_TABLE]:
            successful = True

        if not successful:
            err_msg = "ERROR executing %s!\nCheck the command: %s" % (program, _command)
            self.show_print(err_msg, [self.LOG_FILE], font = self.YELLOW)
            self.show_print(self.finish_time(start, "Interrupted after"), [self.LOG_FILE], font = self.YELLOW)
            self.show_print("", [self.LOG_FILE])
            sys.exit(1)
        else:
            if program == self.PROGRAM_FASTQC and self.KEY_OS_TYPE.lower() == self.OS_TYPE_WINDOWS:
                input_path = os.path.dirname(command[1])
                input_file = os.path.basename(command[1])

                if input_path != self.KEY_OUTPUT_PATH:
                    _prefix, _ = os.path.splitext(input_file)

                    zip_file = '%s_fastqc.zip' % _prefix
                    html_file = '%s_fastqc.html' % _prefix

                    out_zip_file = os.path.join(self.KEY_OUTPUT_PATH, zip_file)
                    if os.path.isfile(out_zip_file):
                        os.remove(out_zip_file)

                    out_html_file = os.path.join(self.KEY_OUTPUT_PATH, html_file)
                    if os.path.isfile(out_html_file):
                        os.remove(out_html_file)

                    os.rename(os.path.join(input_path, zip_file), out_zip_file)
                    os.rename(os.path.join(input_path, html_file), out_html_file)

        self.show_print(self.finish_time(start, "Elapsed time"), [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

        return output_run

    def run_get_primers(self):
        primers_file = os.path.join(self.KEY_DATABASE_PATH, self.KEY_PRIMERS_FILE)

        primer_fwd = ''
        primer_rev = ''
        for index, record in enumerate(SeqIO.parse(primers_file, 'fasta')):
            if index == 0:
                # forward-primer
                primer_fwd = record.seq
            elif index == 1:
                # reverse-primer
                primer_rev = record.seq

        arr_cmd = ['%s %s' % (self.KEY_PYTHON_VERSION, os.path.join(self.KEY_UTIL_PATH, self.PROGRAM_RC)),
                   '%s' % primer_rev]

        words = 'Reverse-complement'

        primer_rev_rc = self.run_program(program = self.PROGRAM_RC,
                                         command = arr_cmd,
                                         success_words = words)

        return primer_fwd, primer_rev_rc

    def run_vsearch(self, params, step = None, extra_info = None):
        arr_cmd = []
        words = ''
        if step == 'fastq_mergepairs':
            arr_cmd = ['%s' % os.path.join(self.UTILITIES_PATH, self.PROGRAM_VSEARCH),
                       '--fastq_mergepairs %s' % params['r1'],
                       '--reverse %s' % params['r2'],
                       '--threads %s' % self.KEY_THREADS,
                       '--fastqout %s' % params['output'],
                       '--fastq_eeout']

            words = 'Statistics of merged reads'
        elif step == 'fastq_filter':
            arr_cmd = ['%s' % os.path.join(self.UTILITIES_PATH, self.PROGRAM_VSEARCH),
                       '--fastq_filter %s' % params['input'],
                       '--fastq_maxee %s' % params['fastq_maxee'],
                       '--fastq_minlen %s' % params['fastq_minlen'],
                       '--eeout',
                       '--fastqout %s' % params['fastqout'],
                       '--fastaout %s' % params['fastaout'],
                       '--fasta_width %s' % params['fasta_width']]

            words = 'Reading input file 100%'
        elif step == 'derep_fulllength':
            arr_cmd = ['%s' % os.path.join(self.UTILITIES_PATH, self.PROGRAM_VSEARCH),
                       '--derep_fulllength %s' % params['input'],
                       '--minuniquesize %s' % params['minuniquesize'],
                       '--strand %s' % params['strand'],
                       '--output %s' % params['output'],
                       '--sizeout',
                       '--uc %s' % params['uc'],
                       '--relabel %s' % params['relabel'],
                       '--fasta_width %s' % params['fasta_width']]

            words = 'Writing FASTA output file 100%'
        elif step == 'derep_fulllength_all':
            arr_cmd = ['%s' % os.path.join(self.UTILITIES_PATH, self.PROGRAM_VSEARCH),
                       '--derep_fulllength %s' % params['input'],
                       '--minuniquesize %s' % params['minuniquesize'],
                       '--output %s' % params['output'],
                       '--sizein',
                       '--sizeout',
                       '--uc %s' % params['uc'],
                       '--fasta_width %s' % params['fasta_width']]

            words = 'Writing FASTA output file 100%'
        elif step == 'cluster_size':
            arr_cmd = ['%s' % os.path.join(self.UTILITIES_PATH, self.PROGRAM_VSEARCH),
                       '--cluster_size %s' % params['input'],
                       '--threads %s' % self.KEY_THREADS,
                       '--id %s' % params['id'],
                       '--strand %s' % params['strand'],
                       '--sizein',
                       '--sizeout',
                       '--fasta_width %s' % params['fasta_width'],
                       '--uc %s' % params['uc'],
                       '--centroids %s' % params['centroids']]

            words = 'Clustering 100%'
        elif step == 'uchime_denovo':
            arr_cmd = ['%s' % os.path.join(self.UTILITIES_PATH, self.PROGRAM_VSEARCH),
                       '--uchime_denovo %s' % params['input'],
                       '--sizein',
                       '--sizeout',
                       '--fasta_width %s' % params['fasta_width'],
                       '--nonchimeras %s' % params['nonchimeras']]

            words = 'Detecting chimeras 100%'
        elif step == 'uchime_ref':
            arr_cmd = ['%s' % os.path.join(self.UTILITIES_PATH, self.PROGRAM_VSEARCH),
                       '--uchime_ref %s' % params['input'],
                       '--threads %s' % self.KEY_THREADS,
                       '--db %s' % params['db'],
                       '--sizein',
                       '--sizeout',
                       '--fasta_width %s' % params['fasta_width'],
                       '--nonchimeras %s' % params['nonchimeras']]

            words = 'Detecting chimeras 100%'
        elif step == 'cluster_size_otu_table':
            arr_cmd = ['%s' % os.path.join(self.UTILITIES_PATH, self.PROGRAM_VSEARCH),
                       '--cluster_size %s' % params['input'],
                       '--threads %s' % self.KEY_THREADS,
                       '--id %s' % params['id'],
                       '--strand %s' % params['strand'],
                       '--sizein',
                       '--sizeout',
                       '--fasta_width %s' % params['fasta_width'],
                       '--relabel %s' % params['relabel'],
                       '--uc %s' % params['uc'],
                       '--centroids %s' % params['centroids'],
                       '--otutabout %s' % params['otutabout'],
                       '--biomout %s' % params['biomout']]

            words = 'Clustering 100%'

        self.run_program(program = self.PROGRAM_VSEARCH,
                         command = arr_cmd,
                         success_words = words,
                         extra_info = extra_info)

    def run_usearch(self, params, step = None, extra_info = None):
        arr_cmd = []
        words = ''
        if step == 'fastx_subsample':
            arr_cmd = ['%s' % os.path.join(self.UTILITIES_PATH, self.PROGRAM_USEARCH),
                       '-fastx_subsample %s' % params['input'],
                       '-sample_size %s' % params['sample_size'],
                       '-fastqout %s' % params['output']]

            words = 'Sampling'
        elif step == 'search_oligodb':
            arr_cmd = ['%s' % os.path.join(self.UTILITIES_PATH, self.PROGRAM_USEARCH),
                       '-search_oligodb %s' % params['input'],
                       '-db %s' % params['db'],
                       '-strand %s' % params['strand'],
                       '-userout %s' % params['output'],
                       '-userfields %s' % params['userfields']]

            words = 'matched'

        self.run_program(program = self.PROGRAM_USEARCH,
                         command = arr_cmd,
                         success_words = words,
                         extra_info = extra_info)

    def run_cutadapt(self, params, step = None, extra_info = None):
        arr_cmd = []
        words = ''
        if step == 'forward':
            arr_cmd = ['%s' % os.path.join(self.UTILITIES_PATH, self.PROGRAM_CUTADAPT),
                       '-g %s' % params['f_primer'],
                       '--discard-untrimmed',
                       '-o %s' % params['output'],
                       '%s' % params['input']]

            words = 'Overview of removed sequences'
        elif step == 'reverse':
            arr_cmd = ['%s' % os.path.join(self.UTILITIES_PATH, self.PROGRAM_CUTADAPT),
                       '-a %s' % params['r_primer_rc'],
                       '--discard-untrimmed',
                       '-o %s' % params['output'],
                       '%s' % params['input']]

            words = 'Overview of removed sequences'

        self.run_program(program = self.PROGRAM_CUTADAPT,
                         command = arr_cmd,
                         success_words = words,
                         extra_info = extra_info)

    def run_map(self, params, extra_info = None):
        arr_cmd = ['%s %s' % (self.KEY_PYTHON_VERSION, os.path.join(self.KEY_UTIL_PATH, self.PROGRAM_MAP)),
                   '%s' % params['fasta1'],
                   '%s' % params['uc'],
                   '%s' % params['fasta2'],
                   '%s' % params['output']]

        primer_rev_rc = self.run_program(program = self.PROGRAM_MAP,
                                         command = arr_cmd,
                                         extra_info = extra_info)

    def run_blastn(self, params, extra_info = None):
        arr_cmd = ['%s' % os.path.join(self.UTILITIES_PATH, self.PROGRAM_BLASTN),
                   '-db %s' % params['db'],
                   '-query %s' % params['query'],
                   '-perc_identity %s' % params['perc_identity'],
                   '-qcov_hsp_perc %s' % params['qcov_hsp_perc'],
                   '-outfmt "%s"' % params['outfmt'],
                   '-out %s' % params['out']]

        self.run_program(program = self.PROGRAM_BLASTN,
                         command = arr_cmd,
                         extra_info = extra_info)

    def run_get_abundances_table(self, params, extra_info = None):
        arr_cmd = ['%s %s' % (self.KEY_PYTHON_VERSION, os.path.join(self.KEY_UTIL_PATH, self.PROGRAM_ABUNDANCE_TABLE)),
                   '%s' % params['db_type'],
                   '%s' % params['blast_file'],
                   '%s' % params['otutab_file'],
                   '%s' % params['output']]

        primer_rev_rc = self.run_program(program = self.PROGRAM_ABUNDANCE_TABLE,
                                         command = arr_cmd,
                                         extra_info = extra_info)

    def run_merge_all(self, fasta_file, prefix):
        self.show_print("---------------------------------------------------------------------------------", [self.LOG_FILE], font = self.IGREEN)
        self.show_print("[Merge all samples]", [self.LOG_FILE], font = self.BIGREEN)
        self.show_print("---------------------------------------------------------------------------------", [self.LOG_FILE], font = self.IGREEN)

        with open(fasta_file, 'w') as fw:
            for subdir, dirs, files in os.walk(self.KEY_OUTPUT_PATH):
                for file in files:
                    if file.endswith('.dereplicated.fa') and not file.startswith(prefix):
                        with open(os.path.join(subdir, file), 'r') as fr:
                            for line in fr:
                                fw.write(line)
                        fr.close()
        fw.close()

        self.show_print("  Output file: %s" % fasta_file, [self.LOG_FILE])
        self.show_print("  Number of sequences: %s" % self.count_sequences(fasta_file), [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

    def run_fastqc(self, params, extra_info = None):
        fastqc_path = os.path.join(os.path.dirname(self.UTILITIES_PATH), 'common', 'FastQC')

        if self.KEY_OS_TYPE.lower() == self.OS_TYPE_WINDOWS:
            jar1_path = os.path.join(fastqc_path, 'sam-1.103.jar')
            jar2_path = os.path.join(fastqc_path, 'jbzip2-0.9.jar')

            program_path_fqc = 'java -Xmx250m -classpath %s;%s;%s uk.ac.babraham.FastQC.FastQCApplication' % (fastqc_path, jar1_path, jar2_path)

            arr_cmd = ['%s' % program_path_fqc,
                       '%s' % params['input']]
        elif self.KEY_OS_TYPE.lower() == self.OS_TYPE_GNULINUX:
            program_path_fqc = os.path.join(fastqc_path, self.PROGRAM_FASTQC)

            arr_cmd = ['%s' % program_path_fqc,
                       '-f %s' % 'fastq',
                       '-o %s' % self.KEY_OUTPUT_PATH,
                       '%s' % params['input']]

        words = 'Analysis complete'

        self.run_program(program = self.PROGRAM_FASTQC,
                         command = arr_cmd,
                         success_words = words,
                         extra_info = extra_info)

    def run_pipeline(self):
        primer_fwd, primer_rev_rc = self.run_get_primers()

        for subdir, dirs, files in os.walk(self.KEY_SAMPLES_PATH):
            for file in files:
                if re.search('[_][Rr][1][_](\w|[-])+\.([Ff][Aa][Ss][Tt][Qq]|[Ff][Qq])$', file):
                    fastq_r1_file = os.path.join(subdir, file)
                    fastq_r2_file = file.replace('_R1_', '_R2_')
                    fastq_r2_file = os.path.join(subdir, fastq_r2_file)
                    prefix = file.split('_R1_')[0]

                    #################################################################################
                    # [Rawdata] Checking the quality of the reads
                    #################################################################################

                    info = '%s: Checking the quality of the reads [R1]' % prefix
                    params = {'input': fastq_r1_file}
                    self.run_fastqc(params, extra_info = info)

                    info = '%s: Checking the quality of the reads [R2]' % prefix
                    params = {'input': fastq_r2_file}
                    self.run_fastqc(params, extra_info = info)

                    #################################################################################
                    # Merge paired-end sequence reads into one sequence
                    #################################################################################

                    output_merged = '%s.merged.fq' % prefix
                    output_merged = os.path.join(self.KEY_OUTPUT_PATH, output_merged)
                    info = '%s: Merge paired-end sequence reads' % prefix

                    params = {'r1': fastq_r1_file,
                              'r2': fastq_r2_file,
                              'output': output_merged}

                    self.run_vsearch(params, step = 'fastq_mergepairs', extra_info = info)

                    #################################################################################
                    # [Merged] Checking the quality of the reads
                    #################################################################################

                    info = '%s: Checking the quality of the reads' % prefix
                    params = {'input': output_merged}
                    self.run_fastqc(params, extra_info = info)

                    #################################################################################
                    # Verification of the primers
                    # Extraction of a subsample of 1000 reads
                    #################################################################################

                    output_subset = '%s.merged_subset_1000.fq' % prefix
                    output_subset = os.path.join(self.KEY_OUTPUT_PATH, output_subset)
                    info = '%s: Extraction of a subsample of 1000 reads' % prefix

                    params = {'input': output_merged,
                              'sample_size': '1000',
                              'output': output_subset}

                    self.run_usearch(params, step = 'fastx_subsample', extra_info = info)

                    #################################################################################
                    # Verification of the position of the primers
                    #################################################################################

                    primers_file = os.path.join(self.KEY_DATABASE_PATH, self.KEY_PRIMERS_FILE)
                    output_oligodb = '%s.merged_primer_hits.txt' % prefix
                    output_oligodb = os.path.join(self.KEY_OUTPUT_PATH, output_oligodb)
                    info = '%s: Verification of the position of the primers' % prefix

                    params = {'input': output_subset,
                              'db': primers_file,
                              'strand': 'both',
                              'userfields': 'query+qlo+qhi+qstrand',
                              'output': output_oligodb}

                    self.run_usearch(params, step = 'search_oligodb', extra_info = info)

                    #################################################################################
                    # Removal of the forward-primer (5')
                    #################################################################################

                    output_trimmed_pfwd = '%s.trimmed_pfwd.fq' % prefix
                    output_trimmed_pfwd = os.path.join(self.KEY_OUTPUT_PATH, output_trimmed_pfwd)
                    info = '%s: Removal of the forward-primer (5\')' % prefix

                    params = {'input': output_merged,
                              'f_primer': primer_fwd,
                              'output': output_trimmed_pfwd}

                    self.run_cutadapt(params, step = 'forward', extra_info = info)

                    #################################################################################
                    # Removal of the reverse-primer (3')
                    #################################################################################

                    output_trimmed_prev = '%s.trimmed_prev.fq' % prefix
                    output_trimmed_prev = os.path.join(self.KEY_OUTPUT_PATH, output_trimmed_prev)
                    info = '%s: Removal of the reverse-primer (3\')' % prefix

                    params = {'input': output_trimmed_pfwd,
                              'r_primer_rc': primer_rev_rc,
                              'output': output_trimmed_prev}

                    self.run_cutadapt(params, step = 'reverse', extra_info = info)

                    #################################################################################
                    # Quality filtering
                    #################################################################################

                    output_filter_fq = '%s.filtered.fq' % prefix
                    output_filter_fq = os.path.join(self.KEY_OUTPUT_PATH, output_filter_fq)
                    output_filter_fa = '%s.filtered.fa' % prefix
                    output_filter_fa = os.path.join(self.KEY_OUTPUT_PATH, output_filter_fa)
                    info = '%s: Quality filtering' % prefix

                    params = {'input': output_trimmed_prev,
                              'fastq_maxee': '0.5',
                              'fastq_minlen': '300',
                              'fasta_width': '0',
                              'fastqout': output_filter_fq,
                              'fastaout': output_filter_fa}

                    self.run_vsearch(params, step = 'fastq_filter', extra_info = info)

                    #################################################################################
                    # [Filtered] Checking the quality of the reads
                    #################################################################################

                    info = '%s: Checking the quality of the reads' % prefix
                    params = {'input': output_filter_fq}
                    self.run_fastqc(params, extra_info = info)

                    #################################################################################
                    # Dereplicate at sample level and relabel with sample.n
                    #################################################################################

                    output_fulllength = '%s.dereplicated.fa' % prefix
                    output_fulllength = os.path.join(self.KEY_OUTPUT_PATH, output_fulllength)
                    output_fulllength_uc = '%s.dereplicated.uc' % prefix
                    output_fulllength_uc = os.path.join(self.KEY_OUTPUT_PATH, output_fulllength_uc)
                    info = '%s: Dereplicate at sample level and relabel with sample.n' % prefix

                    params = {'input': output_filter_fa,
                              'minuniquesize': '2',
                              'strand': 'plus',
                              'fasta_width': '0',
                              'relabel': '%s.' % prefix,
                              'uc': output_fulllength_uc,
                              'output': output_fulllength}

                    self.run_vsearch(params, step = 'derep_fulllength', extra_info = info)

        #################################################################################
        # Merge all samples
        #################################################################################

        prefix = 'all'
        all_fasta_file = '%s.fa' % prefix
        all_fasta_file = os.path.join(self.KEY_OUTPUT_PATH, all_fasta_file)

        self.run_merge_all(all_fasta_file, prefix)

        #################################################################################
        # Dereplicate across samples and remove singletons
        #################################################################################

        output_dereplicated_all = '%s.dereplicated.fa' % prefix
        output_dereplicated_all = os.path.join(self.KEY_OUTPUT_PATH, output_dereplicated_all)
        output_dereplicated_all_uc = '%s.dereplicated.uc' % prefix
        output_dereplicated_all_uc = os.path.join(self.KEY_OUTPUT_PATH, output_dereplicated_all_uc)
        info = 'Dereplicate across samples and remove singletons'

        params = {'input': all_fasta_file,
                  'minuniquesize': '2',
                  'fasta_width': '0',
                  'uc': output_dereplicated_all_uc,
                  'output': output_dereplicated_all}

        self.run_vsearch(params, step = 'derep_fulllength_all', extra_info = info)

        self.show_print("Unique non-singleton sequences: %s" % self.count_sequences(output_dereplicated_all), [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

        #################################################################################
        # Precluster at 97% before chimera detection
        #################################################################################

        output_preclustered_all = '%s.preclustered.fa' % prefix
        output_preclustered_all = os.path.join(self.KEY_OUTPUT_PATH, output_preclustered_all)
        output_preclustered_all_uc = '%s.preclustered.uc' % prefix
        output_preclustered_all_uc = os.path.join(self.KEY_OUTPUT_PATH, output_preclustered_all_uc)
        info = 'Precluster at 97% before chimera detection'

        params = {'input': output_dereplicated_all,
                  'id': '0.97',
                  'strand': 'plus',
                  'fasta_width': '0',
                  'uc': output_preclustered_all_uc,
                  'centroids': output_preclustered_all}

        self.run_vsearch(params, step = 'cluster_size', extra_info = info)

        self.show_print("Unique sequences after preclustering: %s" % self.count_sequences(output_preclustered_all), [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

        #################################################################################
        # De novo chimera detection
        #################################################################################

        output_nonchimeras_dn = '%s.denovo.nonchimeras.fa' % prefix
        output_nonchimeras_dn = os.path.join(self.KEY_OUTPUT_PATH, output_nonchimeras_dn)
        info = 'De novo chimera detection'

        params = {'input': output_preclustered_all,
                  'fasta_width': '0',
                  'nonchimeras': output_nonchimeras_dn}

        self.run_vsearch(params, step = 'uchime_denovo', extra_info = info)

        self.show_print("Unique sequences after de novo chimera detection: %s" % self.count_sequences(output_nonchimeras_dn), [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

        #################################################################################
        # Reference chimera detection
        #################################################################################

        database_fasta = os.path.join(self.KEY_DATABASE_PATH, self.KEY_DATABASE_FASTA)
        output_nonchimeras_ref = '%s.ref.nonchimeras.fa' % prefix
        output_nonchimeras_ref = os.path.join(self.KEY_OUTPUT_PATH, output_nonchimeras_ref)
        info = 'Reference chimera detection'

        params = {'input': output_nonchimeras_dn,
                  'db': database_fasta,
                  'strand': 'plus',
                  'fasta_width': '0',
                  'nonchimeras': output_nonchimeras_ref}

        self.run_vsearch(params, step = 'uchime_ref', extra_info = info)

        self.show_print("Unique sequences after reference-based chimera detection: %s" % self.count_sequences(output_nonchimeras_ref), [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

        #################################################################################
        # Extract all non-chimeric, non-singleton sequences, dereplicated
        #################################################################################

        output_map_1 = '%s.nonchimeras.dereplicated.fa' % prefix
        output_map_1 = os.path.join(self.KEY_OUTPUT_PATH, output_map_1)
        info = 'Extract all non-chimeric, non-singleton sequences, dereplicated'

        params = {'fasta1': output_dereplicated_all,
                  'uc': output_preclustered_all_uc,
                  'fasta2': output_nonchimeras_ref,
                  'output': output_map_1}

        self.run_map(params, extra_info = info)

        self.show_print("Unique non-chimeric, non-singleton sequences: %s" % self.count_sequences(output_map_1), [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

        #################################################################################
        # Extract all non-chimeric, non-singleton sequences in each sample
        #################################################################################

        output_map_2 = '%s.nonchimeras.fa' % prefix
        output_map_2 = os.path.join(self.KEY_OUTPUT_PATH, output_map_2)
        info = 'Extract all non-chimeric, non-singleton sequences in each sample'

        params = {'fasta1': all_fasta_file,
                  'uc': output_dereplicated_all_uc,
                  'fasta2': output_map_1,
                  'output': output_map_2}

        self.run_map(params, extra_info = info)

        self.show_print("Sum of unique non-chimeric, non-singleton sequences in each sample: %s" % self.count_sequences(output_map_1), [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

        #################################################################################
        # Cluster at 97% and relabel with OTU_n, generate OTU table
        #################################################################################

        output_cluster_fa = '%s.otus.fa' % prefix
        output_cluster_fa = os.path.join(self.KEY_OUTPUT_PATH, output_cluster_fa)
        output_cluster_uc = '%s.clustered.uc' % prefix
        output_cluster_uc = os.path.join(self.KEY_OUTPUT_PATH, output_cluster_uc)
        output_cluster_otutab = '%s.otutab.txt' % prefix
        output_cluster_otutab = os.path.join(self.KEY_OUTPUT_PATH, output_cluster_otutab)
        output_cluster_biom = '%s.otutab.biom' % prefix
        output_cluster_biom = os.path.join(self.KEY_OUTPUT_PATH, output_cluster_biom)
        info = 'Cluster at 97% and relabel with OTU_n, generate OTU table'

        params = {'input': output_map_2,
                  'id': '0.97',
                  'strand': 'plus',
                  'fasta_width': '0',
                  'relabel': 'OTU_',
                  'uc': output_cluster_uc,
                  'centroids': output_cluster_fa,
                  'otutabout': output_cluster_otutab,
                  'biomout': output_cluster_biom}

        self.run_vsearch(params, step = 'cluster_size_otu_table', extra_info = info)

        self.show_print("Number of OTUs: %s" % self.count_sequences(output_cluster_fa), [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

        #################################################################################
        # Identification of OTUs using BLAST
        #################################################################################

        database_bin = os.path.join(self.KEY_DATABASE_PATH, self.KEY_DATABASE_BIN)
        output_blastn = 'taxonomy.blast'
        output_blastn = os.path.join(self.KEY_OUTPUT_PATH, output_blastn)
        info = 'Identification of OTUs using BLAST'

        params = {'db': database_bin,
                  'query': output_cluster_fa,
                  'perc_identity': '97.0',
                  'qcov_hsp_perc': '90.0',
                  'outfmt': '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qcovs',
                  'out': output_blastn}

        self.run_blastn(params, extra_info = info)

        self.show_print("Blast file: %s" % output_blastn, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

        #################################################################################
        # Get table of abundances of OTUs with taxonomy
        #################################################################################

        output_abundances_table = 'abundance_table_otu.csv'
        output_abundances_table = os.path.join(self.KEY_OUTPUT_PATH, output_abundances_table)
        info = 'Get table of abundances of OTUs with taxonomy'

        params = {'db_type': 'silva',
                  'blast_file': output_blastn,
                  'otutab_file': output_cluster_otutab,
                  'output': output_abundances_table}

        self.run_get_abundances_table(params, extra_info = info)

        self.show_print("Abundance file: %s" % output_abundances_table, [self.LOG_FILE], font = opipe.IGREEN)
        self.show_print("", [self.LOG_FILE])

def main(args):
    try:
        menu(args)
        start = opipe.start_time()

        opipe.show_print("#################################################################################", [opipe.LOG_FILE], font = opipe.BIGREEN)
        opipe.show_print("###################################### RUN ######################################", [opipe.LOG_FILE], font = opipe.BIGREEN)
        opipe.show_print("#################################################################################", [opipe.LOG_FILE], font = opipe.BIGREEN)

        opipe.run_pipeline()

        opipe.show_print(opipe.finish_time(start, "Elapsed time [Total]"), [opipe.LOG_FILE])
        opipe.show_print("Done!", [opipe.LOG_FILE])
    except Exception as e:
        opipe.show_print("\n%s" % traceback.format_exc(), [opipe.LOG_FILE], font = opipe.RED)
        opipe.show_print(opipe.finish_time(start, "Elapsed time [Total]"), [opipe.LOG_FILE])
        opipe.show_print("Done!", [opipe.LOG_FILE])

if __name__ == '__main__':
    opipe = Pipeline()
    main(sys.argv)
