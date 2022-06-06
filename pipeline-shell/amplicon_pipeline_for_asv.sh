# Parameters ----------------------------------------------------------------------------
# Paths
samples_path=
database_path=
util_path=
output_path=

# Taxonomy database (files must be in database_path)
database_fasta=

# Primers file (file must be in database_path)
primers_file=

# Threads
threads=10

# For quality filtering (maxee default: 0.8 | filter_maxlen is optional)
filter_maxee=0.8
filter_minlen=
filter_maxlen=

# High identity to count ASVs (default: 99)
high_identity_asv=99

# For taxonomic assignment (default: 0.8)
sintax_cutoff=0.8
# ---------------------------------------------------------------------------------------

USEARCH=$(which usearch)
VSEARCH=$(which vsearch)
CUTADAPT=$(which cutadapt)
FASTQC=$(which fastqc)
PYTHON=$(which python3)

# Read primers
primer_fwd=`awk 'NR==2' ${database_path}/${primers_file}`
primer_rev=`awk 'NR==4' ${database_path}/${primers_file}`

rc=`$PYTHON ${util_path}/reverse_complement.py ${primer_rev}`
primer_rev_rc=`echo ${rc} | awk -F' ' '{print $2}'`

mkdir -p ${output_path}

# Process samples
cd ${output_path}

for R in ${samples_path}/*_R1_001.fastq ; do
  prefix=`echo $R | awk -F'_R1_' '{print $1}'`
  prefix=./`basename ${prefix}`
  s=$(cut -d. -f1 <<< "$R")

  echo "========================================================================================="
  echo "Processing sample: `basename ${s}`"
  echo "========================================================================================="

  echo ""
  echo "[Rawdata] Checking the quality of the reads"

  $FASTQC -f fastq -o . ${R}
  $FASTQC -f fastq -o . ${R/_R1/_R2}

  echo ""
done

echo ""
echo "========================================================================================="
echo "Processing all samples together"
echo "========================================================================================="

echo ""
echo "Merge all samples into one fastq file"

# The program automatically recognizes the R2
$USEARCH -fastq_mergepairs ${samples_path}/*_R1_001.fastq \
         -fastqout all_samples_merged.fq \
         -relabel @

echo ""
echo "[Merged] Checking the quality of the reads"

$FASTQC -f fastq -o . all_samples_merged.fq

echo ""
echo "Verification of the primers"

echo ""
echo "Extraction of a subsample of 1000 reads"

$USEARCH -fastx_subsample all_samples_merged.fq \
         -sample_size 1000 \
         -fastqout subset_1000_samples_merged.fq

echo ""
echo "Verification of the position of the primers"

$USEARCH -search_oligodb subset_1000_samples_merged.fq \
         -db ${database_path}/${primers_file} \
         -strand both \
         -userout primer_hits.txt \
         -userfields query+qlo+qhi+qstrand

echo ""
echo "Removal of the forward-primer (5')"

$CUTADAPT -g ${primer_fwd} \
          --discard-untrimmed \
          -o all_samples_trimmed_pfwd.fq \
          all_samples_merged.fq

echo ""
echo "Removal of the reverse-primer (3')"

$CUTADAPT -a ${primer_rev_rc} \
          --discard-untrimmed \
          -o all_samples_trimmed_prev.fq \
          all_samples_trimmed_pfwd.fq

echo ""
echo "Quality filtering"

$VSEARCH --fastq_filter all_samples_trimmed_prev.fq \
         --fastq_maxee ${filter_maxee} \
         --fastq_minlen ${filter_minlen} \
         --fastq_maxlen ${filter_maxlen} \
         --eeout \
         --fastqout all_samples_filtered.fq \
         --fastaout all_samples_filtered.fa \
         --fasta_width 0 \
         --fastq_qmax 45

echo ""
echo "[Filtered] Checking the quality of the reads"

$FASTQC -f fastq -o . all_samples_filtered.fq

echo ""
echo "Dereplicate reads"

$VSEARCH --derep_fulllength all_samples_filtered.fa \
         --strand plus \
         --sizein \
         --sizeout \
         --fasta_width 0 \
         --uc all_samples_dereplicated.uc \
         --output all_samples_dereplicated.fa

echo ""
echo "Generating ASVs"

# Already disregards the chimeras
$USEARCH -unoise3 all_samples_dereplicated.fa \
         -zotus ASVs.fa \
         -tabbedout unoise3.txt

# unoise3.txt: Text filename which reports the processing done for each sequence, 
# e.g. if it is classified as noisy or chimeric.

echo ""
echo "Generating a count table"

sed -i 's/Zotu/ASV_/' ASVs.fa
$VSEARCH --usearch_global all_samples_filtered.fa \
         --threads $threads \
         --db ASVs.fa \
         --id `bc -l <<< "scale=2; ${high_identity_asv}/100"` \
         --otutabout ASV_counts.txt

########################################
# With SILVA
########################################

echo ""
echo "Assigning taxonomy"

$USEARCH -sintax ASVs.fa \
         -db ${database_path}/${database_fasta} \
         -tabbedout ASV_taxonomy.txt \
         -strand both \
         -sintax_cutoff ${sintax_cutoff}

echo ""
echo "Get table of abundances of ASVs with taxonomy"

$PYTHON ${util_path}/get_abundances_table_asv.py \
        ASV_taxonomy.txt \
        ASV_counts.txt \
        abundance_table_asv.csv

echo "Abundance file: abundance_table_asv.csv"

echo ""
echo "done."
