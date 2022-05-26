# Parameters ----------------------------------------------------------------------------
# Paths
samples_path=
util_path=
database_path=

# Primers
primers_file=illumina.primers.fa

# Taxonomy database
silva_db_fa=SILVA_138.1_SSURef_NR99_tax_silva_for_asv.fasta

# Threads
THREADS=10
# ---------------------------------------------------------------------------------------

USEARCH=$(which usearch)
VSEARCH=$(which vsearch)
FASTQC=$(which fastqc)
CUTADAPT=$(which cutadapt)
PYTHON=$(which python3)

# Read primers
primer_fwd=`awk 'NR==2' ${util_path}/${primers_file}`
primer_rev=`awk 'NR==4' ${util_path}/${primers_file}`

rc=`$PYTHON ${util_path}/reverse_complement.py ${primer_rev}`
primer_rev_rc=`echo ${rc} | awk -F' ' '{print $2}'`

# Process samples
cd ${samples_path}

for R in *_R1_001.fastq ; do
  prefix=`echo $R | awk -F'_R1_' '{print $1}'`
  s=$(cut -d. -f1 <<< "$R")

  echo "========================================================================================="
  echo "Processing sample: ${s}"
  echo "========================================================================================="

  echo ""
  echo "[Rawdata] Checking the quality of the reads"

  $FASTQC -f fastq -o ${samples_path} ${R}
  $FASTQC -f fastq -o ${samples_path} ${R/_R1/_R2}

  echo ""
done

echo ""
echo "========================================================================================="
echo "Processing all samples together"
echo "========================================================================================="

echo ""
echo "Merge all samples into one fastq file"

# The program automatically recognizes the R2
$USEARCH -fastq_mergepairs *_R1_001.fastq \
         -fastqout all_samples_merged.fq \
         -relabel @

echo ""
echo "[Merged] Checking the quality of the reads"

$FASTQC -f fastq -o ${samples_path} all_samples_merged.fq

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
         -db ${util_path}/${primers_file} \
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
         --fastq_maxee 0.5 \
         --fastq_minlen 300 \
         --eeout \
         --fastqout all_samples_filtered.fq \
         --fastaout all_samples_filtered.fa \
         --fasta_width 0 \
         --fastq_qmax 42

echo ""
echo "[Filtered] Checking the quality of the reads"

$FASTQC -f fastq -o ${samples_path} all_samples_filtered.fq

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
         --threads $THREADS \
         --db ASVs.fa \
         --id 0.99 \
         --otutabout ASV_counts.txt

########################################
# With SILVA
########################################

echo ""
echo "Assigning taxonomy"

$USEARCH -sintax ASVs.fa \
         -db ${database_path}/${silva_db_fa} \
         -tabbedout ASV_taxonomy.txt \
         -strand both \
         -sintax_cutoff 0.8

echo ""
echo "Get table of abundances of ASVs with taxonomy"

$PYTHON ${util_path}/get_abundances_table_asv.py \
        ASV_taxonomy.txt \
        ASV_counts.txt \
        abundance_table_asv.csv

echo "Abundance file: abundance_table_asv.csv"

echo ""
echo "done."
