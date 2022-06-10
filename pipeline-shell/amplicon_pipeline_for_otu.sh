# Parameters ----------------------------------------------------------------------------
# Paths
samples_path=
database_path=
util_path=
output_path=

# Taxonomy database (files must be in database_path | database_type can be silva, rdp or unite | database_bin only for OTUs)
database_type=
database_fasta=
database_bin=

# Primers file (file must be in database_path)
primers_file=

# Threads
threads=10

# For quality filtering (maxee default: 0.8 | filter_maxlen is optional)
filter_maxee=0.8
filter_minlen=
filter_maxlen=

# For clustering (default: 97)
cluster_identity=97

# For taxonomic alignment (default: 97)
blast_identity=97
# ---------------------------------------------------------------------------------------

USEARCH=$(which usearch)
VSEARCH=$(which vsearch)
CUTADAPT=$(which cutadapt)
BLASTN=$(which blastn)
FASTQC=$(which fastqc)
PYTHON=$(which python3)

# Read primers
primer_fwd=`awk 'NR==2' ${database_path}/${primers_file}`
primer_rev=`awk 'NR==4' ${database_path}/${primers_file}`

rc=`$PYTHON ${util_path}/reverse_complement.py ${primer_rev}`
primer_rev_rc=`echo ${rc} | awk -F' ' '{print $2}'`

# echo "Checking FASTQ format version for one file"
# $VSEARCH --fastq_chars $(ls -1 ${samples_path}/*.fastq | head -1)

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
  echo "Merge paired-end sequence reads into one sequence"

  $VSEARCH --fastq_mergepairs ${R} \
           --threads $threads \
           --reverse ${R/_R1/_R2} \
           --fastqout ${prefix}.merged.fq \
           --fastq_eeout

  echo ""
  echo "[Merged] Checking the quality of the reads"

  $FASTQC -f fastq -o . ${prefix}.merged.fq

  # echo ""
  # echo "Calculate quality statistics"

  # $VSEARCH --fastq_eestats ${prefix}.merged.fq \
  #          --output ${prefix}.stats

  echo ""
  echo "Verification of the primers"

  echo ""
  echo "Extraction of a subsample of 1000 reads"

  $USEARCH -fastx_subsample ${prefix}.merged.fq \
           -sample_size 1000 \
           -fastqout ${prefix}.merged_subset_1000.fq

  echo ""
  echo "Verification of the position of the primers"

  $USEARCH -search_oligodb ${prefix}.merged_subset_1000.fq \
           -db ${database_path}/${primers_file} \
           -strand both \
           -userout ${prefix}.merged_primer_hits.txt \
           -userfields query+qlo+qhi+qstrand

  echo ""
  echo "Removal of the forward-primer (5')"

  $CUTADAPT -g ${primer_fwd} \
            --discard-untrimmed \
            -o ${prefix}.trimmed_pfwd.fq \
            ${prefix}.merged.fq

  echo ""
  echo "Removal of the reverse-primer (3')"

  $CUTADAPT -a ${primer_rev_rc} \
            --discard-untrimmed \
            -o ${prefix}.trimmed_prev.fq \
            ${prefix}.trimmed_pfwd.fq

  echo ""
  echo "Quality filtering"

  $VSEARCH --fastq_filter ${prefix}.trimmed_prev.fq \
           --fastq_maxee ${filter_maxee} \
           --fastq_minlen ${filter_minlen} \
           --fastq_maxlen ${filter_maxlen} \
           --eeout \
           --fastqout ${prefix}.filtered.fq \
           --fastaout ${prefix}.filtered.fa \
           --fasta_width 0 \
           --relabel `basename ${prefix}`.

  echo ""
  echo "[Filtered] Checking the quality of the reads"

  $FASTQC -f fastq -o . ${prefix}.filtered.fq

  echo ""
done

echo "Sum of sequences in each sample:" $(cat *.filtered.fa | grep -c "^>")

# At this point there should be one fasta file for each sample
# It should be quality filtered and dereplicated.

echo ""
echo "========================================================================================="
echo "Processing all samples together"
echo "========================================================================================="

echo ""
echo "Merge all samples"

cat *.filtered.fa > all.fa

echo ""
echo "Dereplicate across samples and remove singletons"

$VSEARCH --derep_fulllength all.fa \
         --minuniquesize 2 \
         --sizein \
         --sizeout \
         --fasta_width 0 \
         --uc all.dereplicated.uc \
         --output all.dereplicated.fa

echo "Unique non-singleton sequences:" $(grep -c "^>" all.dereplicated.fa)

echo ""
echo "Precluster at 97% before chimera detection"

$VSEARCH --cluster_size all.dereplicated.fa \
         --threads $threads \
         --id `bc -l <<< "scale=2; ${cluster_identity}/100"` \
         --strand plus \
         --sizein \
         --sizeout \
         --fasta_width 0 \
         --uc all.preclustered.uc \
         --centroids all.preclustered.fa

echo "Unique sequences after preclustering:" $(grep -c "^>" all.preclustered.fa)

echo ""
echo "De novo chimera detection"

$VSEARCH --uchime_denovo all.preclustered.fa \
         --sizein \
         --sizeout \
         --fasta_width 0 \
         --nonchimeras all.denovo.nonchimeras.fa

echo "Unique sequences after de novo chimera detection:" $(grep -c "^>" all.denovo.nonchimeras.fa)

echo ""
echo "Reference chimera detection"

$VSEARCH --uchime_ref all.denovo.nonchimeras.fa \
         --threads $threads \
         --db ${database_path}/${database_fasta} \
         --sizein \
         --sizeout \
         --fasta_width 0 \
         --nonchimeras all.ref.nonchimeras.fa

echo "Unique sequences after reference-based chimera detection:" $(grep -c "^>" all.ref.nonchimeras.fa)

echo ""
echo "Extract all non-chimeric, non-singleton sequences, dereplicated"

$PYTHON ${util_path}/map.py \
        all.dereplicated.fa \
        all.preclustered.uc \
        all.ref.nonchimeras.fa \
        all.nonchimeras.dereplicated.fa

echo "Unique non-chimeric, non-singleton sequences:" $(grep -c "^>" all.nonchimeras.dereplicated.fa)

echo ""
echo "Extract all non-chimeric, non-singleton sequences in each sample"

$PYTHON ${util_path}/map.py \
        all.fa \
        all.dereplicated.uc \
        all.nonchimeras.dereplicated.fa \
        all.nonchimeras.fa

echo "Sum of unique non-chimeric, non-singleton sequences in each sample:" $(grep -c "^>" all.nonchimeras.fa)

echo ""
echo "Cluster at 97% and relabel with OTU_n, generate OTU table"

$VSEARCH --cluster_size all.nonchimeras.fa \
         --threads $threads \
         --id `bc -l <<< "scale=2; ${cluster_identity}/100"` \
         --strand plus \
         --sizein \
         --sizeout \
         --fasta_width 0 \
         --relabel OTU_ \
         --uc all.clustered.uc \
         --centroids all.otus.fa \
         --otutabout all.otutab.txt \
         --biomout all.otutab.biom

echo "Number of OTUs:" $(grep -c "^>" all.otus.fa)

########################################
# With SILVA
########################################

echo ""
echo "Identification of OTUs using BLAST"

$BLASTN -db ${database_path}/${database_bin} \
        -query all.otus.fa \
        -perc_identity ${blast_identity} \
        -qcov_hsp_perc 90.0 \
        -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qcovs" \
        -out taxonomy.blast

echo "Blast file: taxonomy.blast"

echo ""
echo "Get table of abundances of OTUs with taxonomy"

$PYTHON ${util_path}/get_abundances_table_otu.py \
        ${database_type} \
        taxonomy.blast \
        all.otutab.txt \
        abundance_table_otu.csv

echo "Abundance file: abundance_table_otu.csv"

echo ""
echo "done."
