Pipeline for amplicon analysis
======================
[![License](https://poser.pugx.org/badges/poser/license.svg)](./LICENSE)

_Pipeline_ para análise metagenômico de _amplicon_ 16S rRNA utilizando as abordabens de agrupamento por OTUs (Operational Taxonomic Unit) e geração de ASVs (Amplicon Sequence Variant). Para ambas as abordagens, se disponibilizam scripts em [_Python 3_](#pipeline-em-python-3) que podem ser executados nas plataformas **GNU/Linux** e **Windows** e [_Shell script_](#pipeline-em-shell-script) que pode ser executado apenas em **GNU/Linux** (ou Mac OS).

## Recursos

- [Banco de dados SILVA](#banco-de-dados-silva)
    - [Baixar e criar os binários do banco de dados](#baixar-e-criar-os-binários-do-banco-de-dados)
    - [Formatação do banco de dados para uso de ASVs](#formatação-do-banco-de-dados-para-uso-de-asvs)
- [Pré-requisitos](#pré-requisitos)
    - [Programas](#programas)
    - [Bibliotecas Python](#bibliotecas-python)
- [Scripts](#scripts)
- [Pipeline](#pipeline)
    - [Arquivo de configuração](#arquivo-de-configuração)
- [Como usar o pipeline](#como-usar-o-pipeline)
- [Credits](#credits)
- [Author](#author)
- [Organization](#organization)
- [License](#license)

## Banco de dados SILVA

### Baixar e criar os binários do banco de dados
Para análise 16S rRNA pode-se usar o banco de dados [SILVA 138 SSU NR](https://www.arb-silva.de/no_cache/download/archive/current/Exports). Pode baixar o arquivo FASTA e criar os binários com os seguintes comandos:

```sh
  wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
  gzip -d SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz

  makeblastdb -in SILVA_138.1_SSURef_NR99_tax_silva.fasta \
              -dbtype nucl \
              -out silva_db
```

### Formatação do banco de dados para uso de ASVs
O banco de dados SILVA não pode ser utilizado diretamente no _pipeline_ com ASVs, entretanto pode ser formatado com o script _rename_silva.py_, da seguinte forma.

```sh
  python3 rename_silva.py SILVA_138.1_SSURef_NR99_tax_silva.fasta
```

## Pré-requisitos

### Programas
O _pipeline_ precisa dos seguintes requisitos que estão inclusos no código fonte (pasta **pipeline-python/bin**).
- [VSEARCH](https://github.com/torognes/vsearch)
- [USEARCH 32-bit](https://drive5.com/usearch)
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)
  - Para o correto funcionamento do FastQC préviamente precisa instalar o [JDK de JAVA](https://www.oracle.com/br/java/technologies/javase/jdk11-archive-downloads.html)
- [Python 3](https://www.python.org)

### Bibliotecas Python

```sh
  sudo pip3 install pandas
  sudo pip3 install biopython
  sudo pip3 install colorama
  sudo pip3 install cutadapt # Only for GNU/Linux
```

## _Scripts_

- **util/map.py**: _Script_ para mapear leituras _non-singletons_ e _non-chimeras_ (adaptado de [map.pl](https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline)).
- **util/rename_database.py**: _Script_ para renomear os cabeçalhos do arquivo FASTA dos banco de dados [SILVA 138 SSU NR](https://www.arb-silva.de/no_cache/download/archive/current/Exports) e [UNITE](https://unite.ut.ee/repository.php) para serem utilizados com o _pipeline_ com ASVs.
- **util/reverse_complement.py**: _Script_ para obter a reversa-complementar de uma sequência (_forward-primer_).
- **util/get_abundances_table_otu.py**: _Script_ para obter a tabela de abundâncias dos OTUs com dados taxonômicos.
- **util/get_abundances_table_asv.py**: _Script_ para obter a tabela de abundâncias dos ASVs com dados taxonômicos.

## _Pipeline_

Se desenvolveu um _script_ em _Python 3_ para ambas as unidades de medidas (OTUs e ASVs). Este _script_ pode ser utilizado tanto na plataforma **GNU/Linux** (ou Mac OS) quanto no **Windows**. O _script_ já têm incorporado os [programas](#programas) requeridos na pasta **pipeline-python/bin** (arquivos: _common.zip_, _gnulinux.zip_ e _win.zip_), os quais serão descompactados no primeiro uso do _Pipeline_. No entanto, o usuário precisa configurar os parâmetros do arquivo [**config.txt**](#arquivo-de-configuração).

#### _Pipeline_ para as abordagens de ASVs e OTUs
- **pipeline-python/amplicon_pipeline.py**: Fluxo (_pipeline_) para a geração de uma tabela de abundâncias de ASVs ou OTUs com dados taxonômicos utilizando o banco de dados SILVA, a partir de dados de sequenciamentos de _amplicon_ 16S rRNA.

#### Arquivo de configuração

O arquivo **config.txt** contém os parâmetros necessários para o funcionamento do _pipeline_. Exemplo de configuração:

```sh
  [PARAMETERS]

  # Type of approach to analysis (ASV or OTU)
  approach_type = asv

  # Paths
  samples_path = D:\my_workstation\amplicon\dataset
  database_path = D:\my_workstation\amplicon\database
  output_path = D:\my_workstation\amplicon\results-otu

  # Taxonomy database (files must be in database_path | database_type can be silva, rdp or unite | database_bin only for OTUs)
  database_type = silva
  database_fasta = SILVA_138.1_SSURef_NR99_tax_silva.fasta
  database_bin = silva_db

  # Primers file (file must be in database_path)
  primers_file = illumina.primers.fa

  # Multiprocessing
  threads = 10

  # Platform type (gnulinux: for GNU/Linux | win: for Windows)
  platform_type = win

  # Python version (python3: for Python 3.x in GNU/Linux | python: for Python 3.x in Windows)
  python_version = python

  # [ASVs/OTUs] For quality filtering (maxee default: 0.8 | filter_maxlen is optional)
  filter_maxee = 0.8
  filter_minlen = 350
  filter_maxlen = 

  # [Only OTUs] For clustering (default: 97)
  cluster_identity = 97

  # [Only OTUs] For taxonomic alignment (default: 97)
  blast_identity = 97

  # [Only ASVs] High identity to count ASVs (default: 99)
  high_identity_asv = 99

  # [Only ASVs] For taxonomic assignment (default: 0.8)
  sintax_cutoff = 0.8
```
Descrição de parâmetros, que também se aplicam para os [_Shell Script_](#exemplo-de-configuração-os-parâmetros-internos-dos-scripts):

| Parameter             | Description |
| --------------------- | ----------- |
| **approach_type**     | Tipo de abotrdagem (ASV ou OTU). |
| **samples_path**      | Caminho absoluto da pasta que contem os arquivos FASTQ. |
| **database_path**     | Caminho absoluto da pasta que contem o banco de dados SILVA (FASTA e binários). |
| **output_path**       | Caminho absoluto da pasta de saída. |
| **database_type**     | Tipo de banco de dados taxonômico: **silva**, **rdp** e **unite**. |
| **database_fasta**    | Nome do arquivo FASTA do banco de dados SILVA. |
| **database_bin**      | Prefixo dos arquivos binários do banco de dados SILVA. (Usado apenas com **OTUs**) |
| **primers_file**      | Nome do arquivo FASTA que contém os _primers_ _forward_ e _reverse_ (o arquivo debe estar em **database_path**). |
| **threads**           | Número de _threads_ para multiprocessamento. |
| **platform_type**     | Tipo de plataforma: **gnulinux** para GNU/Linux ou **win** para Windows. |
| **python_version**    | Tipo de executable do Python 3: **python3** geralmente usado em GNU/Linux ou **python** geralmente usado em Windows. |
| **filter_maxee**      | Máximo valor do erro esperado (E_max) das leituras. Se descartam as leituras com > E_max (_default_: 0.8). |
| **filter_minlen**     | Tamanho mínimo de leitura (minlen). Se descartam as leituras com tamanho < minlen. |
| **filter_maxlen**     | Tamanho máximo de leitura (maxlen). Se descartam as leituras com tamanho > maxlen. |
| **cluster_identity**  | Valor de identidade do alinhamento a ser usado para a geração dos _clusters_ (_default_: 97). (Usado apenas com **OTUs**) |
| **blast_identity**    | Valor de identidade (blastn) para a atribuição taxonômica com o banco de dados taxonômico (_default_: 97). (Usado apenas com **OTUs**) |
| **high_identity_asv** | Valor de identidade para o mapeamento dos ASVs (_default_: 99). (Usado apenas com **ASVs**) |
| **sintax_cutoff**     | Valor do _cutoff_ para a atribuição taxonômica dos ASVs com o banco de dados taxonômico (_default_: 0.8). (Usado apenas com **ASVs**) |

> **Nota**: Do exemplo de configuração, o arquivo FASTA **illumina.primers.fa** deve conter primeiro o _forward-primer_ e depois o _reverse-primer_. É extremamente importante saber os _primers_ que foram utilizados na amplificação (PCR) dos seus dados. Aqui mostramos um exemplo do conteúdo do arquivo de configuração com os _primers_ universais **341F** e **806R**:
```sh
>341F 
CCTACGGGRSGCAGCAG
>806R
GGACTACHVGGGTWTCTAAT
```

> **Nota**: Os parâmetros **database_bin**, **cluster_identity** e **blast_identity** são utilizados apenas para a abordagem com OTUs. Os parâmetros **high_identity_asv** e **sintax_cutoff** são utilizados apenas para a abordagem com ASVs.

## Como usar o _pipeline_

Para executar o _pipeline_, primeiramente devem ser configurados os parâmetros do arquivo [config.txt](#arquivo-de-configuração-config.txt)

Para ver a forma de usar (_help_)

```sh
  $ python3 amplicon_pipeline.py --help
  usage: amplicon_pipeline.py [-h] -c FILE [--version]

  Pipeline for analysis of 16s rRNA amplicons, using ASVs (Amplicon Sequence Variant) or OTUs (Operational Taxonomic Unit)

  optional arguments:
    -h, --help            show this help message and exit
    -c FILE, --config_file FILE
                        Configuration file
    --version             show program's version number and exit

  Thank you!
```

Forma de executar:

```sh
  python3 amplicon_pipeline.py -c config.txt
```

## Credits

- O _pipeline_ com a abordagem de OTUs foi baseado no _pipeline_ do VSEARCH proposto [aqui](https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline).
- O _pipeline_ com a abordagem de ASVs foi baseado no _pipeline_ proposto [aqui](https://astrobiomike.github.io/amplicon/workflow_ex).

## Author

* [Glen Jasper](https://github.com/glenjasper)

## Organization
* [Molecular and Computational Biology of Fungi Laboratory](https://sites.icb.ufmg.br/lbmcf/index.html) (LBMCF, ICB - **UFMG**, Belo Horizonte, Brazil).

## License

This project is licensed under the MIT License - see the [LICENSE](./LICENSE) file for details.

## Mandatory Attribution

**Any use, distribution, or modification of this software must include the following attribution:**

> This software was developed by Glen Jasper (https://github.com/glenjasper),  
> originally available at https://github.com/glenjasper/pipeline-for-amplicon-analysis
