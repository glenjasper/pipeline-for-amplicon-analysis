Curso de treinamento para _amplicon_ 16S rRNA
======================
[![License](https://poser.pugx.org/badges/poser/license.svg)](./LICENSE)

Curso de treinamento para análise de _amplicon_ 16S rRNA utilizando as abordabens de agrupamento por OTUs (Operational Taxonomic Unit) e geração de ASVs (Amplicon Sequence Variant). Para ambas as abordagens, se disponibilizam scripts em [_Python 3_](#pipeline-em-python-3) que podem ser executados nas plataformas **GNU/Linux** e **Windows** e [_Shell script_](#pipeline-em-shell-script) que pode ser executado apenas em **GNU/Linux** (ou Mac OS).

## Recursos

- [Dados de treinamento](#dados-de-treinamento)
- [Banco de dados SILVA](#banco-de-dados-silva)
    - [Baixar e criar os binários do banco de dados](#baixar-e-criar-os-binários-do-banco-de-dados)
    - [Formatação do banco de dados para uso de ASVs](#formatação-do-banco-de-dados-para-uso-de-asvs)
- [Pré-requisitos](#pré-requisitos)
    - [Programas](#programas)
    - [Bibliotecas Python](#bibliotecas-python)
    - [Bibliotecas R](#bibliotecas-r)
- [Scripts](#scripts)
- [Pipelines](#pipelines)
    - [Pipeline em Python 3](#pipeline-em-python-3)
    - [Arquivo de configuração config.txt](#arquivo-de-configuração-config-txt)
    - [Pipeline em Shell Script](#pipeline-em-shell-script)
- [Como usar o pipeline](#como-usar-o-pipeline)
- [Credits](#credits)
- [Author](#author)
- [Organization](#organization)
- [License](#license)

## Dados de treinamento

**training-files.zip**: Arquivos _pair-end_ FASTQ de treinamento. Neste treinamento utilizaremos dados de 4 tipo de amostras com suas réplicas biológicas que estão representados em 11 arquivos FASTQ pareados (3 x BRS, 3 x BPA, 3 x BANHT e 2 x ENV), podem ser baixados desde [aqui](https://drive.google.com/file/d/1cvn8NVWhU0C5dbOj9gWKsPrt9G58kbfR/view?usp=sharing).

## Banco de dados SILVA

### Baixar e criar os binários do banco de dados
Usaremos o banco de dados [SILVA 138 SSU NR](https://www.arb-silva.de/no_cache/download/archive/current/Exports). Pode baixar o arquivo FASTA e criar os binários com os seguintes comandos:

```sh
  wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
  gzip -d SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz

  makeblastdb -in SILVA_138.1_SSURef_NR99_tax_silva.fasta \
              -dbtype nucl \
              -out silva_db
```

> **NOTA**: Para o treinamento podemos usar uma versão _small_ do banco de dados. Pode baixá-lo desde [aqui](https://drive.google.com/file/d/1N7PL1VVn6xOyYA2DVkwvs2Nkxh2Dio13/view?usp=sharing) ou criá-lo com o seguintes comandos: 

```sh
  sed -n 1,1500032p SILVA_138.1_SSURef_NR99_tax_silva.fasta > SILVA_138.1_SSURef_NR99_tax_silva_small.fasta

  makeblastdb -in SILVA_138.1_SSURef_NR99_tax_silva_small.fasta \
              -dbtype nucl \
              -out silva_db_small
```

### Formatação do banco de dados para uso de ASVs
O banco de dados SILVA não pode ser utilizado diretamente no _pipeline_ com ASVs, entretanto pode ser formatado com o script _rename_silva.py_, da seguinte forma. Se estiver utilizando a versão _small_ do banco de dados, também deve ser formatado.

```sh
  python3 rename_silva.py SILVA_138.1_SSURef_NR99_tax_silva.fasta
```

## Pré-requisitos

### Programas
Os _pipelines_ precisam dos seguintes programas (ou linguagen de programação). Para o [_pipeline_ em _Shell Script_](#pipeline-em-shell-script) estes requerimentos necessáriamente tem que estar instalados, no entanto, se utilizar o [_pipeline_ em _Python 3_](#pipeline-em-python-3), estes requerimentos já estão inclusos no código fonte (pasta **pipeline-python/utilities**).
- [VSEARCH](https://github.com/torognes/vsearch)
- [USEARCH 32-bit](https://drive5.com/usearch)
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)
- [Python 3](https://www.python.org)

Para as análises estatísticas se precisam ter instalados:
- [R-4](https://cran.r-project.org)
- [RStudio](https://www.rstudio.com/products/rstudio/download)

### Bibliotecas Python

```sh
  sudo pip3 install pandas
  sudo pip3 install biopython
  sudo pip3 install colorama
```

### Bibliotecas R

```sh
  install.packages('ggplot2')
  install.packages('heatmaply')
  install.packages('dendextend')
  install.packages('RColorBrewer')
  install.packages('viridis')
```

## _Scripts_

- **util/map.py**: _Script_ para mapear leituras _non-singletons_ e _non-chimeras_ (adaptado de [map.pl](https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline)).
- **util/rename_silva.py**: _Script_ para renomear os cabeçalhos do arquivo FASTA do banco de dados [SILVA 138 SSU NR](https://www.arb-silva.de/no_cache/download/archive/current/Exports) para ser utilizado com o _pipeline_ com ASVs.
- **util/reverse_complement.py**: _Script_ para obter a reversa-complementar de uma sequência (_forward-primer_).
- **util/get_abundances_table_otu.py**: _Script_ para obter a tabela de abundâncias dos OTUs com dados taxonômicos.
- **util/get_abundances_table_asv.py**: _Script_ para obter a tabela de abundâncias dos ASVs com dados taxonômicos.
- **util/get_abundances_by_tax.py**: _Script_ para obter tabelas de dados para gerar _Heatmaps_, Diagramas de _Venn_ e Bar-plots.

## _Pipelines_

### _Pipeline_ em _Python 3_

Se desenvolveu um _script_ em _Python 3_ para ambas as unidades de medidas (OTUs e ASVs). Este _script_ podem ser utilizados tanto na plataforma **GNU/Linux** (ou Mac OS) quanto no **Windows**. O _script_ já têm incorporado os [programas](#programas) requeridos na pasta **pipeline-python/utilities** (arquivos: _common.zip_, _gnulinux.zip_ e _win.zip_), os quais devem ser descompactados pelo usuário. O usuário também precisa configurar os parâmetros do arquivo **config.txt**.

#### _Pipeline_ para as abordagens de ASVs e OTUs
- **pipeline-python/amplicon_pipeline.py**: Fluxo (_pipeline_) para a geração de uma tabela de abundâncias de ASVs ou OTUs com dados taxonômicos utilizando o banco de dados SILVA, a partir de dados de sequenciamentos de _amplicon_ 16S rRNA.

#### Arquivo de configuração config.txt:

Exemplo de configuração do arquivo **config.txt**:

```sh
  [PARAMETERS]

  # Type of approach to analysis (ASV or OTU)
  approach_type = asv

  # Paths
  samples_path = D:\my_workstation\amplicon\dataset
  database_path = D:\my_workstation\amplicon\database
  util_path = D:\my_workstation\amplicon\util
  output_path = D:\my_workstation\amplicon\results-otu

  # Taxonomy database (files must be in database_path | database_bin only for OTUs)
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
  filter_minlen = 
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
| **util_path**         | Caminho absoluto da pasta que contem os _scripts_ utilitários. |
| **output_path**       | Caminho absoluto da pasta de saída. |
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

### _Pipeline_ em _Shell Script_

Se desenvolveram dois _scripts_ em _Shell Script_ para ambas as unidades de medidas (OTUs e ASVs). Estes _scripts_ podem ser utilizados apenas na plataforma **GNU/Linux** (ou Mac OS). Cada _script_ apresenta internamente parâmetros a serem configurados pelo usuário. 

#### _Pipeline_ com abordagem de OTUs
- **pipeline-shell/amplicon_pipeline_for_otu.sh**: Fluxo (_pipeline_) para a geração de uma tabela de abundâncias de OTUs com dados taxonômicos utilizando o banco de dados SILVA, a partir de dados de sequenciamentos de _amplicon_ 16S rRNA.

#### _Pipeline_ com abordagem de ASVs
- **pipeline-shell/amplicon_pipeline_for_asv.sh**: Fluxo (_pipeline_) para a geração de uma tabela de abundâncias de ASVs com dados taxonômicos utilizando o banco de dados SILVA, a partir de dados de sequenciamentos de _amplicon_ 16S rRNA.

#### Exemplo de configuração os parâmetros internos dos _scripts_:
```sh
  # Paths
  samples_path=D:\my_workstation\amplicon\dataset
  database_path=D:\my_workstation\amplicon\database
  util_path=D:\my_workstation\amplicon\util
  output_path=D:\my_workstation\amplicon\results-otu

  # Taxonomy database (files must be in database_path)
  database_fasta=SILVA_138.1_SSURef_NR99_tax_silva.fasta
  database_bin=silva_db

  # Primers file (file must be in database_path)
  primers_file=illumina.primers.fa

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

  # High identity to count ASVs (default: 99)
  high_identity_asv=99

  # For taxonomic assignment (default: 0.8)
  sintax_cutoff=0.8
```

> **Nota**: Os parâmetros **database_bin**, **cluster_identity** e **blast_identity** são utilizados apenas para a abordagem com OTUs. Os parâmetros **high_identity_asv** e **sintax_cutoff** são utilizados apenas para a abordagem com ASVs.

## Como usar o _pipeline_

Para executar qualquer _pipeline_, primeiramente devem ser configurados os parâmetros do arquivo [config.txt](#arquivo-de-configuração-config.txt)

- Executar o _pipeline_ em _Python 3_:

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

- Executar o _pipeline_ em _Shell Script_:

```sh
  ./amplicon_pipeline_for_asv.sh
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
