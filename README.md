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
    - [Pipeline em Shell Script](#pipeline-em-shell-script)
- [Author](#author)
- [Organization](#organization)
- [License](#license)

## Dados de treinamento

**training-files.zip**: Arquivos _pair-end_ FASTQ de treinamento. Neste treinamento utilizaremos dados de 4 tipo de amostras com suas réplicas biológicas que estão representados em 11 arquivos FASTQ pareados (3xBRS, 3xBPA, 3xBANHT e 2xENV), podem ser baixados desde [aqui](https://drive.google.com/file/d/1cvn8NVWhU0C5dbOj9gWKsPrt9G58kbfR/view?usp=sharing).

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

- **util/map.py**: Script para mapear leituras non-singletons e non-chimeras (adaptado de [map.pl](https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline)).
- **util/rename_silva.py**: Script para renomear os cabeçalhos do arquivo FASTA do banco de dados [SILVA 138 SSU NR](https://www.arb-silva.de/no_cache/download/archive/current/Exports) para ser utilizado com o _pipeline_ para a geração de ASVs.
- **util/reverse_complement.py**: Script para obter a reversa-complementar de um _primer_.
- **util/get_abundances_table_otu.py**: Script para obter a tabela de abundâncias dos OTUs com dados taxonômicos.
- **util/get_abundances_table_asv.py**: Script para obter a tabela de abundâncias dos ASVs com dados taxonômicos.
- **util/get_abundances_by_tax.py**: Script para obter tabelas de dados para gerar Heatmaps, Diagramas de Venn e Bar-plots.

## _Pipelines_

### _Pipeline_ em _Python 3_

Se desenvolveram dois _scripts_ em _Python 3_ para ambas as unidades de medidas (OTUs e ASVs). Estes _scripts_ podem ser utilizados tanto na plataforma **GNU/Linux** (ou Mac OS) quanto no **Windows**. Estes _script_ já têm incorporado os [programas](#programas) requeridos na pasta **pipeline-python/utilities** (arquivos: _common.zip_, _gnulinux.zip_ e _win.zip_), os quais devem ser descompactados pelo usuário. O usuário também precisa configurar os parâmetros do arquivo **config.txt**.

#### _Pipeline_ com abordagem de OTUs
- **pipeline-python/amplicon_pipeline_for_otu.py**: Fluxo (_pipeline_) para a geração de uma tabela de abundâncias de OTUs com dados taxonômicos utilizando o banco de dados SILVA, a partir de dados de sequenciamentos de _amplicon_ 16S rRNA.

#### _Pipeline_ com abordagem de ASVs
- **pipeline-python/amplicon_pipeline_for_asv.py**: Fluxo (_pipeline_) para a geração de uma tabela de abundâncias de ASVs com dados taxonômicos utilizando o banco de dados SILVA, a partir de dados de sequenciamentos de _amplicon_ 16S rRNA.

#### Exemplo de configuração do arquivo **config.txt**:

```sh
  [PARAMETERS]
  # Paths
  samples_path = D:\my_workstation\amplicon\dataset
  database_path = D:\my_workstation\amplicon\database
  util_path = D:\my_workstation\amplicon\util
  output_path = D:\my_workstation\amplicon\results-otu

  # Taxonomy database (files must be in database_path)
  database_fasta = SILVA_138.1_SSURef_NR99_tax_silva.fasta
  database_bin = silva_db

  # Primers file (file must be in database_path)
  primers_file = illumina.primers.fa

  # Multiprocessing
  threads = 4

  # OS type (gnulinux: for GNU/Linux | win: for Windows)
  os_type = win

  # Python version (python3: for Python 3.x in GNU/Linux | python: for Python 3.x in Windows)
  python_version = python
```

> - samples_path
> - database_path
> - util_path
> - output_path
> - database_fasta
> - database_bin
> - primers_file
> - threads
> - os_type
> - python_version

| Syntax      | Description |
| ----------- | ----------- |
| Header      | Title       |
| Paragraph   | Text        |

> **Nota**: Para a abordagem com ASV não se precisa do parâmetro **database_bin**.

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
  THREADS=4
```

> **Nota**: Para a abordagem com ASV não se precisa do parâmetro **database_bin**.

## Author

* [Glen Jasper](https://github.com/glenjasper)

## Organization
* [Molecular and Computational Biology of Fungi Laboratory](https://sites.icb.ufmg.br/lbmcf/index.html) (LBMCF, ICB - **UFMG**, Belo Horizonte, Brazil).

## License

This project is licensed under the MIT License - see the [LICENSE](./LICENSE) file for details.
