Curso de treinamento para _amplicon_ 16S rRNA
======================
[![License](https://poser.pugx.org/badges/poser/license.svg)](./LICENSE)

Curso de treinamento para análise de _amplicon_ 16S rRNA utilizando as abordabens de agrupamento por OTUs (Operational Taxonomic Unit) e geração de ASVs (Amplicon Sequence Variant). Se disponibilizam scripts em [_Python 3_](#pipeline-em-python-3) e [_Shell script_](#pipeline-em-shell-script), para ambas as abordagens.

## Recursos

- [Dados de treinamento](#arquivos-fastq)
- [Banco de dados SILVA](#banco-de-dados-silva)
    - [Formatação do banco de dados SILVA](#formatação-do-banco-de-dados-silva)
- [Pré-requisitos](#pré-requisitos)
    - [Programas](#programas)
    - [Bibliotecas Python](#bibliotecas-python)
    - [Bibliotecas R](#bibliotecas-r)
- [Scripts](#scripts)
- [Pipelines](#pipelines)
    - [Pipeline para clusterização de OTUs](#pipeline-para-clusterização-de-otus)
    - [Pipeline para geração de ASVs](#pipeline-para-geração-de-asvs)
- [Author](#author)
- [Organization](#organization)
- [License](#license)

## Arquivos-FASTQ
**training-files.zip**: Arquivos de treinamento. Usaremos 4 tipo de amostras (BRS, BPA, BANHT e ENV), pode baixar os dados desde [aqui](https://drive.google.com/file/d/1cvn8NVWhU0C5dbOj9gWKsPrt9G58kbfR/view?usp=sharing).

## Banco de dados SILVA
Usaremos o banco de dados [SILVA 138 SSU NR](https://www.arb-silva.de/no_cache/download/archive/current/Exports). Pode-se baixar e criar os binários com o seguinte comando:

```sh
  wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
  gzip -d SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz

  makeblastdb -in SILVA_138.1_SSURef_NR99_tax_silva.fasta \
              -dbtype nucl \
              -out silva_db
```

**NOTA**: Para o treinamento podemos usar uma versão _small_ do banco de dados. Pode baixá-lo desde [aqui](https://drive.google.com/file/d/1N7PL1VVn6xOyYA2DVkwvs2Nkxh2Dio13/view?usp=sharing) ou criá-lo com o seguinte comando: 

```sh
  sed -n 1,1500032p SILVA_138.1_SSURef_NR99_tax_silva.fasta > SILVA_138.1_SSURef_NR99_tax_silva_small.fasta

  makeblastdb -in SILVA_138.1_SSURef_NR99_tax_silva_small.fasta \
              -dbtype nucl \
              -out silva_db_small
```

### Formatação do banco de dados SILVA
O banco de dados SILVA não pode ser utilizado diretamente no _pipeline_ com ASVs, entretanto pode ser formatado com o script _rename_silva.py_, da seguinte forma. Se estiver utilizando a versão _small_ do banco de dados, também deve ser formatado.

```sh
  python3 rename_silva.py SILVA_138.1_SSURef_NR99_tax_silva.fasta
```

## Pré-requisitos
### Programas
Os _pipelines_ precisam dos seguintes programas (ou linguagen de programação). Os pipelines em 
- [VSEARCH](https://github.com/torognes/vsearch)
- [USEARCH 32-bit](https://drive5.com/usearch)
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)
- [Python 3](https://www.python.org)


- [R-4](https://cran.r-project.org)
- [RStudio](https://www.rstudio.com/products/rstudio/download)
- 
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
- **map.py**: Script para mapear leituras non-singletons e non-chimeras (adaptado de [map.pl](https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline)).
- **rename_silva.py**: Script para renomear os cabeçalhos do arquivo FASTA do banco de dados [SILVA 138 SSU NR](https://www.arb-silva.de/no_cache/download/archive/current/Exports) para ser utilizado com o _pipeline_ para a geração de ASVs.
- **reverse_complement.py**: Script para obter a reversa-complementar de um _primer_.
- **get_abundances_table_otu.py**: Script para obter a tabela de abundâncias dos OTUs com dados taxonômicos.
- **get_abundances_table_asv.py**: Script para obter a tabela de abundâncias dos ASVs com dados taxonômicos.
- **get_abundances_by_tax.py**: Script para obter tabelas de dados para gerar Heatmaps, Diagramas de Venn e Bar-plots.

## _Pipelines_

### _Pipeline_ em _Python 3_

#### _Pipeline_ para clusterização de OTUs
- **amplicon_pipeline_otu.sh**: Fluxo (_pipeline_) para a geração de uma tabela de abundâncias de OTUs com dados taxonômicos utilizando o banco de dados SILVA, a partir de dados de sequenciamentos de _amplicon_ 16S rRNA.

#### _Pipeline_ para geração de ASVs
- **amplicon_pipeline_asv.sh**: Fluxo (_pipeline_) para a geração de uma tabela de abundâncias de ASVs com dados taxonômicos utilizando o banco de dados SILVA, a partir de dados de sequenciamentos de _amplicon_ 16S rRNA.

### _Pipeline_ em _Shell Script_

## Author

* [Glen Jasper](https://github.com/glenjasper)

## Organization
* [Molecular and Computational Biology of Fungi Laboratory](https://sites.icb.ufmg.br/lbmcf/index.html) (LBMCF, ICB - **UFMG**, Belo Horizonte, Brazil).

## License

This project is licensed under the MIT License - see the [LICENSE](./LICENSE) file for details.
