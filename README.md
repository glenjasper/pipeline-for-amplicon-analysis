Curso de treinamento para _amplicon_ 16S rRNA
======================
[![License](https://poser.pugx.org/badges/poser/license.svg)](./LICENSE)

Curso para treinamento em _amplicon_ 16S rRNA

## Recursos

- [Dados de treinamento](#arquivos-fastq)
- [Pré-requisitos](#pré-requisitos)
    - [Programas](#programas)
    - [Bibliotecas Python](#bibliotecas-python)
    - [Bibliotecas R](#bibliotecas-r)
- [Scripts](#scripts)
- [Author](#author)
- [Organization](#organization)
- [License](#license)

## Arquivos-FASTQ
**training-files.zip**: Arquivos de treinamento. Usaremos 4 tipo de amostras (BRS, BPA, BANHT e ENV), pode baixar os dados desde [aqui](https://drive.google.com/file/d/1cvn8NVWhU0C5dbOj9gWKsPrt9G58kbfR/view?usp=sharing).

## Pré-requisitos
### Programas
Deve ter instalado os seguintes programas:
- [VSEARCH](https://github.com/torognes/vsearch)
- [USEARCH 32-bit](https://drive5.com/usearch)
- [R-4](https://cran.r-project.org)
- [RStudio](https://www.rstudio.com/products/rstudio/download)
- [Python 3](https://www.python.org)
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

### Bibliotecas Python

```sh
  $ sudo pip3 install pandas
  $ sudo pip3 install biopython
  $ sudo pip3 install colorama
  $ sudo pip3 install cutadapt
```

### Bibliotecas R

```sh
  install.packages('ggplot2')
  install.packages('heatmaply')
  install.packages('dendextend')
  install.packages('RColorBrewer')
  install.packages('viridis')
```

## Scripts
- **map.py**: Script para mapear leituras non-singletons e non-chimeras (adaptado de [map.pl](https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline)).
- **rename_silva.py**: Script para renomear os cabeçalhos do arquivo FASTA do bando de dados [SILVA 138 SSU NR](https://www.arb-silva.de/no_cache/download/archive/current/Exports) para ser utilizado com o _pipeline_ para a geração de ASVs.
- **reverse_complement.py**: Script para obter a reversa-complementar de um _primer_.
- **get_abundances_table_otu.py**: Script para obter a tabela de adundâncias dos OTUs com dados taxonômicos.
- **get_abundances_table_asv.py**: Script para obter a tabela de adundâncias dos ASVs com dados taxonômicos.
- **get_abundances_by_tax.py**: Script para obter tabelas de dados para gerar Heatmaps, Diagramas de Venn e Bar-plots.

## Author

* [Glen Jasper](https://github.com/glenjasper)

## Organization
* [Molecular and Computational Biology of Fungi Laboratory](https://sites.icb.ufmg.br/lbmcf/index.html) (LBMCF, ICB - **UFMG**, Belo Horizonte, Brazil).

## License

This project is licensed under the MIT License - see the [LICENSE](./LICENSE) file for details.
