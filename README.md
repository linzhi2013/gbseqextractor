# gbseqextractor

## 1 Introduction

`gbseqextractor` is a tool to extract any CDS or rNRA or tRNA DNA sequences of genes from Genbank file. with `Biopython` (http://www.biopython.org/)

## 2 Installation

    pip install gbseqextractor

There will be a command `gbseqextractor` created under the same directory as your `pip` command.

## 3 Usage

    $ gbseqextractor
    usage: gbseqextractor.py [-h] -f <STR> -prefix <STR> [-seqPrefix <STR>]
                             [-types {CDS,rRNA,tRNA,wholeseq} [{CDS,rRNA,tRNA,wholeseq} ...]] [-gi] [-p] [-t] [-s] [-l] [-rv]
                             [-F]

    Extract any CDS or rNRA or tRNA DNA sequences of genes from Genbank file.

    Seqid will be the value of '/gene=' or '/product=', if they both were not
    present, the gene will not be output!

    Please cite:
    Guanliang Meng, Yiyuan Li, Chentao Yang, Shanlin Liu,
    MitoZ: a toolkit for animal mitochondrial genome assembly, annotation
    and visualization, Nucleic Acids Research, https://doi.org/10.1093/nar/gkz173



    optional arguments:
      -h, --help            show this help message and exit
      -f <STR>              Genbank file
      -prefix <STR>         prefix of output file.
      -seqPrefix <STR>      prefix of each seq id. default: None
      -types {CDS,rRNA,tRNA,wholeseq} [{CDS,rRNA,tRNA,wholeseq} ...]
                            what kind of genes you want to extract? wholeseq for whole fasta seq.[CDS]
      -gi                   use gi number as sequence ID instead of accession number when " gi number is present. (default:
                            accession number)
      -p                    output the position information on the ID line. Warning: the position on ID line is 0 left-most!
                            [False]
      -t                    output the taxonomy lineage on ID line [False]
      -s                    output the species name on the ID line [False]
      -l                    output the seq length on the ID line [False]
      -rv                   reverse and complement the sequences if the gene is on minus strand [False]
      -F                    only output full length genes,i.e., exclude the genes with '>' or '<' in their location [False]

## Author
Guanliang MENG

## Citation
This script is part of the package `MitoZ`, when you use the script in your work, please cite:

    Guanliang Meng, Yiyuan Li, Chentao Yang, Shanlin Liu,
    MitoZ: a toolkit for animal mitochondrial genome assembly, annotation and visualization, Nucleic Acids Research, https://doi.org/10.1093/nar/gkz173


Meanwhile, since `gbseqextractor` makes use of `Biopython`, you should alos cite it if you use `gbseqextractor` in your work:

    Peter J. A. Cock, Tiago Antao, Jeffrey T. Chang, Brad A. Chapman, Cymon J. Cox, Andrew Dalke, Iddo Friedberg, Thomas Hamelryck, Frank Kauff, Bartek Wilczynski, Michiel J. L. de Hoon: “Biopython: freely available Python tools for computational molecular biology and bioinformatics”. Bioinformatics 25 (11), 1422–1423 (2009). https://doi.org/10.1093/bioinformatics/btp163

Please go to `http://www.biopython.org/` for more details.






