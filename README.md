Single-cell RNA analysis of immune cells from melanoma tumors
=============================================================

This repository holds the code that reproduce the analysis done in the
[Li et al,
2018](https://www.cell.com/cell/fulltext/S0092-8674(18)31568-X) paper.
It also downloads the required processed and auxilary data.

The core analysis is done with the
[metacell](https://tanaylab.github.io/metacell) R
package.

#### Quick Links

-   Metacell paper: Baran et al. 2018
    ([bioarxiv](https://www.biorxiv.org/content/early/2018/10/08/437665)).
-   Metacell R package
    [homepage](https://tanaylab.bitbucket.io/metacell-r/index.html)
    (with functions reference and usage vignettes).
-   Raw data is available under EGA accession EGAS00001003363 (access is
    restricted, follow guidelines at the EGA site)
-   Processed data is avaliable under GEO accession
    [GSE123139](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123139)
-   Raw UMI counts for all cells with metacell membership information:
    [Li2018\_umi\_counts.tar.gz](https://li-melanoma-2018.s3.eu-west-1.amazonaws.com/Li2018_umi_counts.tar.gz).
    This gzipped archive contains a gene-cells UMI table in
    matrix-market format and a metadata table per cell. The metadata
    table also contains the metacell membership and the annotated group
    per metacell model - all cells, T/NK cells and non-T/NK cells.
-   Code repository for the TCR sequences analysis is available at
    [TCRseq](https://github.com/DiklaGelbard/TCRseq) github page.

#### Requirements

R with these packages:

-   Matrix
-   pheatmap
-   flowCore
-   dplyr
-   glmnet
-   plotrix
-   MASS
-   data.table
-   [metacell](https://tanaylab.github.io/metacell/index.html)

**Note**: Metacell is implemented in R and C++. In particular it uses
the Tanay group tgstat library that utilizes shared memory and
distributed computing (as well as some specific optional CPU features).
The package is tested on linux and macbooks, and is currently not
compatible on Windows. A typical application of metacell requires at
least 16G RAM. For the current dataset we recommend a dual CPU
multi-core workstation with 128GM RAM or more.

#### Usage

In an R session opened on the repository root directory

``` r
# Loading code and downloading required data files
source("pipe.r")

# Building the metacells
build_metacells()

# Generate figures
generate_figs()

# Reproduce the Guo et al 2018 (lung cancer scRNA data) analysis
source("guo2018.r")
```

**Note**: Metacell generation was sensitive to the initial random seed,
so the exact metacell solution cannot be exactly reproduced.
build\_metacells function is included for reference, to see how the
metacells were produced. In order to reproduce the exact figures as in
the paper (e.g. supporting the manually selected metacell IDs for
annotation), the metacell object used in the paper is supplied and used
by the generate\_figs function. This function also does a cross-metacell
comparison to compare the cell membership between the 2 metacell
solutions, which demonstrate that the differences are not substantial.

#### Contact

For help, please contact <yaniv.lubling@weizmann.ac.il>
