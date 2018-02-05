Readme
======

This repository contains data and scripts allowing to reproduce the figures of
[Universal rule for the symmetric division of plant cells](https://doi.org/10.1073/pnas.1011866108).

Prerequisites
-------------

*   MATLAB version R2007a or later
*   the API allowing to read the cell structures available at
    https://github.com/sbesson/plant-tissue/releases/tag/v0.1.0

Files
-----

Raw images used for performing the cell division analysis are available on request.

This repository contains

*   the results of the analysis of dividing cells for multiple systems saved as MAT files:

    | File | Cell type |
    |--------|----------------|
    | [coleochaeteresults.mat](coleochaeteresults.mat) | Coleochaete cells |
    | [quadrantcells.mat](quadrantcells.mat) | Dionaea, quadrant cells |
    | [dionaea_results.mat](dionaea_results.mat) | Dionaea, triangular cells |


*   the scripts used to produce the figures and movies of the publication

    | Script | Paper resource |
    |--------|----------------|
    | [Figure3A.m](Figure3A.m) | [Figure 3, panel A](http://www.pnas.org/content/108/15/6294#F3) |
    | [Figure3B.m](Figure3B.m) | [Figure 3, panel B](http://www.pnas.org/content/108/15/6294#F3) |
    | [Figure3B.m](Figure3C.m) | [Figure 3, panel C](http://www.pnas.org/content/108/15/6294#F3) |
    | [Figure6.m](Figure6.m) | [Figure 6](http://www.pnas.org/content/108/15/6294#F6) |
    | [MovieS1.m](MovieS3.m) | [Movie S1](http://www.pnas.org/content/108/15/6294/tab-figures-data) |
    | [MovieS3.m](MovieS3.m) | [Movie S3](http://www.pnas.org/content/108/15/6294/tab-figures-data) |
