# rqt: utilities for gene-level meta-analysis

## Installation

### Development version

```rqt``` is currently under consideration for acceptance into Bioconductor:  https://github.com/Bioconductor/Contributions/issues/212
and hence requires the development version of R (>=3.4) and the development version of Bioconductor (3.5).
If you have these installed, then ```rqt``` can be installed from Github using biocLite:

```
source("https://bioconductor.org/biocLite.R")
biocLite("izhbannikov/rqt", dependencies=TRUE, build_vignettes=TRUE)
```

### Release version

The last version of rqt that is compatible with the current version of R (3.3), 
which can be downloaded using devtools:

```
devtools::install_github("izhbannikov/rqt@develop", buildVignette=TRUE)
```

There are no significant changes to functionality (except R >=3.3).


## Usage

###Single dataset

```
library(rqt)
# Loading data and constructing the objects #
data <- data.matrix(read.table(system.file("extdata/test.bin1.dat",
    package="rqt"), header=TRUE))
pheno <- data[,1]
geno <- data[, 2:dim(data)[2]]
colnames(geno) <- paste(seq(1, dim(geno)[2]))
geno.obj <- SummarizedExperiment(geno)
obj <- rqt(phenotype=pheno, genotype=geno.obj)
# Analysis #
res <- geneTest(obj, method="pca", out.type = "D")
print(res)
```

### Multiple datasets (meta analysis)
```
library(rqt)
data1 <- data.matrix(read.table(system.file("extdata/phengen2.dat",
                                            package="rqt"), skip=1))
pheno <- data1[,1]
geno <- data1[, 2:dim(data1)[2]]
colnames(geno) <- paste(seq(1, dim(geno)[2]))
geno.obj <- SummarizedExperiment(geno)
obj1 <- rqt(phenotype=pheno, genotype=geno.obj)

data2 <- data.matrix(read.table(system.file("extdata/phengen3.dat",
                                            package="rqt"), skip=1))
pheno <- data2[,1]
geno <- data2[, 2:dim(data2)[2]]
colnames(geno) <- paste(seq(1, dim(geno)[2]))
geno.obj <- SummarizedExperiment(geno)
obj2 <- rqt(phenotype=pheno, genotype=geno.obj)

data3 <- data.matrix(read.table(system.file("extdata/phengen.dat",
                                            package="rqt"), skip=1))
pheno <- data3[,1]
geno <- data3[, 2:dim(data3)[2]]
colnames(geno) <- paste(seq(1, dim(geno)[2]))
geno.obj <- SummarizedExperiment(geno)
obj3 <- rqt(phenotype=pheno, genotype=geno.obj)

res.meta <- geneTestMeta(list(obj1, obj2, obj3))
print(res.meta)
```
