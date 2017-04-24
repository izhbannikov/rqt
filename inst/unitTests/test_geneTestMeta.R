test_geneTestMeta <- function() {
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
    checkEqualsNumeric(res.meta$final.pvalue, 0.005502565, tolerance=1.0e-5)
}
