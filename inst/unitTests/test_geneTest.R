test_geneTest <- function() {
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
    checkEqualsNumeric(results(res)$Qstatistic$Q1, 0.3798991, tolerance=1.0e-4)
    checkEqualsNumeric(results(res)$Qstatistic$Q2, 13.83897, tolerance=1.0e-4)
    checkEqualsNumeric(results(res)$Qstatistic$Q3, 0.3798991, tolerance=1.0e-4)
}
