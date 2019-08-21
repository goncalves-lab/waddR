library("testthat")
library("waddR")
library("SingleCellExperiment")

ts.names <- c(  "d.wass", "d.wass^2", "d.comp^2", "d.comp", "location",
                "size", "shape", "rho", "pval", "p.ad.gdp", "N.exc",
                "perc.loc", "perc.size", "perc.shape", "decomp.error",
                "p.zero", "p.combined", "p.adj.nonzero","p.adj.zero",
                "p.adj.combined")

os.names <-  c( "d.wass", "d.wass^2", "d.comp^2", "d.comp", "location",
                "size", "shape", "rho", "pval", "p.ad.gdp", "N.exc",
                "perc.loc", "perc.size", "perc.shape", "decomp.error",
                "pval.adj")

# input data: x, y represent the expression levels for 16 individuals in two
# conditions of one gene
x <- c( 0.0000000, 1.3389541, 0.9190596, 1.5532290, 0.0000000, 1.1014202,
        1.1996965, 1.6332709, 1.0306684, 0.8707476, 1.8230579, 0.0000000,
        1.7195922, 1.1879964, 1.2522275, 0.9516880, 1.0937032, 0.0000000,
        1.3789891, 1.6415560, 0.9871225, 0.0000000, 0.0000000, 0.0000000,
        1.2619454, 1.9289606, 1.7133513, 0.6006942, 0.0000000, 0.7728401,
        0.0000000, 1.2871939, 1.1820060, 1.2520491, 1.6542059, 1.7451944,
        1.3761065, 1.7967770, 0.0000000, 1.1881549, 0.0000000, 0.0000000,
        0.0000000, 2.5839898, 1.0142013, 0.0000000, 0.0000000, 0.0000000,
        1.0225624, 1.3279117, 2.0216712, 1.4012923, 1.7473436, 2.1056459,
        1.5037761, 0.7913705, 1.7510320, 1.0444939, 1.0175730, 1.4748563,
        1.8437568, 1.0464078, 1.6420597, 0.0000000, 1.2362845, 1.3117466,
        1.6209499, 2.8094403)

y <- c( 0.0000000, 2.0119764, 0.0000000, 1.9831597, 2.2152026, 0.0000000,
        0.0000000, 0.0000000, 1.1296844, 1.8360126, 0.0000000, 1.4104730,
        0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
        1.1391284, 0.0000000, 2.1396251, 1.2516926, 0.0000000, 2.2685640,
        1.4992348, 0.0000000, 1.4731691, 0.0000000, 0.0000000, 1.0289898,
        2.3849062, 0.0000000, 0.9660725, 2.2188470, 1.2731303, 2.3051090,
        1.4992348, 0.0000000, 1.4731691, 0.0000000, 1.1537754, 2.0288085,
        2.1383034, 0.0000000, 1.5264164, 0.0000000 )

z <- rep(0, 40)
z2 <- rep(0,34)

# matrix input
dat <- matrix(c(x, y), nrow=1)
dat2 <- matrix(c(x, z), nrow=1)
dat3 <- matrix(c(z, z2), nrow=1)
condition1 <- c(rep(0, length(x)), rep(1, length(y)))
condition2 <- c(rep(0, length(x)), rep(1, length(z)))
condition3 <- c(rep(0, length(z)), rep(1, length(z2)))

# SingleCellExperiment input
sce.a <- SingleCellExperiment(assays=list(counts=matrix(x, nrow=1)))
sce.b <- SingleCellExperiment(assays=list(counts=matrix(y, nrow=1)))
sce.a2 <- SingleCellExperiment(assays=list(counts=matrix(z, nrow=1)))
sce.b2 <- SingleCellExperiment(assays=list(counts=matrix(z2, nrow=1)))

test_that("Correctness of wasserstein single cell output", {

    # "normal" expression values -- Two Stage test
    res <- wasserstein.sc(dat, condition1, 10, "TS")
    res.dup <- wasserstein.sc(sce.a, sce.b, 10, "TS")

    # reference results
    res.values <- matrix(c(0.3237935, 0.1048422, 0.1060877, 0.3257111,
                           0.07527902, 0.0007892882, 0.0300194, 0.9262569,
                           0.09090909, NA, NA, 70.96, 0.74, 28.3,
                           0.01173997, 0.9434837, 0.2964316, 0.09090909,
                           0.9434837, 0.2964316), nrow=1)
    colnames(res.values) <- ts.names

    # compare
    expect_equal(res.dup, res.values, tolerance=0.00001)
    expect_equal(res, res.values, tolerance=0.00001)
    expect_equal(res, res.dup)

    # "normal" expression values -- One Stage test
    res.os1 <- wasserstein.sc(dat, condition1, 10, "OS")
    res.os.dup1 <- wasserstein.sc(sce.a, sce.b, 10, "OS")

    # reference results
    res.os.values1 <- matrix(c(0.4926601, 0.242714, 0.2474003, 0.4973935,
                               0.03169485, 0.03660025, 0.1791052, 0.8637736,
                               0.09090909, NA, NA, 12.81, 14.79, 72.39,
                               0.01894254, 0.09090909), nrow=1)
    colnames(res.os.values1) <- os.names

    expect_equal(res.os1, res.os.values1, tolerance=0.0000001)
    expect_equal(res.os.dup1, res.os.values1, tolerance=0.0000001)
    expect_equal(res.os.dup1, res.os1)

    # case of a zero-expressed gene -- Two Stage
    res.ts <- wasserstein.sc(dat3, condition3, 10, "TS")
    res.ts.dup <- wasserstein.sc(sce.a2, sce.b2, 10, "TS")
    res.ts.values2 <- matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                               NA, NA, NA, 1, 1, NA, 1, 1), nrow=1)
    colnames(res.ts.values2) <- ts.names
    expect_equal(res.ts, res.ts.values2)
    expect_equal(res.ts.dup, res.ts.values2)
    expect_equal(res.ts, res.ts.dup)

    # case of a zero-expressed gene -- One Stage
    res.os2 <- wasserstein.sc(dat3, condition3, 10, "OS")
    res.os.dup2 <- wasserstein.sc(sce.a2, sce.b2, 10, "OS")
    res.os.values2 <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 1, NA, NA, NaN, NaN,
                               NaN, NaN, 1), nrow=1)
    colnames(res.os.values2) <- os.names
    expect_equal(res.os2, res.os.values2)
    expect_equal(res.os.dup2, res.os.values2)
    expect_equal(res.os2, res.os.dup2)
})


test_that("Example run wasserstein single cell", {

    expect_equal(   colnames(wasserstein.sc(dat, condition1, 10, "TS")), 
                    ts.names)
    expect_equal(   colnames(wasserstein.sc(dat2, condition2, 10, "OS")),
                    os.names)
    expect_equal(   colnames(wasserstein.sc(sce.a, sce.b, 10, "OS")),
                    os.names)
    expect_equal(   colnames(wasserstein.sc(sce.a, sce.b2, 10, "TS")),
                    ts.names)
})
