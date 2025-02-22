library(SuperLearner)
library(AIPW)
load("out/NYC.clustered.Rdata")
source("utils.R")

set.seed(123)
rtmle_res <- rtmle(W = W,
                   A = A.outlet,
                   Y = Y.yr,
                   a = 0,
                   b = 0.085,
                   wt = rep(1, nrow(W)))
tmle_res <- rtmle(W = W,
                  A = A.outlet,
                  Y = Y.yr,
                  a = 0,
                  b = 1,
                  wt = rep(1, nrow(W)))
