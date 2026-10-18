#!/usr/bin/env Rscript
# Cross-language metrics parity fixture, deterministic numbers exported by Python.
source('R/functions/metrics.R')
a <- read.csv('data/synthetic/metric_fixture.csv')
ref <- read.csv('data/synthetic/metric_fixture_expected.csv')
stopifnot(abs(rmse(a$observed,a$predicted)-ref$rmse)<1e-10)
stopifnot(abs(mae(a$observed,a$predicted)-ref$mae)<1e-10)
cat('Fixture metric parity: R matches Python reference.\n')
