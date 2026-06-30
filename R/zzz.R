# Quiet the R CMD check note for the ggplot2 .data pronoun used in plot.maki_test.
if (getRversion() >= "2.15.1") utils::globalVariables(".data")
