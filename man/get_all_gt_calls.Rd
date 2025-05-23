% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geno_read_base.R
\name{get_all_gt_calls}
\alias{get_all_gt_calls}
\title{Get all possible genotype calls given a unique set of alleles}
\usage{
get_all_gt_calls(alleles, ploidity, sep = "")
}
\arguments{
\item{alleles}{List[String]. A list of unique alleles, e.g., c("A", "B", "C").}

\item{ploidity}{Integer. The ploidity level of the organism.}
}
\value{
List[String]. A list of all possible genotype calls for the given alleles
        and ploidity level.
}
\description{
This function generates all possible genotype calls for a given set of alleles
and ploidity level. The genotype calls are represented as strings of characters,
with each allele being a single character.
}
\details{
The function uses a recursive approach to generate all possible combinations of
alleles for the specified ploidity level. For example, with two alleles "A" and "B",
and a ploidity of 2 (diploid), the function would generate the following genotype
calls: "AA", "AB", "BA", "BB".
}
\examples{
get_all_gt_calls(c("A", "B"), 2)  # Returns c("AA", "AB", "BA", "BB")
get_all_gt_calls(c("A", "B", "C"), 3)  # Returns all 27 possible triploid calls
get_all_gt_calls(c("A"), 1)  # Returns c("A")
}
