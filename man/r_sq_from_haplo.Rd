% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_r_squared.R
\name{r_sq_from_haplo}
\alias{r_sq_from_haplo}
\title{Calculate the r-squared matrix from control haplotype data from Hapgen2.}
\usage{
r_sq_from_haplo(control_haps)
}
\arguments{
\item{control_haplo_file}{A Hapgen2 control haplotype file (usually with a
.haps file extension).}
}
\value{
A square matrix with dimension equal to the number of snps in the
   supplied in \code{control_haplo_file} haplotype file.
}
\description{
Calculate the r-squared matrix from control haplotype data from Hapgen2.
}
