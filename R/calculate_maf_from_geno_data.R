#' Calculate the minor allele frequency from control genotype data in
#' num_controls by num_snp format.
#' @param control_geno_data Control genotype data.
#' @return A vector of minor allele frequencies of length equal to the number of SNPs.
#' @export

calculate_maf_from_geno_data <- function(control_geno_data) {
  num_controls <- dim(control_geno_data)[1]
  allele_freq <- colSums(control_geno_data)/(2* num_controls)
  minor_allele_freq <- sapply(allele_freq, function(x) min(x, 1- x))
  return(minor_allele_freq)
}

