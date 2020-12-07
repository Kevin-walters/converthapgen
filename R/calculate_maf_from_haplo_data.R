#' Calculate the MAF from a single Hapgen2 control haplotype data file.
#' @param control_haplo_file A file containing Hapgen2 control haplotype data.
#' @return A vector of minor allele frequencies of length equal to the number of SNPs.
#' @export

calculate_maf_from_haplo_data <- function(control_haplo_file) {
  control_haps <- utils::read.table(control_haplo_file)
  num_controls <- dim(control_haps)[2] / 2
  transposed_control_geno_mat <-
    control_haps[, seq(1, 2 * num_controls -1, by = 2)] +
    control_haps[, 2 * seq(1, num_controls)]
  control_geno_mat <- t(transposed_control_geno_mat)
  allele_freq <- colSums(control_geno_mat)/(2* num_controls)
  minor_allele_freq <- sapply(allele_freq, function(x) min(x, 1- x))
  return(minor_allele_freq)
}

#mafs <- calculate_maf_from_haplo_data("snp7rep1.controls.haps")
