#' Calculate the r-squared matrix from control haplotype data from Hapgen2.
#' @param control_haplo_file A Hapgen2 control haplotype file (usually with a
#'   .haps file extension).
#' @return A square matrix with dimension equal to the number of snps in the
#'    supplied in \code{control_haplo_file} haplotype file.
#' @export

r_sq_from_haplo <- function(control_haps){
  #control_haps <- utils::read.table(control_haplo_file)
  num_snps <- dim(control_haps)[1]
  num_haps <- dim(control_haps)[2]
  r_sq <- diag(num_snps)

  hapmat <- as.matrix(control_haps)
  pA <- rowSums(hapmat) / num_haps # proportion of 1s by SNP
  for(snp in 1:num_snps)
  {
    pB <- pA[snp]
    #add_hap <- matrix(hapmat[snp,], nrow(hapmat), ncol(hapmat), byrow = T)
    pAB <- colSums((t(hapmat)+hapmat[snp,]) == 2)  / num_haps
    r2 <- (pAB - pA * pB) ^ 2 / (pA * (1 - pA) * pB * (1 - pB))
    r2[is.na(r2)] <- 0
    r_sq[snp, ] <- r2
  }
  return(r_sq)
}#end of function
