% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_haps_to_geno.R
\name{create_geno_from_haplo_data}
\alias{create_geno_from_haplo_data}
\title{Make genotype files from Hapgen .haps files in usual n by p format. SNPs with
MAF below \code{maf_thresh} are removed. The r-squared matrix can be returned
for each dataset if required.}
\usage{
create_geno_from_haplo_data(
  case_haps_file,
  cont_haps_file,
  num_datasets = 1,
  maf_thresh = 0.01,
  causal_snp = 67,
  calc_r2 = T
)
}
\arguments{
\item{num_datasets}{The number of datasets supplied.}

\item{maf_thresh}{The threshold to apply when removing rarer SNPs. If the
threshold is zero it will just remove monomorphic SNPs.}

\item{causal_snp}{An interger indicating the causal SNP number (between 1 and
the number of snps).}

\item{case_haps}{A vector of strings of length two. The case files are
assumed to be numbered from 1 to \code{num_datsets}. The first element of
the vector specifies the string before the number (including the path if
needed) and the second element specifies the string after the number. One
of the strings can be NULL if necessary.}

\item{cont_haps}{A string vector as for \code{case_haps} but for
control file names.}
}
\value{
A list. The first list element is a 3-dimensional array of dimension
  'number of cases plus controls' by 'number of snps' by \code{num_datasets}
  containing the genotypes. They are coded as 0, 1 and 2 where 1 is
  heterozygous and 0 and 2 are homozygous (either could be the rare SNP). The
  second list element is the new causal SNP number (since SNPs with minor
  allele frequency below the\code{maf_thresh} have been removed). The third
  list element is the numbers of the SNPs that have been retained (MAF above
  the threshold). The fourth, optional, list element is the matrix of
  r-squared values.
}
\description{
Make genotype files from Hapgen .haps files in usual n by p format. SNPs with
MAF below \code{maf_thresh} are removed. The r-squared matrix can be returned
for each dataset if required.
}
