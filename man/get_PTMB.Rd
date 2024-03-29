% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{get_PTMB}
\alias{get_PTMB}
\title{Calculate the Pathway-based Tumor Mutational Burden.}
\usage{
get_PTMB(freq_matrix, genesmbol, path_mut_cutoff = 0, gene_path)
}
\arguments{
\item{freq_matrix}{The mutations matrix,generated by `get_mut_matrix`.}

\item{genesmbol}{The genes' length matrix,generated by `get_gene_length`.}

\item{path_mut_cutoff}{A threshold value(zero percent as the default value).Pathways with a given mutation frequency equal or greater than the threshold value are retained for the following analysis.}

\item{gene_path}{User input pathways geneset list.}
}
\value{
Return the Pathway-based Tumor Mutational Burden matrix.
}
\description{
The function `get_PTMB` uses to calculate the Pathway-based Tumor Mutational Burden (PTMB). PTMB is defined as pathway-based tumor mutational burden corrected by genes’ length and number.
}
\examples{
#get the path of the mutation annotation file and samples' survival data
\donttest{maf<-system.file("extdata","data_mutations_extended.txt",package = "pathwayTMB")
sur_path<-system.file("extdata","sur.csv",package = "pathwayTMB")
sur<-read.csv(sur_path,header=TRUE,row.names = 1)
#perform the function 'get_mut_matrix'
mut_matrix<-get_mut_matrix(maffile=maf,mut_fre = 0.01,is.TCGA=FALSE,sur=sur)}
#perform the function `get_PTMB`
PTMB_matrix<-get_PTMB(freq_matrix=mut_matrix,genesmbol=genesmbol,gene_path=gene_path)
}
