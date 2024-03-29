% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plotKMcurves}
\alias{plotKMcurves}
\title{Drawing Kaplan Meier Survival Curves Using the final survival-related PTMB.}
\usage{
plotKMcurves(
  sig_PTMB,
  sur,
  method = "Multivariate",
  returnAll = TRUE,
  pval = TRUE,
  color = NULL,
  plots = TRUE,
  palette = NULL,
  linetype = 1,
  conf.int = FALSE,
  pval.method = FALSE,
  test.for.trend = FALSE,
  surv.median.line = "none",
  risk.table = FALSE,
  cumevents = FALSE,
  cumcensor = FALSE,
  tables.height = 0.25,
  add.all = FALSE,
  ggtheme = theme_survminer()
)
}
\arguments{
\item{sig_PTMB}{The signature's PTMB matrix,which rows are samples and columns are pathways.}

\item{sur}{A nx2 data frame of samples' survival data,the first line is samples' survival event and the second line is samples' overall survival.}

\item{method}{Method must be one of "Univariate" and "Multivariate".}

\item{returnAll}{Logicalvalue.Default is TRUE. If TRUE, return the riskscore and the coefficient of cox regression.}

\item{pval}{Logical value, a numeric or a string. If logical and TRUE, the p-value is added on the plot. If numeric, than the computet p-value is substituted with the one passed with this parameter. If character, then the customized string appears on the plot.}

\item{color}{Color to be used for the survival curves.If the number of strata/group (n.strata) = 1, the expected value is the color name. For example color = "blue".If n.strata > 1, the expected value is the grouping variable name. By default, survival curves are colored by strata using the argument color = "strata", but you can also color survival curves by any other grouping variables used to fit the survival curves. In this case, it's possible to specify a custom color palette by using the argument palette.}

\item{plots}{logical value.Default is TRUE.If TRUE,plot the Kaplan Meier Survival Curves.}

\item{palette}{the color palette to be used. Allowed values include "hue" for the default hue color scale; "grey" for grey color palettes; brewer palettes e.g. "RdBu", "Blues", ...; or custom color palette e.g. c("blue", "red"); and scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty". See details section for more information. Can be also a numeric vector of length(groups); in this case a basic color palette is created using the function palette.}

\item{linetype}{line types. Allowed values includes i) "strata" for changing linetypes by strata (i.e. groups); ii) a numeric vector (e.g., c(1, 2)) or a character vector c("solid", "dashed").}

\item{conf.int}{logical value. If TRUE, plots confidence interval.}

\item{pval.method}{whether to add a text with the test name used for calculating the pvalue, that corresponds to survival curves' comparison - used only when pval=TRUE}

\item{test.for.trend}{logical value. Default is FALSE. If TRUE, returns the test for trend p-values. Tests for trend are designed to detect ordered differences in survival curves. That is, for at least one group. The test for trend can be only performed when the number of groups is > 2.}

\item{surv.median.line}{character vector for drawing a horizontal/vertical line at median survival. Allowed values include one of c("none", "hv", "h", "v"). v: vertical, h:horizontal.}

\item{risk.table}{Allowed values include:(1)TRUE or FALSE specifying whether to show or not the risk table. Default is FALSE.(2)"absolute" or "percentage". Shows the absolute number and the percentage of subjects at risk by time, respectively.(3)"abs_pct" to show both absolute number and percentage.(4)"nrisk_cumcensor" and "nrisk_cumevents". Show the number at risk and, the cumulative number of censoring and events, respectively.}

\item{cumevents}{logical value specifying whether to show or not the table of the cumulative number of events. Default is FALSE.}

\item{cumcensor}{logical value specifying whether to show or not the table of the cumulative number of censoring. Default is FALSE.}

\item{tables.height}{numeric value (in [0 - 1]) specifying the general height of all tables under the main survival plot.}

\item{add.all}{a logical value. If TRUE, add the survival curve of pooled patients (null model) onto the main plot.}

\item{ggtheme}{function, ggplot2 theme name. Default value is theme_survminer. Allowed values include ggplot2 official themes: see theme.}
}
\value{
Return a list of riskscore and coefficient of cox regression.
}
\description{
The function 'plotKMcurves' uses to draw Kaplan-Meier Survival Curves based on PTMB-related riskscore.The riskscore is generated by the signature's PTMB and the coefficient of "Univariate" or "Multivariate" cox regression.
}
\examples{
#get the path of the mutation annotation file and samples' survival data
\donttest{maf<-system.file("extdata","data_mutations_extended.txt",package = "pathwayTMB")
sur_path<-system.file("extdata","sur.csv",package = "pathwayTMB")
sur<-read.csv(sur_path,header=TRUE,row.names = 1)
#perform the function 'get_mut_matrix'
mut_matrix<-get_mut_matrix(maffile=maf,mut_fre = 0.01,is.TCGA=FALSE,sur=sur)
#perform the function `get_PTMB`
PTMB_matrix<-get_PTMB(freq_matrix=mut_matrix,genesmbol=genesmbol,gene_path=gene_path)
set.seed(1)
final_character<-get_final_signature(PTMB=PTMB_matrix,sur=sur)
#plot the K-M survival curve
plotKMcurves(t(PTMB_matrix[final_character,]),sur=sur,risk.table = TRUE)}

}
