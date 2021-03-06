% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadPancreasEsets.R
\name{loadPancreasEsets}
\alias{loadPancreasEsets}
\title{Function to load pancreas cancer expression sets from the Experiment Hub}
\usage{
loadPancreasEsets(removeDuplicates = TRUE, quantileCutoff = 0,
  rescale = FALSE, minNumberGenes = 0, minNumberEvents = 0,
  minSampleSize = 0, removeSeqSubset = TRUE, keepCommonOnly = FALSE,
  imputeMissing = FALSE)
}
\arguments{
\item{removeDuplicates}{remove patients with a Spearman correlation greater than or equal to 0.98 with other patient expression profiles (default TRUE)}

\item{quantileCutoff}{A numeric between 0 and 1 specifying to remove genes with standard deviation below the required quantile (default 0)}

\item{rescale}{apply centering and scaling to the expression sets (default FALSE)}

\item{minNumberGenes}{an integer specifying to remove expression sets with less genes than this number (default 0)}

\item{minNumberEvents}{an integer specifying how man survival events must be in the dataset to keep the dataset (default 0)}

\item{minSampleSize}{an integer specifying the minimum number of patients required in an eset (default 0)}

\item{removeSeqSubset}{currently only removes the ICGSSEQ dataset as it contains the same patients as the ICGS microarray dataset (defeault TRUE, currently just ICGSSEQ)}

\item{keepCommonOnly}{remove probes not common to all datasets (default FALSE)}

\item{imputeMissing}{remove patients from datasets with missing expression values}
}
\value{
a list with 2 elements. The First element named esets contains the datasets. The second element named duplicates contains
a vector with patient IDs for the duplicate patients (those with  Spearman correlation greater than or equal to 0.98 with other patient expression profiles).
}
\description{
This function returns pancreas cancer datasets from the hub and a vector of patients from the datasets that are most likely duplicates
}
\examples{

esetsAndDups = loadPancreasEsets()
}
