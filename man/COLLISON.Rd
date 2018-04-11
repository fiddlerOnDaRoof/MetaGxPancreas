\name{COLLISON}
\alias{COLLISON}
\docType{data}
\title{
 COLLISON Study Pancreatic Cancer Dataset
}
\description{
An expresion set object containing pancreatic cancer data
}
\format{
    eSet object
}
\details{
    more details can be found in the experimentData section of the object
}
\source{
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17891
}
\references{
Collisson et al.,Nat Med 2011
}
\examples{
COLLISONEset = loadPancreasEsets()$esets$COLLISON
experimentData(COLLISONEset)

}
\keyword{datasets}
