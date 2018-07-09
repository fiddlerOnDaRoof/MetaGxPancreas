\name{VANDENBROECK}
\alias{VANDENBROECK}
\docType{data}
\title{
 VANDENBROECK Study Pancreatic Cancer Dataset
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
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42952
}
\references{
NA
}
\examples{
VANDENBROECKEset = loadPancreasEsets()$esets$VANDENBROECK
experimentData(VANDENBROECKEset) 
}
\keyword{datasets}