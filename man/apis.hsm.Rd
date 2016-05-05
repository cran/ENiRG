\name{apis.hsm}

\alias{apis.hsm}

\docType{data}

\title{
List with values returned from \code{\link{enirg.predict}}. 
}

\description{
\code{apis.hsm} is a list which stores output of prediction for the African lineage of honeybees in the Iberian Peninsula (Canovas \emph{et al}., 2002, 2008, 2014).
}

\usage{data(apis.hsm)}

\details{
A list of:
\itemize{
\item predictions. A data frame with observed and predicted values.

\item validation. A data frame with the distribution of predicted values for both the observed data and the entire predicted map.

\item map. Prediction map for the African lineage of honeybees in the Iberian Peninsula.
}
}

\references{
Canovas, F., De la Rua, P., Serrano, J. and Galian, J. (2002). Mitochondrial DNA variability in Apis mellifera iberica populations from Galicia. Archivos de Zootecnia, 51: 441-448

Canovas, F., De la Rua, P., Serrano, J. and Galian, J. (2008). Geographical patterns of mitochondrial DNA variation in Apis mellifera iberiensis (Hymenoptera: Apidae). Journal of Zoological Systematics and Evolutionary Research, 46(1): 23-30

Canovas, F., De la Rua, P., Serrano, J. and Galian, J. (2014). Analysing and describing a contact area of two distinct evolutionary units. An ecological perspective. Journal of Insect Conservation, 18(5): 927-937
}

\examples{
data(apis.hsm)

require(raster)

plot(apis.hsm$African_predicted_hsm)
contour(apis.hsm$African_predicted_hsm, add = TRUE)
}

\keyword{datasets}
