\name{apis.enirg}

\alias{apis.enirg}

\docType{data}

\title{
List of an object \code{enirg}. 
}

\description{
\code{apis.enirg} is a list which stores the output of an analysis performed with \code{\link{enirg}} over occurrence records of the \emph{Apis mellifera} species across the Iberian Peninsula (Canovas et al., 2002, 2008, 2014).
}

\usage{data(apis.enirg)}

\details{
The source dataset only refers to the African lineage occurrence records.
}

\references{
Canovas, F., De la Rua, P., Serrano, J. and Galian, J. (2002). Mitochondrial DNA variability in Apis mellifera iberica populations from Galicia. Archivos de Zootecnia, 51: 441-448

Canovas, F., De la Rua, P., Serrano, J. and Galian, J. (2008). Geographical patterns of mitochondrial DNA variation in Apis mellifera iberiensis (Hymenoptera: Apidae). Journal of Zoological Systematics and Evolutionary Research, 46(1): 23-30

Canovas, F., De la Rua, P., Serrano, J. and Galian, J. (2014). Analysing and describing a contact area of two distinct evolutionary units. An ecological perspective. Journal of Insect Conservation, 18(5): 927-937

Canovas, F., Magliozzi, C., Mestre, F., Palazon-Ferrando, J.A. and Gonzalez-Wanguemert,
M. (2016). ENiRG: R-GRASS interface for efficiently characterizing the ecological niche of
species and predicting habitat suitability. Ecography 39: 593-598.
}

\examples{
require(raster)

# loading map of Iberian Peninsula
mask.file <- paste(system.file(package = "ENiRG"),
                   "/ext/", "mask.asc", sep="")
plot(raster(mask.file), col="brown", legend=FALSE, xlab = "longitude",
     ylab = "latitude")

data(apis.enirg)
points(apis.enirg$presences[, 1:2],
       cex=apis.enirg$presences[, 3] / 30, pch=19)
}

\keyword{datasets}
