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
Canovas, F., De la Rua, P. Serrano, J. and Galian, J. (2002). Mitochondrial DNA variability in Apis mellifera iberica populations from Galicia. Archivos de Zootecnia, 51: 441-448

Canovas, F., De la Rua, P. Serrano, J. and Galian, J. (2008). Geographical patterns of mitochondrial DNA variation in Apis mellifera iberiensis (Hymenoptera: Apidae). Journal of Zoological Systematics and Evolutionary Research, 46(1): 23-30

Canovas, F., De la Rua, P. Serrano, J. and Galian, J. (2014). Analysing and describing a contact area of two distinct evolutionary units. An ecological perspective. Journal of Insect Conservation. DOI: 10.1007/s10841-014-9701-1.
}

\examples{
data(apis.enirg)
plot(apis.enirg$presences[, 1:2],
     cex=apis.enirg$presences[, 3] / 10, pch=19)
}

\keyword{datasets}
