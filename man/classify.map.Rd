\name{classify.map}
\alias{classify.map}

\title{
Classification of the HSM
}

\description{
\code{classify.map} classifies the Habitat Suitability Map (HSM) using suitability classes.}

\usage{
classify.map(enirg.results, suit.classes, output.name, load.map = FALSE)
}

\arguments{
  \item{enirg.results}{
object of class \code{\link{enirg}}
}
  \item{suit.classes}{
object of class \code{\link{boyce}}, giving the suitability classes.
}
  \item{output.name}{
string. Name for the classified map. If set to \code{NULL} (by default), suffix "_hsm_class" will be automatically added to the species name stored into the \code{\link{enirg}} object.
}
  \item{load.map}{
logical. Whether the map should be uploaded as an object of class \code{\link[raster]{raster}}
}
}


\value{
The function \code{classify.map} allows to classify the HSM according to the breaks, which are the limits for the intervals in the four suitability classes (unsuitable, marginal, suitable and optimal; Hirzel et al., 2006). Those values can be obtained by using CBI \code{\link{boyce}}. HSM classified map layer in GRASS, which can be also uploaded into R by setting \code{load.map} to TRUE.
}

\references{
Hirzel, A.H., Le Lay, G., Helfer, V., Randin, C., Guisan, A. (2006). Evaluating the ability of the habitat suitability models to predict species presences. Ecological Modelling 199, 142-152.
}


\author{
Fernando Canovas \email{fcgarcia@ualg.pt}
}


\seealso{
\code{\link{enirg}}, \code{\link{enirg.predict}}
}


\examples{
\dontrun{
# starting GRASS session
initGRASS("/usr/bin/grass-6.4.0", home=tempdir())
initGRASS("C:/GRASS", home=tempdir())

data(apis.enirg)

# presences table
lina <- apis.enirg$presences

# loading the environmental information in batch
predictor.names <- c("tann", "mxtwm", "mntcm", "rfann", "rfdm", "rfwm")
predictor.maps <- paste("std_", predictor.names, sep="")
file.names <- paste(system.file(package = "ENiRG"),
                    "/ext/", predictor.names, ".asc", sep="")

import.egvs(file.names, predictor.names)

# standardization
stdz.maps(predictor.names, predictor.maps)

# performing the Ecologigal Niche Factor Analysis (ENFA)
enirg(presences.table = lina, qtegv.maps = predictor.maps,
      species.name = "African", nf = 1, scannf = FALSE,
      method = "normal") -> apis.enfa

enirg.predict(apis.enfa, load.map = TRUE, method = "normal") -> apis.hsm

data(cbi.results)

classify.map(apis.enfa, cbi.results,
             output.name="African_hsm_class") -> apis.hsm.class
}
}

\keyword{ HSM}
\keyword{predict}
