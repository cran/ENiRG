\name{enirg.plot}
\alias{enirg.plot}

\title{
Ordination plot of \code{\link{enirg}} results 
}

\description{
Performs the ordination diagram for objects of class \code{\link{enirg}}.
}

\usage{
enirg.plot(enirg.results, mar.col = "grey", spe.col = "black",
           method = "extended", plot.egvs = TRUE)
}

\arguments{
 
  \item{enirg.results}{
an object of class \code{enirg}
}
 
  \item{mar.col}{
the color for filling the available space (habitat)
}

  \item{spe.col}{
the color for filling the used space (presences data)
}

  \item{method}{
a string character. "Extended" represents a polygon of the available space, "simplified" represents a clouds of points for representing the used space, with size proportional to the abundance.
}
  \item{plot.egvs}{
logical. If TRUE, ecogeographical variables will be also plotted in the ordination diagram.
}
}

\details{
\code{enirg.plot} displays the factorial map of pixels, by projection of the row coordinates of ENFA analysis. The available and the used habitat are displayed as polygon area or clouds of points. The marginality axis is the X-axis and the specialization is the Y-axis. The quantitative environmental predictors used for the analysis are represented using blue arrows and qualitative ones are represented by text labels.
}

\references{
Basille, M., Calenge, C., Marboutin, E., Andersen, R. \& Gaillard, J.M. (2008) Assessing habitat selection using multivariate statistics: Some refinements of the ecological-niche factor analysis. Ecological Modelling, 211, 233--240.
}

\author{
Fernando Canovas \email{fcgarcia@ualg.pt}, Chiara  Magliozzi \email{chiara.magliozzi@libero.it}
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

enirg.plot(apis.enirg, method = "simplified")

}
}


\seealso{
\code{\link{enirg}}, \code{\link{import.egvs}}, \code{\link[adehabitat]{enfa}}, \code{\link[adehabitat]{scatter.enfa}}, \code{\link[spgrass6]{initGRASS}}
}

\keyword{ordination}
\keyword{enirg}
\keyword{ENFA}
