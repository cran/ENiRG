\name{enirg.plot}
\alias{enirg.plot}

\title{
Ordination plot of \code{\link{enirg}} results 
}

\description{
Performs an ordination diagram for objects of class \code{\link{enirg}}.
}

\usage{
enirg.plot(enirg.results, mar.col = "grey", spe.col = "black",
           method = "extended", plot.egvs = TRUE, asp = FALSE,
           title = NULL)
}

\arguments{
 
  \item{enirg.results}{
an object of class \code{enirg}
}
 
  \item{mar.col}{
color for filling available space (habitat)
}

  \item{spe.col}{
color for filling used space (presences data)
}

  \item{method}{
a string character. "extended" represents available space by using the convex hull, "simplified" represents available space by a polygon, representing the range.
}
  \item{plot.egvs}{
logical. If TRUE, ecogeographical variables will be also plotted in the ordination diagram.
}
  \item{asp}{
logical. If TRUE, marginality and specialization axes will be scaled.
}
  \item{title}{
string. A title for plot.
}
}

\details{
\code{enirg.plot} displays the factorial map of pixels, by projection of the row coordinates of ENFA analysis. Available and used habitat are displayed as polygon area and clouds of points, respectively. Marginality axis is the X-axis and specialization is the Y-axis. Quantitative environmental predictors used for the analysis are represented using blue arrows and qualitative ones are represented by text labels. White dot shows the centroid of the used habitat.
}

\references{
Basille, M., Calenge, C., Marboutin, E., Andersen, R. and Gaillard, J.M. (2008) Assessing habitat selection using multivariate statistics: Some refinements of the ecological-niche factor analysis. Ecological Modelling, 211, 233--240.
}

\author{
Fernando Canovas \email{fcgarcia@ualg.pt}, Chiara Magliozzi \email{chiara.magliozzi@libero.it}
}

\examples{
\dontrun{
# starting GRASS session
initGRASS("/usr/bin/grass-7.0.0", home=tempdir())
initGRASS("C:/GRASS", home=tempdir())

data(apis.enirg)

# presences table
lina <- apis.enirg$presences

# loading the environmental information in batch
predictor.names <- c("tann", "mxtwm", "mntcm", "rfdm", "rfseas")
predictor.maps <- paste("std_", predictor.names, sep="")
file.names <- paste(system.file(package = "ENiRG"),
                    "/ext/", predictor.names, ".asc", sep="")

import.egvs(file.names, predictor.names)

# standardization
stdz.maps(predictor.names, predictor.maps)

# performing the Ecologigal Niche Factor Analysis (ENFA)
enirg(presences.table = lina, qtegv.maps = predictor.maps,
      species.name = "African", nf = 3, scannf = FALSE,
      method = "normal") -> apis.enfa

enirg.plot(apis.enfa)

}
}


\seealso{
\code{\link{enirg}}, \code{\link{import.egvs}}, \code{\link[adehabitatHS]{enfa}}, \code{\link[adehabitatHS]{scatter.enfa}}, \code{\link[rgrass7]{initGRASS}}
}

\keyword{ordination}
\keyword{enirg}
\keyword{ENFA}
