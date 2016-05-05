\name{enirg.predict}
\alias{enirg.predict}

\title{
Predicting species' suitable habitat
}

\description{
\code{enirg.predict} computes habitat suitability maps (HSM) using the Ecological Niche Factor Analysis (see \code{\link{enirg}}), and Mahalanobis distances method.
}

\usage{
enirg.predict(enirg.results, qtegv.maps = NULL, qlegv.maps = NULL,
              load.map = FALSE, method = "normal",
              prediction.name = "predicted")
}

\arguments{
  \item{enirg.results}{
object of class \code{\link{enirg}}
}
  \item{qtegv.maps}{
vector giving names of quantitative environmental variables raster maps. If set to NULL, automatically uses the same variables as the one used to perform \code{\link{enirg}}. 
}
  \item{qlegv.maps}{
by default is set to NULL. vector of strings, giving names of raster maps, containing qualitative environmental variables (see details). If set to NULL, automatically uses the same variables as the one used to perform \code{\link{enirg}}.
}
  \item{load.map}{
logical. Whether map should be uploaded as an object of class \code{\link[raster]{raster}}.
}
  \item{method}{
string. "normal" or "large".
}
  \item{prediction.name}{
string. A suffix for naming derived maps.
}

}
\details{
Function \code{enirg.predict}, bases prediction on the results obtained from \code{enirg} function. User can choose to upload new environmental variables (i.e., to make predictions under different scenarios), or to use the same predictors, which were provided to compute the ENFA analysis. In both cases position of the niche is determined computing Mahalanobis distances for each pixel from the barycentre of the niche using the row coordinates of the ENFA analysis. Computed HSM has values ranging from 0 (complete absence) to 1 (complete presence).

Results should be later classified by using \code{\link{boyce}}.

Depending on the extension of the study area and the environmental variables used to performed the analysis, User can choose one of two available methods: "normal" and "large". The first, strongly relies on the \code{rgrass7} package and thus is limited by R storage capacity and computations; on the contrary, "large" method directly interface with GRASS, allowing calculations over large areas or high resolution maps with huge amount of data, but it is only available for Linux/Unix OS, at the moment.
}

\value{
It computes HSM and stores it as raster a map layer in GRASS. It also returns a list of:
\itemize{
\item predictions. A data frame with observed and predicted values.

\item validation. A data frame with distribution of predicted values for both the observed data and the entire predicted map.

\item map. Prediction map can be also uploaded into R by setting \code{load.map} to TRUE
}
}

\references{
Hirzel, A.H., Hausser, J., Chessel, D. \& Perrin, N. (2002) Ecological-niche factor analysis: How to compute habitat-suitability maps without absence data? Ecology, 83, 2027-2036.
}


\author{
Fernando Canovas \email{fcgarcia@ualg.pt}, Chiara  Magliozzi \email{chiara.magliozzi@libero.it}
}


\seealso{
\code{\link{enirg}}, \code{\link{import.egvs}}, \code{\link[adehabitatHS]{predict.enfa}}, \code{\link[rgrass7]{initGRASS}}
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
      species.name = "African", nf = 1, scannf = FALSE,
      method = "normal") -> apis.enfa

enirg.predict(apis.enfa, load.map = TRUE, method = "normal") -> apis.hsm

}

require(raster)

# Results can be directly loaded from:
data(apis.hsm)

plot(apis.hsm$African_predicted_hsm)
contour(apis.hsm$African_predicted_hsm, add = TRUE)
}

\keyword{HSM}
\keyword{prediction}
