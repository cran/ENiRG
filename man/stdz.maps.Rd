\name{stdz.maps}
\alias{stdz.maps}

\title{
Standardization of raster maps
}

\description{
This function standardizes according to the mean of the values contained in a raster map, previously imported into GRASS.
}

\usage{
stdz.maps(map.names, output.names = NULL)
}

\arguments{
  \item{map.names}{
character string. Name of a raster maps to be standardized.
}
 
  \item{output.names}{
character string. Output name of the standardized raster map. 
}
}

\details{
When working with multiple rasters, it is often necessary to measures them at the same numerical scale. The standardization allows to easier find maximum likelihood estimates and thus compare one raster to the others. It is important to retrieve the mean and standard deviation for each raster individually and then normalize all the different raster maps. This function calculates univariate statistics (minimum and maximum cell values, range, arithmetic mean, population variance, standard deviation, and coefficient of variation) from the non-null cells of a raster map.

Normalized raster map is then calculated following this normalized equation:

\eqn{Z_{(i,j)} = \frac{x_{(i,j)} - \bar{X}}{\sigma _{X}}}

, where \eqn{Z_{(i,j)}} is the normalized value in a particular cell of cartesian coordinates (i,j) in the resolution of the map, \eqn{x_{(i,j)}} is the current value for that cell, \eqn{\bar{X}} is the mean of the map and \eqn{\sigma _{X}} is the standard deviation of the map.
}

\value{
This function returns a new raster map resulted from the normalization. When no output names are provided, obtained maps will be the same as the input, but adding a prefix "std_"
}

\author{
Chiara  Magliozzi \email{chiara.magliozzi@libero.it}, Fernando Canovas \email{fcgarcia@ualg.pt}
}

\seealso{
\code{\link{import.egvs}}, \code{\link{enirg}}, \code{\link[rgrass7]{initGRASS}}
}

\examples{
\dontrun{
# starting GRASS session
# initGRASS("/usr/bin/grass-7.0.0", home=tempdir())
# initGRASS("C:/GRASS", home=tempdir())

## Path to the file to be imported
file1 <- paste(system.file(package = "ENiRG"), "/ext/tann.asc", sep="")

import.egvs(file1, "tann")

stdz.maps("tann","std_tann")

}
}

\keyword{standardization}
\keyword{raster}
