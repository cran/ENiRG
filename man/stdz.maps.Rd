\name{stdz.maps}
\alias{stdz.maps}

\title{
Standardization of raster maps
}

\description{
This function stardandize according to the mean of the values contained in a raster map, previously imported into GRASS.
}

\usage{
stdz.maps(map.names, output.names = NULL)
}

\arguments{
  \item{map.names}{
character string. The name of the raster maps to be standardized.
}
 
  \item{output.names}{
character string. The output name of the standardized raster map. 
}
}

\details{
When working with multiple rasters, it is often necessary to measures them at the same numerical scale. The standardization allows to easier find the  maximum likelihood estimates and thus compare one raster to the others. It is important to retrieve the mean and standard deviation for each raster individually and then normalize all the different raster maps. The function calculates the univariate statistics (minimum and maximum cell values, range, arithmetic mean, population variance, standard deviation, and coefficient of variation) from the non-null cells of a raster map.

The normalized raster map is calculated following the normalized equation:

\eqn{Z_{(i,j)} = \frac{x_{(i,j)} - \bar{X}}{\sigma _{X}}}

, where \eqn{Z_{(i,j)}} is the normalized value in a particular cell of cartesian coordinates (i,j) in the resolution of the map, \eqn{x_{(i,j)}} is the current value for the cell, \eqn{\bar{X}} is the mean of the map and \eqn{\sigma _{X}} is the standard deviation of the map.
}

\value{
The function returns a new raster map resulted from the normalization with the arithmetic expression. When no output names are provided, the obtained maps will be the same as the input, but adding a prefix "std_"
}

\references{
\url{http://grass.osgeo.org/grass65/manuals/r.univar.html}
\url{http://grass.osgeo.org/grass64/manuals/r.mapcalc.html}
}

\author{
Chiara  Magliozzi \email{chiara.magliozzi@libero.it}, Fernando Canovas \email{fcgarcia@ualg.pt}
}

\seealso{
\code{\link{import.egvs}}, \code{\link{enirg}}, \code{\link[spgrass6]{initGRASS}}
}

\examples{
\dontrun{
# starting GRASS session
# initGRASS("/usr/bin/grass-6.4.0", home=tempdir())
# initGRASS("C:/GRASS", home=tempdir())

## Path of the file to be imported
file1 <- paste(system.file(package = "ENiRG"), "/ext/tann.asc", sep="")

import.egvs(file1, "tann")

stdz.maps("tann","std_tann")

}
}

\keyword{standardization}
\keyword{raster}
