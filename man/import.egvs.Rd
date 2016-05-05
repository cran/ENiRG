\name{import.egvs}
\alias{import.egvs}

\title{
Raster file importation
}

\description{
This function allows to import raster files into raster maps layers supported in GRASS. Available formats are those provided by GDAL library (\url{http://www.gdal.org/formats_list.html}).
}

\usage{
import.egvs(filenames, output.names)
}

\arguments{
  \item{filenames}{
character string. Path to raster file to be imported.
}

  \item{output.names}{
character string. Name of output raster map that will be created in GRASS environment.
}
}

\details{
Datum and projection of source dataset must match projection of current GRASS location, otherwise a report error message can be generated.
}

\value{
A raster map layer will be uploaded into the opened GRASS environment.
}

\references{
\url{http://grass.osgeo.org/grass64/manuals/r.in.gdal.html}
}

\author{
Chiara  Magliozzi \email{chiara.magliozzi@libero.it}, Fernando Canovas \email{fcgarcia@ualg.pt}
}

\seealso{
\code{\link{list.maps}}, \code{\link[rgrass7]{initGRASS}}
}

\examples{
\dontrun{
# starting GRASS session
initGRASS("/usr/bin/grass-7.0.0", home=tempdir())
initGRASS("C:/GRASS", home=tempdir())

## Path to file
file1 <- paste(system.file(package = "ENiRG"), "/ext/tann.asc", sep="")

import.egvs(file1, "tann")

list.maps()

}
}

\keyword{raster}
\keyword{map}
