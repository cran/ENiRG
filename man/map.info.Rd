\name{map.info}
\alias{map.info}

\title{
Retrieve metadata of a raster file from the current GRASS location.
}

\description{
Metadata about GRASS environment and about the map (i.e., raster map resolution, type and range) will be returned as a list. The user can choose the GRASS format to store the information, i.e. in a .txt file, or to have the output as a list in R.
}

\usage{
map.info(map.name, format = "grass")
}

\arguments{

\item{map.name}{character string. The name of the raster map stored in GRASS as returned from \code{\link{list.maps}}.}
  
\item{format}{character string, by default "grass". The preferred output format between GRASS output format ("grass") and R list ("r").}
}

\value{
Return either a list or text.
}

\references{
\url{http://grass.osgeo.org/grass65/manuals/r.info.html}
}

\author{
Chiara  Magliozzi \email{chiara.magliozzi@libero.it}, Fernando Canovas \email{fcgarcia@ualg.pt} 
}

\seealso{
\code{\link{import.egvs}}, \code{\link[spgrass6]{initGRASS}}
}

\examples{
\dontrun{
# starting GRASS session
initGRASS("/usr/bin/grass-6.4.0", home=tempdir())
initGRASS("C:/GRASS", home=tempdir())

## Path of the file to be read if we want just one of the variables
file1 <- paste(system.file(package = "ENiRG"),"/ext/tann.asc", sep="")
import.egvs(file1, "tann")

map.info("tann", format="r")

}
}

\keyword{metadata}
\keyword{raster}
\keyword{map}
