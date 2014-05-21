\name{list.maps}
\alias{list.maps}

\title{
List available raster and/or vector maps
}

\description{
The function gives the list of maps in the current GRASS session, matching with a user-specified pattern.
}

\usage{
list.maps(type = c("rast", "vect"), prefix = "*")
}

\arguments{

\item{type}{
character string from "rast" (raster maps), or "vect" (vectorial maps).
}

\item{prefix}{
character string. It states the pattern to search for.
}

}

\details{
The argument \code{prefix} allows to list all maps matching a regular expression or pattern, i.e., "std_*" will search for all the maps starting with "std_".
}

\value{
Return a vector of all the maps available in the current GRASS session or matching a regular expression.
}

\references{
\url{http://grass.osgeo.org/grass65/manuals/g.mlist.html}
}


\author{
Chiara  Magliozzi \email{chiara.magliozzi@libero.it}, Fernando Canovas \email{fcgarcia@ualg.pt} 
}

\seealso{
\code{\link{map.info}}, \code{\link{import.egvs}}, \code{\link[spgrass6]{initGRASS}}
}

\examples{
\dontrun{
# starting GRASS session
initGRASS("/usr/bin/grass-6.4.0", home=tempdir())
initGRASS("C:/GRASS", home=tempdir())

## Path of the files to be read
file1 <- paste(system.file(package = "ENiRG"),"/ext/tann.asc", sep="")
import.egvs(file1, "tann")

list.maps("rast")

}
}

\keyword{vector}
\keyword{raster}
\keyword{metadata}
\keyword{GRASS}
