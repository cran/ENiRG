\name{list.maps}
\alias{list.maps}

\title{
List available raster and/or vector maps
}

\description{
This function gives a list of maps in the current GRASS session, matching with a user-specified pattern.
}

\usage{
list.maps(prefix = "*")
}

\arguments{

\item{prefix}{
character string. It states a pattern to search for.
}

}

\details{
Argument \code{prefix} allows to list all maps matching a regular expression or pattern, i.e., "std_*" will search for all maps starting with "std_". See also \url{https://en.wikipedia.org/wiki/Regular_expression} regular expressions (aka regex) - from Wikipedia, the free encyclopedia. 
}

\value{
Return a list of all both vectorial and raster maps available in the current GRASS session, matching a regular expression.
}

\references{
\url{https://grass.osgeo.org/grass70/manuals/g.list.html}
}


\author{
Chiara  Magliozzi \email{chiara.magliozzi@libero.it}, Fernando Canovas \email{fcgarcia@ualg.pt} 
}

\seealso{
\code{\link{map.info}}, \code{\link{import.egvs}}, \code{\link[rgrass7]{initGRASS}}
}

\examples{
\dontrun{
# starting GRASS session
initGRASS("/usr/bin/grass-7.0.0", home=tempdir())
initGRASS("C:/GRASS", home=tempdir())

## Path of the files to be read
predictor.names <- c("tann", "mntcm", "mxtwm", "rfann", "rfdm", "rfwm")
file.names <- paste(system.file(package = "ENiRG"),
                    "/ext/", predictor.names, ".asc", sep="")

import.egvs(file.names, predictor.names)

list.maps("rf*")
list.maps("*ann")
}
}

\keyword{vector}
\keyword{raster}
\keyword{metadata}
\keyword{GRASS}
