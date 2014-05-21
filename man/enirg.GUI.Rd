\name{enirg.GUI}
\alias{enirg.GUI}

\title{
Graphical User Interface of the \emph{ENiRG} package 
}

\description{
\code{enirg.GUI} creates a friendly GUI interface to set the functions available in the package. The interface relies on the \emph{tcltk} and \emph{fgui} routines to parse the functions and create widgets.
}

\usage{
enirg.GUI()
}

\value{
Displays the graphic user interface.
}

\author{
Fernando Canovas \email{fcgarcia@ualg.pt}, Chiara  Magliozzi \email{chiara.magliozzi@libero.it}
}


\examples{
\dontrun{
# starting GRASS session
initGRASS("/usr/bin/grass-6.4.0", home=tempdir())
initGRASS("C:/GRASS", home=tempdir())
}

enirg.GUI()
}

\keyword{enirg}
\keyword{gui}
