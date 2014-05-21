\name{enirg}
\alias{enirg}

\title{
Ecological Niche in R-Grass
}

\description{
\code{enirg} performs the Ecological Niche Factor Analysis (ENFA) by coupling R and GRASS softwares and following Hirzel \emph{et al}. (2002). 
}

\usage{
enirg(presences.table, qtegv.maps, qlegv.maps = NULL, col.w = NULL,
      scannf = TRUE, nf = 1, method = "normal", load.maps = FALSE,
      species.name = "species")
}

\arguments{
  \item{qtegv.maps}{
vector of strings, giving the names of the raster maps, containing quantitative environmental variables. 
}

  \item{qlegv.maps}{
by default is set to NULL. vector of strings, giving the names of the raster maps, containing qualitative environmental variables (see details). 
}

  \item{presences.table}{
data frame with species presence records in three columns, containing:
  \itemize{
    \item 1 - X coordinate or longitude.  
    \item 2 - Y coordinate or latitude.
    \item 3 - number of presences of observations.
  }
  }

  \item{col.w}{
vector with column weights with the same length as the number of quantitative and qualitative maps.
}

  \item{scannf}{
logical. Whether the number of factors should be introduced by the user (TRUE) or not (FALSE). A barplot with the eigenvalues of the analysis is displayed.
}

  \item{nf}{
integer, indicating the number of kept specialization axes when \code{scannf=FALSE}.
}

  \item{method}{
string character. Selection of the method to perform the calculations, "normal" for normal data sets or "large" for large data sets with huge amount of data. See also details. 
} 

  \item{load.maps}{
logical. Whether the produced maps should be uploaded into R. Note that if you are working with large data sets (\code{method="large"}), this option should be FALSE in order to avoid memory allocation problems.
}

  \item{species.name}{
string, indicating the name of the modelled species. 
}
}

\details{
The function \code{enirg} performs the Ecological Niche Factor Analysis (ENFA), following Hirzel \emph{et al.} (2002). The analysis can be carried out for both quantitative variables and qualitative variables. The function requires occurrences records of the modelled species and environmental predictors of the study area.

The user can choose the number of factors that better describes the variance and the eigenvalue computation for the factor extraction, by setting \code{scannf=TRUE}.

Depending on the extension of the study area and the environmental variables used to performed the analysis, the user can choose one of the two available methods: "normal" and "large". The first, strongly relies on the \code{spgrass6} package and thus is limited by R storage capacity and computations; on the contrary the "large" method directly interface with GRASS, allowing calculations over large areas or high resolution maps with huge amount of data (NOTE: it is only available for Linux/Unix OS, at the moment).

Global marginality is calculated as:


\eqn{ M = \frac{ \sqrt{ \sum _{i=1}^V m_i^2 } }{ 1.96 } }


where \eqn{m_i} represents the marginality of the focal species on each EGV (up to \emph{V} number of EGVs), in units of standards deviations of the global distribution.

Global specialization is calculated as:


\eqn{ \frac{ \sqrt{ \sum _{i=1}^V \lambda _i } }{ V } } 


where the eigenvalue \eqn{\lambda _i} associated to any factor expresses the amount of specialization it accounts for, i.e., the ratio of the variance of the global distribution to that of the species distribution on this axis.

If \code{load.maps} is set to TRUE, then the maps of marginality and \emph{nf} number of specialization maps will be loaded into R, by using \code{raster} library. Note that this is not recommended for large data sets.

}

\value{
\code{enirg} returns a list object of class "\code{enirg}" containing the following components:
\itemize{
\item call. Original call

\item nf. An integer indicating the number of kept specialization axes

\item cw. Column weights, a vector with p components

\item species. Species' name

\item egvs. A list of the used environmental predictiors

\item qt.egvs. A list of the used quantitative environmental predictiors

\item ql.egvs. A list of the used quantitative environmental predictiors

\item presences. A data frame with the table of observations

\item total.marginality. Global marginality

\item marginalities. Vector of marginalities for axis

\item total.specialization. Global specialization

\item specializations. Vector of specializations by axis

\item co. Column coordinates, data frame with \emph{V} rows and nf columns

\item mar. Coordinates of the marginality vector

\item m. Marginality

\item s. Vector of the eigenvalues

\item obs.li. Marginality and specialization axis for observations
}
}

\references{
Basille, M., Calenge, C., Marboutin, E., Andersen, R. \& Gaillard, J.M. (2008). Assessing habitat selection using multivariate statistics: Some refinements of the ecological-niche factor analysis. Ecological Modelling, 211, 233-240. 

Hirzel, A.H., Hausser, J., Chessel, D. \& Perrin, N. (2002). Ecological-niche factor analysis: How to compute habitat-suitability maps without absence data? Ecology, 83, 2027-2036.
}


\author{
Chiara Magliozzi \email{chiara.magliozzi@libero.it}, Fernando Canovas \email{fcgarcia@ualg.pt}, Jose Antonio Palazon-Ferrando \email{palazon@um.es}
}

\seealso{
\code{\link{stdz.maps}}, \link[adehabitatHS]{enfa}, \code{\link[spgrass6]{initGRASS}}
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

\keyword{ENFA}
\keyword{GRASS}
\keyword{niche}
