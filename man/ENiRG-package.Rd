\name{ENiRG-package}
\alias{ENiRG-package}

\docType{package}

\title{
Ecological Niche in R-Grass
}


\description{
The package ENiRG has been designed to characterize the realized niche of the species by interfacing R software with GRASS geographical information system in order to overcome issues when working with large data sets (i.e., wide areas or high resolution). This package uses classes defined in the \code{spgrass6} package to deal with spatial data and to interface R and Grass.
}

\details{
\tabular{ll}{
Package: \tab ENiRG \cr
Type: \tab Package \cr
Version: \tab 0.1_1 \cr
Date: \tab 2014-06-02 \cr
License: \tab GPL (>=2) \cr
}

This package is organised in four main parts: 
\itemize{
  \item data preparation (\code{\link{import.egvs}}, \code{\link{list.maps}}, \code{\link{map.info}} and \code{\link{stdz.maps}});

  \item ENFA analysis (\code{\link{enirg}} and \code{\link{enirg.plot}});

  \item prediction of species' niche and HSM calculation (\code{\link{enirg.predict}});

  \item classification of suitability (\code{\link{boyce}} and \code{\link{classify.map}}).
}

\emph{ENiRG} is able to perform the analysis into GRASS, by intitiating a session within R (see also the function \code{\link{initGRASS}} from \emph{spgrass6}). \emph{ENiRG} provides an easy way to import raster maps through \code{\link{import.egvs}}, containing EGV's distribution in any of the formats supported by GDAL library (\href{http://www.gdal.org}{gdal.org}). It also provides several functions to explore \href{http://grass.osgeo.org/}{GRASS} environment, giving the available list of maps (\code{\link{list.maps}}) and metadata (\code{\link{map.info}}). Before run ENFA analysis, standardization of quantitative EGVs is suggested Hirzel \emph{et al}. (2002). \code{\link{stdz.maps}} allows the process by retrieving its univariate statistics.

The characterization of the species's niche can be done by using the main function of the package. \code{\link{enirg}} performs ENFA analysis, while the projection of the ordination diagram of the marginality and specialization (ENFA's principal components), can be computed using \code{\link{enirg.plot}}.

The function \code{\link{enirg.predict}} focuses on computing the prediction of the species' niche accordingly with the ENFA's results.

The main functions of the package offer two methods, accordingly to the amount of data computed. The user can choose one of the two available methods: "normal" and "large". The first, strongly relies on the \code{spgrass6} package and thus is limited by R storage capacity and computations; on the contrary the "large" method directly interface with GRASS, allowing calculations over large areas or high resolution maps with huge amount of data (NOTE: it is only available for Linux/Unix OS, at the moment).

\code{\link{boyce}} function aims to investigate how accurately the map obtained from the function \code{\link{enirg.predict}} is predicting the modelled species presences (Boyce \emph{et al}., 2002). The intervals estimated manually, allow the user to reclassify the predicted niche map, by using \code{\link{classify.map}}), and thus distinguishing unsuitable, marginal, suitable and optimal habitat (Hirzel \emph{et al}., 2006). Evaluation of the habitat suitability model accuracy is made by means of 10-fold cross-validation (Fielding and Bell, 1997).

The user can take advantages of the integration with other R libraries (\emph{raster}, \emph{rasterVIS}), the portability and interoperability within GRASS (i.e. efficient map storage) and can also communicate with other commonly used GIS software, such as \href{http://www.qgis.org/}{QGIS}.

A graphical user interface (GUI) allows better access to the functionalities of \emph{ENiRG} package through the function \code{link{enirg.GUI}}.
}


\author{
Fernando Canovas, Chiara Magliozzi, Jose Antonio Palazon-Ferrando, Frederico Mestre, Mercedes Gonzalez-Wanguemert

Maintainer: Fernando Garcia Canovas \email{fcgarcia@ualg.pt}
}

\references{
Basille, M., Calenge, C., Marboutin, E., Andersen, R. and Gaillard, J.M. (2008). Assessing habitat selection using multivariate statistics: Some refinements of the ecological-niche factor analysis. Ecological Modelling, 211, 233-240. 

Boyce, M.S.,Vernier, P.R.,Nielsen,S.E.,Schmiegelow, F.K.A. (2002). Evaluating resource selection functions. Ecological Modelling 157, 281-300.

Fielding, A., Bell, J. (1997). A review of methods for the assessment of  prediction errors in conservation presence/absence models. Environmental Conservation, 24, 38-49.

Hirzel, A.H., Hausser, J., Chessel, D. and Perrin, N. (2002). Ecological-niche factor analysis: How to compute habitat-suitability maps without absence data? Ecology, 83, 2027-2036.

Hirzel, A.H., Le Lay, G., Helfer, V., Randin, C., Guisan, A. (2006). Evaluating the ability of the habitat suitability models to predict species presences. Ecological Modelling 199, 142-152. 
}

\note{
The package depends on spgrass6, raster, R.utils, stats, gplots, miniGUI, ade4, tcltk2, fgui.
}

\seealso{
adehabitat, adehabitatHS
}

\keyword{ENiRG}
\keyword{ENFA}