\name{boyce}
\alias{boyce}

\title{
Habitat suitability classification using continuous Boyce index (CBI).
}

\description{
\code{boyce.classification} function aims to investigate accuracy of the map obtained from function \code{\link{enirg.predict}}, by using the continuous Boyce index (CBI; Boyce \emph{et al}., 2002).

Curves from P/E ratio give the possibility to reclassify the predicted niche map  (\code{\link{enirg.predict}}) and thus distinguish unsuitable, marginal, suitable and optimal habitats. In accordance with the method proposed by Hirzel \emph{et al}.(2006), it interprets predicted-to-expected ratio (P/E) by partitioning habitat suitability predictions into classes and by calculating frequencies. If model properly delineates suitable areas for the studied species, Spearman rank correlation coefficient of the ratio \eqn{F_i}, will be 1.

Evaluation of habitat suitability model accuracy is made by means of n-fold cross-validation (Fielding and Bell, 1997), partitioning data evenly but randomly into \code{cv.sets} partitions. Once ratio is calculated, Spearman correlation coefficient allows to estimate fitting for the predicted-to-expected relationship.

Categories allows to use function \code{\link{classify.map}}, in order to perform a classification on the HSM (\code{\link{enirg.predict}}).
}

\usage{
boyce(prediction, prediction.map, categories = NULL, cv.sets = 10,
      type = "manual", outcat = "cbi.results")
}


\arguments{

  \item{prediction}{
vector. Predicted suitability values from observations or from a second validation data set. 
}

  \item{prediction.map}{
vector. Predicted suitability values for the entire area of study. 
}

  \item{categories}{
vector with desired categories. NULL if \code{type = "manual"}
}

  \item{cv.sets}{
integer, indicating the number of subsets to use for the cross validation.
}

  \item{type}{
string. If "manual", a GUI assists the classification process. If "none", \code{categories} parameter will be used.
}

  \item{outcat}{
string. Name for object which will contain the results.}
}

\details{
A GUI allows a manual adjustment of suitability classes. An ideal model would give a straight P/E curve. Curve shape and its confidence interval can be used to define boundaries of habitat suitability classes (as suggested by vertical dashed lines).
}

\value{
This function displays predicted/expected ratio curve shapes. Also \code{enirg.predict} returns a list object, containing the following components:
\itemize{
\item coefficients. A vector of two: spearman rank coefficient and adjusted r squared.

\item intervals. Suitability intervals for later using with \code{\link{classify.map}}
}
}


\references{

Boyce, M.S.,Vernier, P.R.,Nielsen,S.E.,Schmiegelow, F.K.A. (2002). Evaluating resource selection functions. Ecological Modelling 157, 281-300.

Fielding, A., Bell, J. (1997). A review of methods for the assessment of  prediction errors in conservation presence/absence models. Environmental Conservation, 24, 38-49.

Hirzel, A.H., Le Lay, G., Helfer, V., Randin, C., Guisan, A. (2006). Evaluating the ability of the habitat suitability models to predict species presences. Ecological Modelling 199, 142-152. 
}

\author{
Fernando Canovas \email{fcgarcia@ualg.pt}
}


\seealso{
\code{\link{enirg.predict}}, 
\code{\link[adehabitatHS]{predict.enfa}} 
}



\examples{
data(apis.hsm)

# vector of predictions for observations:
apis.predictions <- apis.hsm$predictions[, 2]

# vector of predictions for the entire predicted map:
apis.predictions.map <- as.vector(na.exclude(apis.hsm$African_predicted_hsm@data@values))

# Try with intervals:
# unsuitable = 0.65
# marginal = 0.84
# suitable = 0.96
# Note that this species has an optimal distribution
# in the study area, resulting in a wide unsuitable
# interval and narrow suitable ones.

boyce(prediction = apis.predictions,
      prediction.map = apis.predictions.map,
      categories = c(0, 0.65, 0.84, 0.96, 1),
      cv.sets = 10, type = "none")
}

\keyword{classification}
\keyword{HSM}
