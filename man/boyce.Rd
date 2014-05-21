\name{boyce}
\alias{boyce}

\title{
Habitat suitability classification using continuous Boyce index (CBI).
}

\description{
The \code{boyce.classification} function aims to investigate how accurately the map obtained from the function \code{\link{enirg.predict}}, by using the continuous Boyce index (CBI; Boyce \emph{et al}., 2002). 

The curves from the P/E ratio give the possibility to reclassify the predicted niche map  (\code{\link{enirg.predict}}) and thus distinguish unsuitable, marginal, suitable and optimal habitats. In accordance with the method proposed by Hirzel \emph{et al}.(2006), it interprets the predicted-to-expected ratio (P/E) by partitioning the habitat suitability predictions into classes and by calculating frequencies. If the model properly delineates the suitable areas for the studied species, the Spearman rank correlation coefficient of the ratio \eqn{F_i}, will be 1.

Evaluation of the habitat suitability model accuracy is made by means of 10-fold cross-validation (Fielding \& Bell, 1997), partitioning the data evenly but randomly into \code{cv.sets} partitions, defined by the user. Once the ratio is calculated, Spearman correlation coefficient allows to estimate the fitting for the predicted-to-expected relationship.

The categories allows to use the function \code{\link{classify.map}}, in order to perform the classification on the HSM (\code{\link{enirg.predict}}).
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
string. If "manual", a GUI assists the classification process. if "none", \code{categories} parameter will be used.
}

  \item{outcat}{
string. Name for object which will contain the results.}
}

\details{
The GUI allows a manual adjustment of the suitability classes. An ideal model would give a straight P/E curve. The curve shape and its confidence interval can be used to define the boundaries of habitat suitability classes (as suggested by the vertical dashed lines).
}

\value{
The function displays the predicted/expected curve shapes. Asol \code{enirg} returns a list object, containing the following components:
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
\code{\link[adehabitat]{predict.enfa}} 
}



\examples{
data(apis.hsm)

# vector of predictions for observations:
apis.predictions <- apis.hsm$predictions[, 2]

# vector of predictions for the entire predicted map:
apis.predictions.map <- as.vector(na.exclude(apis.hsm$African_hsm@data@values))

# Try with intervals:
# unsuitable = 0.25
# marginal = 0.69
# unsuitable = 0.81

boyce(prediction = apis.predictions,
      prediction.map = apis.predictions.map,
      categories = c(0, 0.25, 0.69, 0.81, 1),
      cv.sets = 10, type = "none")
}

\keyword{classification}
\keyword{HSM}
