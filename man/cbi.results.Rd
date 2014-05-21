\name{cbi.results}

\alias{cbi.results}

\docType{data}

\title{
List with values returned from \code{\link{boyce}}. 
}

\description{
\code{cbi.results} is a list which stores the output of the continuous Boyce classification performed with \code{\link{boyce}}.
}

\usage{data(cbi.results)}

\details{
The source dataset refers to the classification of the niche analysis performed on \code{\link{apis.enirg}} data.
}

\references{
Boyce, M.S.,Vernier, P.R.,Nielsen,S.E.,Schmiegelow, F.K.A. (2002). Evaluating resource selection functions. Ecological Modelling 157, 281-300.

Fielding, A., Bell, J. (1997). A review of methods for the assessment of  prediction errors in conservation presence/absence models. Environmental Conservation, 24, 38-49.

Hirzel, A.H., Le Lay, G., Helfer, V., Randin, C., Guisan, A. (2006). Evaluating the ability of the habitat suitability models to predict species presences. Ecological Modelling 199, 142-152. 

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

data(cbi.results)
cbi.results
}

\keyword{datasets}
