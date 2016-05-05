\name{cbi.results}

\alias{cbi.results}

\docType{data}

\title{
List with values returned from \code{\link{boyce}}. 
}

\description{
\code{cbi.results} is a list which stores output of the continuous Boyce classification performed with \code{\link{boyce}}.
}

\usage{data(cbi.results)}

\details{
Source dataset refers to classification of the niche analysis performed on \code{\link{apis.enirg}} data.
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
      
data(cbi.results)
cbi.results
}

\keyword{datasets}
