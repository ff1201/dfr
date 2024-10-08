% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.sgl}
\alias{print.sgl}
\title{Prints information for one of the following object types: \code{"sgl"}, \code{"sgl_cv"}.}
\usage{
\method{print}{sgl}(x, ...)
}
\arguments{
\item{x}{Object of one of the following classes: \code{"sgl"}, \code{"sgl_cv"}.}

\item{...}{further arguments passed to base function.}
}
\value{
A summary of the model fit(s).
}
\description{
Prints out useful metric from a model fit.
}
\examples{
# specify a grouping structure
groups = c(rep(1:20, each=3),
          rep(21:40, each=4),
          rep(41:60, each=5),
          rep(61:80, each=6),
          rep(81:100, each=7))
# generate data
data = sgs::gen_toy_data(p=500, n=400, groups = groups, seed_id=3)
# run DFR-SGL 
model = dfr_sgl(X = data$X, y = data$y, groups = groups, type="linear", lambda = 1, alpha=0.95, 
standardise = "l2", intercept = TRUE, verbose=FALSE)
# print model
print(model)
}
\seealso{
\code{\link[=dfr_sgl]{dfr_sgl()}}, \code{\link[=dfr_sgl.cv]{dfr_sgl.cv()}}, \code{\link[=dfr_adap_sgl]{dfr_adap_sgl()}}, \code{\link[=dfr_adap_sgl.cv]{dfr_adap_sgl.cv()}}

Other SGL-methods: 
\code{\link{dfr_adap_sgl}()},
\code{\link{dfr_adap_sgl.cv}()},
\code{\link{dfr_sgl}()},
\code{\link{dfr_sgl.cv}()},
\code{\link{plot.sgl}()},
\code{\link{predict.sgl}()}
}
\concept{SGL-methods}
