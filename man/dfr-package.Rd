% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dfr-package.R
\docType{package}
\name{dfr-package}
\alias{dfr}
\alias{dfr-package}
\title{dfr: Dual Feature Reduction for SGL}
\description{
Implementation of the Dual Feature Reduction (DFR) approach for the Sparse Group Lasso (SGL) and the Adaptive Sparse Group Lasso (aSGL) (Feser and Evangelou (2024) \doi{10.48550/arXiv.2405.17094}). The DFR approach is a feature reduction approach that applies strong screening to reduce the feature space before optimisation, leading to speed-up improvements for fitting SGL (Simon et al. (2013) \doi{10.1080/10618600.2012.681250}) and aSGL (Mendez-Civieta et al. (2020) \doi{10.1007/s11634-020-00413-8} and Poignard (2020) \doi{10.1007/s10463-018-0692-7}) models. DFR is implemented using the Adaptive Three Operator Splitting (ATOS) (Pedregosa and Gidel (2018) \doi{10.48550/arXiv.1804.02339}) algorithm, with linear and logistic SGL models supported, both of which can be fit using k-fold cross-validation. Dense and sparse input matrices are supported.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ff1201/dfr}
  \item Report bugs at \url{https://github.com/ff1201/dfr/issues}
}

}
\author{
\strong{Maintainer}: Fabio Feser \email{ff120@ic.ac.uk} (\href{https://orcid.org/0009-0007-3088-9727}{ORCID})

}
\keyword{internal}
