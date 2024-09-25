###############################################################################
#
#    dfr: Dual Feature Reduction for the Sparse Group Lasso and Adaptive Sparse Group Lasso
#    Copyright (C) 2024 Fabio Feser
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#' Prints information for one of the following object types: `"sgl"`, `"sgl_cv"`.
#'
#' Prints out useful metric from a model fit.
#'
#' @param x Object of one of the following classes: \code{"sgl"}, \code{"sgl_cv"}.
#' @param ... further arguments passed to base function.
#' 
#' @seealso [dfr_sgl()], [dfr_sgl.cv()], [dfr_adap_sgl()], [dfr_adap_sgl.cv()]
#' @family SGL-methods
#' 
#' @return A summary of the model fit(s). 
#' 
#' @examples
#' # specify a grouping structure
#' groups = c(rep(1:20, each=3),
#'           rep(21:40, each=4),
#'           rep(41:60, each=5),
#'           rep(61:80, each=6),
#'           rep(81:100, each=7))
#' # generate data
#' data = sgs::gen_toy_data(p=500, n=400, groups = groups, seed_id=3)
#' # run DFR-SGL 
#' model = dfr_sgl(X = data$X, y = data$y, groups = groups, type="linear", lambda = 1, alpha=0.95, 
#' standardise = "l2", intercept = TRUE, verbose=FALSE)
#' # print model
#' print(model)

#' @method print sgl
#' @export
print.sgl <- function(x, ...){ 
  num.nonzero <- if(x$intercept){apply(x$beta,2, function(z){sum(z != 0)-1})}else{apply(x$beta,2, function(z){sum(z != 0)})}
  cat("\n regression type: ", x$type, "\n\n")
  print(cbind(lambda = x$lambda, num.nonzero = num.nonzero, convergence = x$success))
}

#' @method print sgl_cv
#' @export
print.sgl_cv <- function(x, ...){ 
  cat("\n regression type: ", x$type, "\n\n")
  print(cbind(lambda = x$errors$lambda, error = x$errors$error_criteria, num.nonzero = x$errors$num_non_zero))
}
