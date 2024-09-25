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

#' Plot models of the following object types: `"sgl"`, `"sgl_cv"`.
#'
#' Plots the pathwise solution of a cross-validation fit, from a call to one of the following: [dfr_sgl()], [dfr_sgl.cv()], [dfr_adap_sgl()], [dfr_adap_sgl.cv()].
#'
#' @param x Object of one of the following classes: \code{"sgl"}, \code{"sgl_cv"}..
#' @param how_many Defines how many predictors to plot. Plots the predictors in decreasing order of largest absolute value.
#' @param ... further arguments passed to base function.
#' 
#' @seealso [dfr_sgl()], [dfr_sgl.cv()], [dfr_adap_sgl()], [dfr_adap_sgl.cv()]
#' @family SGL-methods
#' 
#' @return A list containing:
#' \item{response}{The predicted response. In the logistic case, this represents the predicted class probabilities.}
#' \item{class}{The predicted class assignments. Only returned if type = "logistic" in the model object.}
#'
#' @examples
#' # specify a grouping structure
#' groups = c(1,1,2,2,3)
#' # generate data
#' data = sgs::gen_toy_data(p=5, n=4, groups = groups, seed_id=3,signal_mean=20,group_sparsity=1)
#' # run DFR-SGL
#' model = dfr_sgl(X = data$X, y = data$y, groups=groups, type = "linear", 
#' path_length = 20, alpha = 0.95, 
#' min_frac = 0.05, standardise="l2",intercept=TRUE,verbose=FALSE)
#' plot(model, how_many = 10)

#' @method plot sgl
#' @export
plot.sgl <- function(x, how_many = 10, ...){ 
  beta_matrix = as.matrix(x$beta)
  plot_path(beta_matrix=beta_matrix,lambdas=x$lambda,how_many=how_many,main="Pathwise solution")
}

#' @method plot sgl_cv
#' @export
plot.sgl_cv <- function(x, how_many = 10, ...){ 
  beta_matrix = as.matrix(x$all_models$beta)
  plot_path(beta_matrix=beta_matrix,lambdas=x$errors$lambda,how_many=how_many,main="Pathwise solution")
}