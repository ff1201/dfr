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

#' Predict using one of the following object types: `"sgl"`, `"sgl_cv"`.
#'
#' Performs prediction from one of the following fits: [dfr_sgl()], [dfr_sgl.cv()], [dfr_adap_sgl()], [dfr_adap_sgl.cv()]. The predictions are calculated for each \code{"lambda"} value in the path.
#'
#' @param object Object of one of the following classes: \code{"sgl"}, \code{"sgl_cv"}.
#' @param x Input data to use for prediction.
#' @param ... further arguments passed to stats function.
#' 
#' @seealso [dfr_sgl()], [dfr_sgl.cv()], [dfr_adap_sgl()], [dfr_adap_sgl.cv()]
#' @family SGL-methods
#' 
#' @return A list containing:
#' \item{response}{The predicted response. In the logistic case, this represents the predicted class probabilities.}
#' \item{class}{The predicted class assignments. Only returned if type = "logistic" in the \code{"sgl"} or \code{"sgl_cv"} object.}
#'
#' @examples
#' # specify a grouping structure
#' groups = c(1,1,1,2,2,3,3,3,4,4)
#' # generate data
#' data = sgs::gen_toy_data(p=10, n=5, groups = groups, seed_id=3,group_sparsity=1)
#' # run DFR-SGL 
#' model = dfr_sgl(X = data$X, y = data$y, groups = groups, type="linear", lambda = 1, alpha=0.95, 
#' standardise = "l2", intercept = TRUE, verbose=FALSE)
#' # use predict function
#' model_predictions = predict(model, x = data$X)

#' @method predict sgl
#' @export
predict.sgl <- function(object, x, ...){
  if (object$type=="linear"){
    if (object$intercept){
      predictions = apply(object$beta, 2, function(z) sgs::arma_mv(Matrix::cbind2(1,x),as.vector(z)))
    } else {
      predictions = apply(object$beta, 2, function(z) sgs::arma_mv(x,as.vector(z)))
    }
  } else if (object$type == "logistic"){
    predictions = c()
    if (object$intercept){
      predictions$response = apply(object$beta, 2, function(z) sigmoid(sgs::arma_mv(Matrix::cbind2(1,x),as.vector(z))))
    } else {
      predictions$response = apply(object$beta, 2, function(z) sigmoid(sgs::arma_mv(x,as.vector(z))))
    }
    predictions$class = ifelse(predictions$response>0.5,1,0)
  } 
  return(predictions)
}

#' @method predict sgl_cv
#' @export
predict.sgl_cv <-  function(object, x, ...){
  if (object$fit$type=="linear"){
    if (object$fit$intercept){
      predictions = apply(object$all_models$beta, 2, function(z) sgs::arma_mv(Matrix::cbind2(1,x),as.vector(z)))
    } else {
      predictions = apply(object$all_models$beta, 2, function(z) sgs::arma_mv(x,as.vector(z)))
    }
  } else if (object$fit$type == "logistic"){
    predictions = c()
    if (object$fit$intercept){
      predictions$response = apply(object$all_models$beta, 2, function(z) sigmoid(sgs::arma_mv(Matrix::cbind2(1,x),as.vector(z))))
    } else {
      predictions$response = apply(object$all_models$beta, 2, function(z) sigmoid(sgs::arma_mv(x,as.vector(z))))
    }
    predictions$class = ifelse(predictions$response>0.5,1,0)
  }
  return(predictions)
}