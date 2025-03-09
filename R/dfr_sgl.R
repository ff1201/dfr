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

#' Fit a DFR-SGL model.
#' 
#' Sparse-group lasso (SGL) with DFR main fitting function. Supports both linear and logistic regression, both with dense and sparse matrix implementations.
#' 
#' \code{dfr_sgl()} fits a DFR-SGL model (Feser and Evangelou (2024)) using Adaptive Three Operator Splitting (ATOS) (Pedregosa and Gidel (2018)). 
#' It solves the convex optimisation problem given by (Simon et al. (2013))
#' \deqn{
#'   \frac{1}{2n} f(b ; y, \mathbf{X}) + \lambda \alpha \sum_{i=1}^{p} |b_i| + \lambda (1-\alpha)\sum_{g=1}^{m}  \sqrt{p_g} \|b^{(g)}\|_2,
#' }
#' where \eqn{f(\cdot)} is the loss function and \eqn{p_g} are the group sizes. In the case of the linear model, the loss function is given by the mean-squared error loss:
#' \deqn{
#'  f(b; y, \mathbf{X}) = \left\|y-\mathbf{X}b \right\|_2^2.
#' }
#' In the logistic model, the loss function is given by 
#' \deqn{
#' f(b;y,\mathbf{X})=-1/n \log(\mathcal{L}(b; y, \mathbf{X})).
#' }
#' where the log-likelihood is given by
#' \deqn{
#'  \mathcal{L}(b; y, \mathbf{X}) = \sum_{i=1}^{n}\left\{y_i b^\intercal x_i - \log(1+\exp(b^\intercal x_i)) \right\}.
#' }
#' SGL can be seen to be a convex combination of the lasso and group lasso, balanced through \code{alpha}, such that it reduces to the lasso for \code{alpha = 1} and to the group lasso for \code{alpha = 0}.
#' By applying both the lasso and group lasso norms, SGL shrinks inactive groups to zero, as well as inactive variables in active groups.
#' DFR uses the dual norm (the \eqn{\epsilon}-norm) and the KKT conditions to discard features at \eqn{\lambda_k} that would have been inactive at \eqn{\lambda_{k+1}}.
#' It applies two layers of screening, so that it first screens out any groups that satisfy
#' \deqn{
#' \|\nabla_g f(\hat{\beta}(\lambda_{k}))\|_{\epsilon_g} \leq \tau_g(2\lambda_{k+1} - \lambda_k)
#' }
#' and then screens out any variables that satisfy
#' \deqn{
#' |\nabla_i f(\hat{\beta}(\lambda_{k}))| \leq \alpha (2\lambda_{k+1} - \lambda_k)
#' }
#' leading to effective input dimensionality reduction. See Feser and Evangelou (2024) for full details.
#' 
#' @param X Input matrix of dimensions \eqn{n \times p}{n*p}. Can be a sparse matrix (using class \code{"sparseMatrix"} from the \code{Matrix} package).
#' @param y Output vector of dimension \eqn{n}. For \code{type="linear"} should be continuous and for \code{type="logistic"} should be a binary variable.
#' @param groups A grouping structure for the input data. Should take the form of a vector of group indices.
#' @param type The type of regression to perform. Supported values are: \code{"linear"} and \code{"logistic"}.
#' @param lambda The regularisation parameter. Defines the level of sparsity in the model. A higher value leads to sparser models: 
#'   - \code{"path"} computes a path of regularisation parameters of length \code{"path_length"}. The path will begin just above the value at which the first predictor enters the model and will terminate at the value determined by \code{"min_frac"}.
#'   - User-specified single value or sequence. Internal scaling is applied based on the type of standardisation. The returned \code{"lambda"} value will be the original unscaled value(s).
#' @param path_length The number of \eqn{\lambda} values to fit the model for. If \code{"lambda"} is user-specified, this is ignored.
#' @param min_frac Smallest value of \eqn{\lambda} as a fraction of the maximum value. That is, the final \eqn{\lambda} will be \code{"min_frac"} of the first \eqn{\lambda} value.
#' @param alpha The value of \eqn{\alpha}, which defines the convex balance between the lasso and group lasso. Must be between 0 and 1. Recommended value is 0.95.
#' @param max_iter Maximum number of ATOS iterations to perform. 
#' @param backtracking The backtracking parameter, \eqn{\tau}, as defined in Pedregosa and Gidel (2018).
#' @param max_iter_backtracking Maximum number of backtracking line search iterations to perform per global iteration.
#' @param tol Convergence tolerance for the stopping criteria.
#' @param standardise Type of standardisation to perform on \code{X}: 
#'   - \code{"l2"} standardises the input data to have \eqn{\ell_2} norms of one. When using this \code{"lambda"} is scaled internally by \eqn{1/\sqrt{n}}.
#'   - \code{"l1"} standardises the input data to have \eqn{\ell_1} norms of one. When using this \code{"lambda"} is scaled internally by \eqn{1/n}.
#'   - \code{"sd"} standardises the input data to have standard deviation of one.
#'   - \code{"none"} no standardisation applied.
#' @param intercept Logical flag for whether to fit an intercept.
#' @param screen Logical flag for whether to apply the DFR screening rules (see Feser and Evangelou (2024)).
#' @param verbose Logical flag for whether to print fitting information.
#' 
#' @return A list containing:
#' \item{beta}{The fitted values from the regression. Taken to be the more stable fit between \code{x} and \code{z}, which is usually the former. A filter is applied to remove very small values, where ATOS has not been able to shrink exactly to zero. Check this against \code{x} and \code{z}.}
#' \item{group_effects}{The group values from the regression. Taken by applying the \eqn{\ell_2} norm within each group on \code{beta}.}
#' \item{selected_var}{A list containing the indicies of the active/selected variables for each \code{"lambda"} value. Index 1 corresponds to the first column in X.}
#' \item{selected_grp}{A list containing the indicies of the active/selected groups for each \code{"lambda"} value. Index 1 corresponds to the first group entry in the \code{groups} vector. You can see the group order by running \code{unique(groups)}.}
#' \item{num_it}{Number of iterations performed. If convergence is not reached, this will be \code{max_iter}.}
#' \item{success}{Logical flag indicating whether ATOS converged, according to \code{tol}.}
#' \item{certificate}{Final value of convergence criteria.} 
#' \item{x}{The solution to the original problem (see Pedregosa and Gidel (2018)).}
#' \item{u}{The solution to the dual problem (see Pedregosa and Gidel (2018)).}
#' \item{z}{The updated values from applying the first proximal operator (see Pedregosa and Gidel (2018)).}
#' \item{screen_set_var}{List of variables that were kept after screening step for each \code{"lambda"} value. (see Feser and Evangelou (2024)).}
#' \item{screen_set_grp}{List of groups that were kept after screening step for each \code{"lambda"} value. (see Feser and Evangelou (2024)).}
#' \item{epsilon_set_var}{List of variables that were used for fitting after screening for each \code{"lambda"} value. (see Feser and Evangelou (2024)).}
#' \item{epsilon_set_grp}{List of groups that were used for fitting after screening for each \code{"lambda"} value. (see Feser and Evangelou (2024)).}
#' \item{kkt_violations_var}{List of variables that violated the KKT conditions each \code{"lambda"} value. (see Feser and Evangelou (2024)).}
#' \item{kkt_violations_grp}{List of groups that violated the KKT conditions each \code{"lambda"} value. (see Feser and Evangelou (2024)).}
#' \item{screen}{Logical flag indicating whether screening was performed.}
#' \item{type}{Indicates which type of regression was performed.}
#' \item{intercept}{Logical flag indicating whether an intercept was fit.}
#  \item{standardise}{Indicates the type of standardisation used.}
#' \item{lambda}{Value(s) of \eqn{\lambda} used to fit the model.}
#'
#' @family SGL-methods
#' 
#' @examples
#' # specify a grouping structure
#' groups = c(1,1,1,2,2,3,3,3,4,4)
#' # generate data
#' data = sgs::gen_toy_data(p=10, n=5, groups = groups, seed_id=3,group_sparsity=1)
#' # run DFR-SGL 
#' model = dfr_sgl(X = data$X, y = data$y, groups = groups, type="linear", path_length = 5, 
#' alpha=0.95, standardise = "l2", intercept = TRUE, verbose=FALSE)
#' @references Feser, F., Evangelou, M. (2024). \emph{Dual feature reduction for the sparse-group lasso and its adaptive variant}, \url{https://arxiv.org/abs/2405.17094}
#' @references Pedregosa, F., Gidel, G. (2018). \emph{Adaptive Three Operator Splitting}, \url{https://proceedings.mlr.press/v80/pedregosa18a.html}
#' @references Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2013). \emph{A Sparse-Group Lasso}, \doi{10.1080/10618600.2012.681250}
#' @export

dfr_sgl <- function(X, y, groups, type="linear", lambda="path", alpha=0.95, max_iter=5000, backtracking=0.7, max_iter_backtracking=100, tol=1e-5, standardise="l2", intercept=TRUE, path_length=20, min_frac=0.05, screen=TRUE, verbose=FALSE){
  # check group ordering to ensure no gaps in group numbers and ordered from 1 to m
  if (check_group_vector(groups)){
    reorder_id = FALSE
    ordered_grp_ids = groups
  } else {
    reorder_id = TRUE
    grp_new = reorder_group(groups)
    order_grp = order(grp_new,decreasing=FALSE)
    ordered_grp_ids = match(groups[order_grp], sort(unique(groups[order_grp])))
    X = X[,order_grp]
  }

  # Run main fitting function
  out = general_fit(X, y, ordered_grp_ids, gen_path_sgl, type, lambda, path_length, alpha, backtracking, max_iter, max_iter_backtracking, 
                    tol, min_frac, standardise, intercept, NULL, NULL, screen, verbose, NULL, NULL, model = "sgl")
  
  # put group ordering back to original
  if (reorder_id){
    out = reorder_output(out, intercept, order_grp, groups)
  }
  
  return(out)
}