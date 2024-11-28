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
# Stop 1 dim matrices becoming vectors
old <- `[`
`[` <- function(...) { old(..., drop=FALSE) }
  
# -------------------------------------------------------------
# generalf functions
# -------------------------------------------------------------
# checks if a variable is binary
is.binary <- function(v) {
    x <- unique(v)
    length(x) - sum(is.na(x)) == 2L
}

sigmoid = function(x) {
   1 / (1 + exp(-x))
}

getGroupID <- function(group) { # from grpSLOPE package, which is no longer available on CRAN
  group.unique <- unique(group)
  n.group <- length(group.unique)
  group.id <- list()
  for (i in 1:n.group){
    id <- as.character(group.unique[i])
    group.id[[id]] <- which(group==group.unique[i])
  }
  return(group.id)
}

norm_vec <- function(x) sqrt(sum(x^2))

plot_path <- function(beta_matrix, lambdas, how_many,main){
  # Plots the fitted beta values along a lambda path
  beta_order = order(abs(beta_matrix[,dim(beta_matrix)[2]]),decreasing=TRUE)
  max_y = max(beta_matrix)/0.99
  min_y = min(beta_matrix)/0.99
  cols = rainbow(n=how_many)
  plot(x=-log(lambdas), y=beta_matrix[beta_order[1],],type='l',col=cols[1],ylim=c(min_y,max_y),lwd=2,xlab=expression(-log(lambda)),ylab="Fitted value",main=main)

  for (i in 2:how_many){
    lines(x=-log(lambdas),y=beta_matrix[beta_order[i],], type='l', col=cols[i],lwd=2)
  }
}

standardise_data <- function(X,y,standardise, intercept,num_obs,type="linear"){
  scale_pen = 1
  standardisation_occured = 0
  y_mean = 0
  X_center = rep(0,ncol(X))
  X_scale = rep(1,ncol(X))

  if (standardise == "l2") { # l2 normalisation
    X_center = apply(X,2,mean)
    X = sapply(1:ncol(X), function(i) X[,i] - X_center[i])
    X_scale = apply(X,2,function(x) norm(x,type="2"))
    if (any(X_scale==0)){
      stop("not able to standardise X as there exists at least one predictor with no variance")
    }
    X = sapply(1:ncol(X), function(i) X[,i]/X_scale[i])
    standardisation_occured = 1
    scale_pen = 1/sqrt(num_obs)
  } else if (standardise == "sd"){ # sd normalisation
    X_center = apply(X,2,mean)
    X = sapply(1:ncol(X), function(i) X[,i] - X_center[i])
    X_scale = apply(X, 2, function(x) sqrt(sum(x^2)/num_obs))
    if (any(X_scale==0)){
      stop("not able to standardise X as there exists at least one predictor with no variance")
    }
    X = sapply(1:ncol(X), function(i) X[,i]/X_scale[i])
    standardisation_occured = 1
  } else if (standardise == "l1"){ # l1 normalisation
    X_center = apply(X,2,mean)
    X = sapply(1:ncol(X), function(i) X[,i] - X_center[i])
    X_scale = colSums(abs(X))
    if (any(X_scale==0)){
      stop("not able to standardise X as there exists at least one predictor with no variance")
    }
    X = sapply(1:ncol(X), function(i) X[,i]/X_scale[i])
    standardisation_occured = 1
    scale_pen = 1/num_obs
  } else { standardisation_occured = 0 } # "none"
  if (intercept) { # center y and X
    if (type == "linear"){
    y_mean = mean(y)
    y = y - y_mean}
    if (standardisation_occured == 0 & !inherits(X,"dgCMatrix")){
      X_center = apply(X,2,mean)
      X = sapply(1:ncol(X), function(i) X[,i] - X_center[i])
    }
  }

out=c()
out$X=X
out$X_scale = X_scale
out$X_center = X_center
out$y = y
out$y_mean = y_mean
out$scale_pen = scale_pen
return(out)
}

which_groups <- function(beta, groups){
# outputs the non-zero group ids and effects from beta values
  num_groups = length(unique(groups))
  group_effects = data.frame(group_id = unique(sort(groups)), effect = rep(0,num_groups))
  grp_counter = 1
  for (group_id in unique(groups)){
    group_inds = which(groups==group_id)
    group_effects[grp_counter,]$effect = norm(beta[group_inds], type="2")
    grp_counter = grp_counter+1
  }
  selected_grp = group_effects[which(group_effects$effect!=0),]$group_id
  group_effects = as.vector(group_effects$effect)
  return(list(selected_grp,group_effects))
}

l2_group_operator = function(x,P, groupIDs,power){
    out = rep(0,length(groupIDs))
    for (g in 1:length(groupIDs)){
        out[g] = (P[g]^power)*norm(x[unlist(groupIDs[g])],type="2")
    }
    return(out)
}

soft_thresholding_operator <- function(x,thres){
    out = sign(x)*ifelse(abs(x) - thres <=0,0, abs(x) - thres)
    return(out)
}

path_shape <- function(lambda_max, path_length, min_frac){
  min_lambda = min_frac*lambda_max
  lambda_seq = exp(seq(log(lambda_max),log(min_lambda), (log(min_lambda) - log(lambda_max))/(path_length-1))) 
  return(lambda_seq)
}

reorder_group <- function(groups){
  max_grp_id = length(unique(groups))
  new_grp = rep(0,length(groups))
  all_grp_indices = as.numeric(names(table(groups)))
  for (i in 1:max_grp_id){
	  var_ids = which(groups == all_grp_indices[i]) 
	  new_grp[var_ids] = i
  }
return(new_grp)
}

# -------------------------------------------------------------
# loss and grad functions
# -------------------------------------------------------------
### gaussian
mse_loss <- function(y, Xbeta, num_obs, crossprod_mat){ # linear loss function
  return(as.double(crossprod_mat(y-Xbeta)/(2*num_obs)))
}

mse_grad <- function(y, Xbeta, num_obs){ # pseudo-gradient of loss, need to multiply by X^T
    return((Xbeta-y)/num_obs)
}

### binomial
log_loss <- function(y,Xbeta,num_obs, crossprod_mat){ # logistic loss function. stable version for y{0,1}
  return(mean((1-y)*Xbeta - logsig(Xbeta)))
}

log_grad <- function(y,Xbeta,num_obs){# stable version for y{0,1}. pseudo-gradient of loss, need to multiply by X^T
  return(expit_b(Xbeta,y)/num_obs)
}

# stable log implementations - from https://fa.bianp.net/blog/2019/evaluate_logistic/
logsig <- function(input){
  out = rep(0,length(input))
  idx_1 = input < -33
  idx_2 = input >= -33 & input < -18
  idx_3 = input >= -18 & input < 37
  idx_4 = input >= 37

  out[idx_1] = input[idx_1]
  out[idx_2] = input[idx_2] - exp(input[idx_2])
  out[idx_3] = -log1p(exp(-input[idx_3]))
  out[idx_4] = -exp(-input[idx_4])

  return(out)
}

expit_b <- function(t,b){
  out = rep(0,length(t))
  idx = t<0
  b_pos = b[idx]
  b_neg = b[!idx]
  exp_pos = exp(t[idx])
  exp_neg = exp(-t[!idx])
  out[idx] = ((1-b_pos)*exp_pos - b_pos)/(1+exp_pos)
  out[!idx] = ((1-b_neg) - b_neg*exp_neg)/(1+exp_neg)
  return(out)
}

# -------------------------------------------------------------
# path functions
# -------------------------------------------------------------
gen_path_sgl <- function(grad_vec_zero, groups, groupIDs, alpha, w_weights=NULL, v_weights=NULL, path_length, min_frac, group_sizes){
  lambda_max = sgl_dual(x=grad_vec_zero,groupIDs,alpha,group_sizes)
  lambda_seq = path_shape(lambda_max,path_length,min_frac)
  return(lambda_seq)
}

gen_path_adap_sgl = function(grad_vec_zero, groups, groupIDs, alpha, w_weights, v_weights, path_length, min_frac, group_sizes){
  group_weights = w_weights*group_sizes
  num_grps = length(groupIDs)
  all_lambda_vals = rep(0,num_grps)
  for (i in 1:num_grps){
    grp_ids = groupIDs[[i]]
    all_lambda_vals[i] = tryCatch(expr = {uniroot(function(lambda) {norm(soft_thresholding_operator(x=grad_vec_zero[grp_ids],thres=v_weights[grp_ids]*lambda*alpha),type="2")^2 - (group_weights[i]^2)*((1-alpha)^2)*lambda^2},
            interval=c(0,100))$root}, error = function(e) {return(0)})
  }
  lambda_max = max(all_lambda_vals)
  lambda_seq = path_shape(lambda_max,path_length,min_frac)
  return(lambda_seq)
}

# -------------------------------------------------------------
# screen functions
# -------------------------------------------------------------
sgl_grp_screen <- function(grad_vec, current_beta, groupIDs, alpha, pen_var_org, pen_grp_org, lambda_new, lambda, group_weights){
  if (all(pen_var_org == 1) & all(pen_grp_org == 1)){ # checks whether to apply adaptive sgl
   tau_g = (alpha+(1-alpha)*group_weights)
   epsilon_g = (tau_g-alpha)/tau_g
   } else {
   tau_g = gamma_g(current_beta,alpha,groupIDs,pen_var_org,pen_grp_org*group_weights)
   epsilon_g = ((1-alpha)*pen_grp_org*group_weights)/tau_g
   }
   lhs_condition = rep(0,length(tau_g))
   for (i in 1:length(tau_g)){
      lhs_condition[i] = epsilon_norm(x=as.vector(grad_vec[groupIDs[[i]]]),alpha=as.numeric(1-epsilon_g[i]),R=as.numeric(epsilon_g[i]))
   }
   rhs_condition = tau_g*(2*lambda_new-lambda)
   screen_set_grp = which(lhs_condition>=rhs_condition)
   return(screen_set_grp)
}   

sgl_var_screen <- function(grad_vec, groupIDs, screen_set_grp, alpha, pen_var_org, lambda_new, lambda, active_set_var){
 screen_set_var_initial = unlist(groupIDs[screen_set_grp])
 screen_set_var_initial = screen_set_var_initial[!(screen_set_var_initial %in% active_set_var)] # remove any from active_set_var as no need to screen these
 lhs_condition = abs(grad_vec[screen_set_var_initial])
 rhs_condition = alpha*pen_var_org[screen_set_var_initial]*(2*lambda_new - lambda)
 screen_set_var = which(lhs_condition >= rhs_condition)
 screen_set_var = screen_set_var_initial[screen_set_var]
 return(screen_set_var)
 }

# kkt check
sgl_kkt_check <- function(grad_vec, current_beta, groups, groupIDs, alpha, pen_var_org, pen_grp_org, lambda, tbl_grps, machine_tol, epsilon_set_var){
  pen_grp_org_wt = sqrt(tbl_grps)*pen_grp_org
  not_in_epsilon_set_var = which(!(1:length(current_beta) %in% epsilon_set_var))
  violations_grp = rep(0,length(groupIDs))
  violations_var = rep(0,length(current_beta))
  var_to_check = as.numeric(unlist(groupIDs))[!(as.numeric(unlist(groupIDs)) %in% epsilon_set_var)]
  lhs_condition = abs(soft_thresholding_operator(grad_vec[var_to_check], thres=(1-alpha)*lambda*pen_grp_org_wt[groups[var_to_check]]))
  rhs_condition = alpha*lambda*pen_var_org[var_to_check]
  violations_var_id = which(lhs_condition - rhs_condition > sqrt(machine_tol))
  violations_var[var_to_check[violations_var_id]] = 1
  return(list(violations_var,violations_grp))
}

# -------------------------------------------------------------
# algorithm functions
# -------------------------------------------------------------
init_lipschitz <- function(f, f_grad, mult_fcn, x0, X, y, num_obs, tX, crossprod_mat){
  L0 = 1e-3
  Xx0 = mult_fcn(X,x0)
  f0 = f(y, Xx0, num_obs, crossprod_mat)
  grad0 = mult_fcn(tX,f_grad(y, Xx0, num_obs))

  x_tilde = x0 - (1 / L0)*grad0
  f_tilde = f(y, mult_fcn(X,x_tilde), num_obs, crossprod_mat) 

  for (i in 1:100){
    if (f_tilde <= f0){
      break
    } else {
      L0 = L0 * 10
      x_tilde = x0 - (1 / L0) * grad0
      f_tilde = f(y, mult_fcn(X,x_tilde), num_obs, crossprod_mat) 
    }
  }
  return(L0)
}

proxLasso <- function(y, lambda){ # Lasso proximal operator
  out = sign(y) * pmax(0, abs(y) - lambda)
  return(out)
}

proxAdapGroupLasso = function(y, lambda, group_id, num_groups){ # adaptive group lasso proximal operator
  out = rep(0,length(y))
  lambda = as.vector(lambda) # stops warning about 1-dim array
  for (i in 1:num_groups){
    #grp_idx = which(group_info == unique(group_info)[i])
    grp_idx = group_id[[i]]
    sel_y = y[grp_idx]
    norm_y = norm(sel_y,type="2")
    if (norm_y > lambda[i]){ 
      out[grp_idx] = sel_y*(1-(lambda[i]/norm_y))
    } else {
      out[grp_idx] = 0
    } 
  }
  return(out)
}

# -------------------------------------------------------------
# sgl-specific functions
# -------------------------------------------------------------
epsilon_norm <- function(x, alpha, R){ # converted to R from the GAP safe rule python code
	if (alpha == 0 & R!=0){
		return(norm(x,"2")/R)
	}
	if (R ==0){
		return(norm(as.matrix(x,ncol=1),"I")/alpha)
	}

	zx = abs(x)
	norm_inf = norm(as.matrix(x,ncol=1),"I")
	I_inf = which(zx > alpha * (norm_inf)/ (alpha+R))
	n_inf = length(I_inf)
	zx = sort(zx[I_inf],decreasing=TRUE)

	if (norm_inf == 0){
		return(0)
	}

	if (n_inf == 1){
		return(zx[1])
	}
	R2 = R^2
	alpha2=alpha^2
	R2onalpha2 = R2/alpha2
	a_k = 0
	S = 0
	S2 = 0
	break_indc = 0
	for (k in 1:(n_inf-1)){
		S = S + zx[k]
		S2 = S2 + zx[k]^2
		b_k = S2 / (zx[k + 1]^2) - 2*S / zx[k + 1] + k

		if (a_k <= R2onalpha2 & R2onalpha2 < b_k){
			j0 = k
			break_indc = 1
			break
		} 
		a_k = b_k
	}
	if (break_indc == 0){
		j0 = n_inf
		S = S+zx[n_inf]
		S2 = S2+zx[n_inf]^2
	}
		
	alpha_S = alpha*S
	j0alpha2_R2 = j0*alpha2 - R2
	if (j0alpha2_R2 == 0){
		return(S2/(2*alpha2))
	}

	delta = alpha_S^2 - S2*j0alpha2_R2
	return((alpha_S - sqrt(delta))/ j0alpha2_R2)
}

sgl_dual = function(x, groupIDs,alpha, group_weights){
    all_vals = rep(0,length(groupIDs))
    epsilon = ((1-alpha)*group_weights)/(alpha+(1-alpha)*group_weights)
    for (i in 1:length(groupIDs)){
        all_vals[i] = (epsilon_norm(x=as.vector(x[groupIDs[[i]]]),alpha=as.numeric(1-epsilon[i]),R=as.numeric(epsilon[i])))/(alpha+(1-alpha)*group_weights[i])
    }
    return(max(all_vals))
}

off_sum_product <- function(vector1,vector2){
total_sum = sum(outer(vector1, vector2))

# Exclude the matching indices products
matching_indices_sum = sum(vector1 * vector2)

result = total_sum - matching_indices_sum
return(result)
}

gamma_g <- function(current_beta, alpha, groupIDs, pen_var_org, pen_grp_org) {
  gamma_vals = rep(0,length(groupIDs))
  for (i in 1:length(groupIDs)) {
    grp_ids = groupIDs[[i]]
    current_beta_grp = current_beta[grp_ids]
    if (all(current_beta_grp == 0)) { # the limit at zero
      #current_beta_grp = rep(sqrt(.Machine$double.eps),length(grp_ids))
      middle_term = (((length(grp_ids) - 1)*alpha)/length(grp_ids))*sum(pen_var_org[grp_ids])
    } else {
      middle_term = (alpha*off_sum_product(pen_var_org[grp_ids],abs(current_beta_grp)))/norm(as.matrix(current_beta_grp),type="1")
    }
    gamma_vals[i] = alpha*norm(as.matrix(pen_var_org[grp_ids]),type="1") - middle_term + (1-alpha)*pen_grp_org[i]
  }
  return(gamma_vals)  
}

check_group_vector <- function(vec) {
  # Check if the vector is sorted and has no gaps
  is_sorted <- all(diff(vec) >= 0)
  has_no_gaps <- all(diff(unique(vec)) == 1)
  
  return(is_sorted && has_no_gaps)
}