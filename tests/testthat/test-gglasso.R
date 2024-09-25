# gglasso doesn't apply standardisation, so can't compare directly with standardisation
test_that("sgl solution reduces to grplasso when alpha=0, with uneven groups, with no intercept and standardisation", {
  skip_if_not_installed("gglasso")
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  p=100
  n=50
  data = sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = TRUE,groups = groups,group_sparsity=0.2,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  lambda = 0.8
  sgl = dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=0, standardise="none", intercept=FALSE)
  gglasso = gglasso::gglasso(x=X,y=y, group=groups, lambda=lambda,intercept=FALSE,eps=1e-5, loss="ls",pf=sqrt(table(groups)))
   
  expect_equivalent(gglasso$beta,
    as.matrix(sgl$beta),
    tol = 1e-3
  )

})
  
test_that("sgl solution reduces to grplasso when alpha=0, with even groups, with no intercept and standardisation", {
  skip_if_not_installed("gglasso")
  groups = rep(1:20,each=5)
  p=100
  n=50
  data = sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = TRUE,groups = groups,group_sparsity=0.2,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  lambda = 0.8
  sgl = dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, standardise="none", alpha=0, intercept=FALSE)
  gglasso = gglasso::gglasso(x=X,y=y, group=groups, lambda=lambda,intercept=FALSE,eps=1e-5, loss="ls",pf=sqrt(table(groups)))
   
  expect_equivalent(gglasso$beta,
    as.matrix(sgl$beta),
    tol = 1e-3
  )
})

test_that("sgl solution reduces to grplasso when alpha=0, with uneven groups, with intercept but no standardisation", {
  skip_if_not_installed("gglasso")
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  p=100
  n=50
  data = sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 4,grouped = TRUE,groups = groups,group_sparsity=0.2,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  lambda = 0.8
  sgl =  dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=0, standardise="none", ,intercept=TRUE)
  gglasso = gglasso::gglasso(x=X,y=y, group=groups, lambda=lambda,intercept=TRUE,eps=1e-5, loss="ls",pf=sqrt(table(groups)))
  
  expect_equivalent(c(gglasso$b0, gglasso$beta),
    as.matrix(sgl$beta),
    tol = 1e-3
  )

})
  
test_that("sgl solution reduces to grplasso when alpha=0, with even groups, with intercept but no standardisation", {
  skip_if_not_installed("gglasso")
  groups = rep(1:20,each=5)
  p=100
  n=50
  data = sgs::gen_toy_data(p=p,n=n,rho = 0,seed_id = 100,grouped = TRUE,groups = groups,group_sparsity=0.2,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  lambda = 0.8
  sgl =  dfr_sgl(X=X,y=y, groups=groups, type="linear", lambda=lambda,standardise="none", alpha=0, intercept=TRUE)
  gglasso = gglasso::gglasso(x=X,y=y, group=groups, lambda=lambda,intercept=TRUE,eps=1e-5, loss="ls",pf=sqrt(table(groups)))

  expect_equivalent(c(gglasso$b0, gglasso$beta),
    as.matrix(sgl$beta),
    tol = 1e-3
  )

})
  