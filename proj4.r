# Group 13: Likang Xu (s2295871), Tongfei Li (s2328507), Yifan Jin (s2434130)

# Github: https://github.com/StaceyTf1999/spproj4_Newton_g13.git

# Team member contributions: 
# Likang Xu: Contributing to part of the code, and modify
# Tongfei Li: Contributing to the main structure of the code
# Yifan Jin: Contributing to part of the code, and wrote overview
# Everyone contributed roughly equally





# The function will return an approximation to the Hessian matrix by finite 
# differencing of the gradient vector when hess=NULL
null.hess <- function(theta,grad,...){
  
  n <- length(theta)
  eps <- 1e-6  # finite difference interval
  A <- matrix(0, nrow = n, ncol = n) # empty matrix to prepare for the Hessian
  
  for (i in 1:n){  # loop over paramaters from 1 to length(theta)
    theta1 <- theta  
    theta1[i] <- theta1[i] + eps   # increase theta1[i] by eps
    # finite differencing of the gradient
    A[,i] <- (grad(theta1,...) - grad(theta,...)) / eps 
    
    # make Hessian matrix symmetric
    hess <- 0.5 * (t(A) + A)  
    
    return(hess)
  }
}


# newt arguments:
# theta: vector of initial optimization parameters
# func: the objective function to minimize, which takes the form func(theta,...)
# grad: the gradient function, has the same arguments as func
# hess: the Hessian matrix function, has the same arguments as func. If 
# hess=NULL, an approximated Hessian would be obtained using function null.hess
# ...: any arguments needs to be passed in func, grad, and hess
# tol: the convergence tolerance
# fscale: an estimate of the magnitude of func near the optimum
# maxit: the maximum number of Newton iterations before giving up
# max.half: the maximum number of times a step should be halved
# eps: the finite difference interval to use for null.hess

# Given the arguments above, the function firstly calculates the values of func
# and grad using theta and arguments provided by "...". Then, judging whether the
# Hessian matrix needs to be estimated, and calculate the Hessian matrix. If the
# objective or derivatives are not finite at the initial theta, stop the function. 

# The next step is iteration: judge whether the iteration is converged at first,
# if the answer is yes, warning if the Hessian matrix is not postive definite, 
# and return a list of:
# f: the minimum value of the objective function
# theta: the value of the parameters at the minimum
# iter: the number of iterations taken to reach the minimum
# g: the gradient vector at the minimum
# Hi: the inverse of the Hessian matrix at the minimum
# In case not converged, if hess is not positive definite, perturb it to be so
# by adding a multiple of the identity matrix to it, the multiple is the ceiling
# of the absolute value of the smallest eigen value. Then, calculate the delta
# for updating theta, and judge whether the delta fails to reduce the objective
# despite trying max.half step halvings, if the answer is yes, stop the function.
# After halvings, theta=theta+delta, and func, grad, hess is updated, and move
# to next ieration.

# If the function did not converge after maxit times of iteration, a warning
# message would be given. 


newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,
                 maxit=100,max.half=20,eps=1e-6){
  # value of objective function of initial theta
  f <- func(theta, ...)
  # gradient of initial theta
  gradient <- grad(theta, ...)
  
  # call the Hessian matrix
  if (is.null(hess)){ # If Hessian is not provided in newt, namely hess=NULL
    hess <- null.hess(theta, eps, ...)  # calling null.hess
  } else {hess <- hess(theta, ...)}
  
  # judge whether objective or derivatives are finite
  if (is.infinite(f) | any(is.infinite(gradient))){
    stop('The objective or derivatives are not finite at the initial theta.')
  }
  
  # iter is the iteration time
  iter <-  0
  # n is the dimension of theta
  n <- length(theta)
  
  while (iter < maxit){
    if (max(abs(gradient)) < (abs(f)+fscale)*tol) {
      # return sentence to tell now is convergent
      cat("This iteration is converged \n")
      
      # Judge whether the Hessian is positive definite
      decom <- try(chol(hess), silent = T)
      # if the try class in error
      if("try-error" %in% class(decom)){
        # if the Hessian is not positive definite, 
        # give the warning
        warning("The Hessian is not positive definite 
              at convergence")
      } 
      
      # Cholesky decomposition of the Hessian 
      R <- chol(hess) 
      # solve the H^-1
      invH <- backsolve(R,forwardsolve(t(R),diag(rep(1,n))))
      # return a list of containing the converged iteration result
      return(list(f=f, theta=theta, iter=iter, g=gradient, Hi=invH))
    } 
    # if this iteration is not converged 
    else {  
      # if hess is not positive definite, perturb it to be so
      if ('try-error' %in% class(try(chol(hess), silent = T))){
        # the ceiling of the absolute value of the minimum eigenvalue of hess
        # in case eigen value might be integer, +0.01
        lambda <- ceiling(abs(eigen(hess)$values[-1])+0.01)
        # H+lambda*I = t(U)*Lambda*U+lambda*I = t(U)*(Lambda+lambda*I)*U
        hess <- hess + lambda*diag(nrow(hess))
      }
        
      # calculate the delta using Cholesky decomposition of 
      # positive definite delta
      R <- chol(hess)
      # H * delta = - gradient
      # t(R) * R * delta = - gradient
      delta <- backsolve(R, forwardsolve(t(R), -gradient))
      
      # If the step fails to reduce the objective despite trying max.half
      # step halvings
      half_count = 0
      # judge whether the value is reduced
      while (func(theta + delta, ...) > f | is.infinite(f) | 
             any(is.infinite(gradient))){
        # havling the delta
        if (half_count < max.half){
          delta <- delta/2
          half_count <- half_count+1
        } else {
          stop(paste('The step fails to reduce the objective despite trying ', 
                     as.character(max.half), ' step halvings'))
        }
      }
        
      # Update the theta
      # theta[new] = theta[old] + delta
      theta <- theta + delta
      # Update the value of objective function
      f <- func(theta, ...)
      # Update the gradient using new theta
      gradient <- grad(theta, ...)
      # Update the Hessian matrix
      hess <- hess(theta, ...)
      # Update the iteration number
      iter <- iter + 1
    }
  }
  
  # Checking the last iteration: whether convergence occurs 
  if (max(abs(gradient)) < (abs(f)+fscale)*tol){
    cat("This iteration is converged \n")
    # Judge whether the Hessian is positive definite
    if("try-error" %in% class(try(chol(hess), silent = T))){
      # if the Hessian is not positive definite, give the warning
      warning("The Hessian is not positive definite 
              at convergence")
    } 
    # Cholesky decomposition of the Hessian 
    R <- chol(hess) 
    # solve the H^-1
    invH <- backsolve(R,forwardsolve(t(R),diag(rep(1,n))))
    # return a list of containing the converged iteration result
    return(list(f=f, theta=theta, iter=iter, g=gradient, Hi=invH))
  } else {
      warning(paste(as.character(maxit), " is reached without convergence"))
  }
}




rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}
gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}
hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}

newt(c(0,0), rb, gb, maxit = 2)




