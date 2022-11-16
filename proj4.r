# Group 13: Likang Xu (s2295871), Tongfei Li (s2328507), Yifan Jin (s2434130)

# Github: https://github.com/StaceyTf1999/spproj4_Newton_g13.git

# Team member contributions: 
# Likang Xu: Contributing to part of the code, and modify
# Tongfei Li: Contributing to the main structure of the code
# Yifan Jin: Contributing to part of the code, and wrote overview
# Everyone contributed roughly equally.


# The practical is to write a function, newt, to implement the Newton minimization 
# method. The basic idea of Newton's method in optimization is that firstly, 
# we should input our initial point theta_0. When given objective function (f), 
# gradient function (g) and Hessian function (h), we can obtain delta_k by solving  
# the equation h(theta_k)*delta_k = -g(theta_k) so that we can keep on iteration  
# by theta_k+1 = theta_k + delta_k until certain stopping criterion is met.
# In newt, we use max(abs(g(theta_k)) < (abs(theta_k)+fscale)*tol as stopping 
# criterion, where fscale is a rough estimate of the magnitude of f near the 
# optimum and tol is the convergence tolerance.

# In our function, both objective function and gradient function are provided 
# while Hessian function is optional. If Hessian matrix is not provided, 
# we use finite differencing of the gradient vector to estimate the Hessian.
# Given a zero n*n matrix A (n is the length of the input theta), for every 
# column Ai in A, we have Ai = (g(theta1) - g(theta)) / eps, where
# theta1[i] = theta[i] + eps, eps is the finite difference intervals. Therefore,
# we can get estimated Hessian A and use (t(A)+A)/2 to make it symmetric.
# Based on that, we build a function 'null.hess' before 'newt' to estimate it.

# Note that we use Cholesky decomposition of h(theta_k) and calculate the inverse
# of h(theta_k) through forwardsolve() and backsolve() to avoid the high 
# computational complexity of solve(). However, h(theta_k) (denotes H below) may
# not be positive definite, so we find a way to add a multiple of the identity 
# matrix to it. More specifically, given that H is real symmetric, we does 
# eigen-decomposition: H = U*Λ*U^{T}. Then, we can obtain the minimal eigenvalue 
# λmin (obviously it's negative). Suppose λ = |λmin|+1e-6 then
# λ+λmin >= 1e-6 >0, namely Λ+λI is positive definite. Hence, 
# H+λI = U*Λ*U^{T} + λI = U*(Λ+λI)*U^{T} is now positive definite.

# What's more, the function also issue warnings in the following 4 cases: 
# 1. If the objective or derivatives are not finite at theta_0
# 2. If the step fails to reduce the objective though trying max.half step halvings 
# 3. If maxit (the maximum number of iterations to try before giving up) 
#    is reached without convergence
# 4. If the Hessian is not positive definite at convergence.



# The function will return an approximation to the Hessian matrix by finite 
# differencing of the gradient vector when hess=NULL.
null.hess <- function(theta,grad,eps,...){
  n <- length(theta)
  A <- matrix(0, nrow = n, ncol = n) # empty matrix to prepare for the Hessian
  
  for (i in 1:n){  # loop over paramaters from 1 to length(theta)
    theta1 <- theta  
    theta1[i] <- theta1[i] + eps   # increase theta1[i] by eps
    # finite differencing of the gradient
    A[,i] <- (grad(theta1,...) - grad(theta,...)) / eps 
  }
  # make Hessian matrix symmetric and return
  return(0.5*(t(A)+A))
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
# if the answer is yes, warning if the Hessian matrix is not positive definite, 
# and return a list of:
# f: the minimum value of the objective function
# theta: the value of the parameters at the minimum
# iter: the number of iterations taken to reach the minimum
# g: the gradient vector at the minimum
# Hi: the inverse of the Hessian matrix at the minimum. If Hi is not positive 
# definite, the function will not return Hi and a warning message will be given.

# In case not converged, if hess is not positive definite, perturb it to be so
# by adding a multiple of the identity matrix to it, the multiple is the smallest
# eigen value plus 1e-6. Then, calculate the delta for updating theta, and judge
# whether the delta fails to reduce the objective despite trying max.half step
# halvings, if the answer is yes, stop the function. After halvings, 
# theta=theta+delta, and func, grad, hess is updated, and move to next iteration.

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
    H <- null.hess(theta, grad, eps, ...)  # calling null.hess
  } else {H <- hess(theta, ...)}
  
  
  # judge whether objective or derivatives are finite
  if (is.infinite(f) | any(is.infinite(gradient))){
    stop('The objective or derivatives are not finite at the initial theta.')
  }
  
  # iter is the iteration time
  iter <-  0
  # n is the dimension of theta
  n <- length(theta)
  
  while (iter <= maxit){
    if (max(abs(gradient)) < (abs(f)+fscale)*tol) {
      # return sentence to tell now is convergent
      cat("This iteration is converged \n")
      
      # Judge whether the Hessian is positive definite
      decom <- try(chol(H), silent = T)
      # if the try class in error
      if("try-error" %in% class(decom)){
        # if the Hessian is not positive definite, give the warning
        warning("The Hessian is not positive definite at convergence")
        # return the list without inverse Hessian
        return(list(f=f, theta=theta, iter=iter, g=gradient))
      } else {
        # Cholesky decomposition of the Hessian 
        R <- chol(H) 
        # solve the H^-1
        invH <- backsolve(R,forwardsolve(t(R),diag(rep(1,n))))
        # return a list of containing the converged iteration result
        return(list(f=f, theta=theta, iter=iter, g=gradient, Hi=invH))
      }
    }
    # if maxit is reached without convergence
    else if (iter == maxit){
      stop(paste('The maximum number of Newton iterations', 
                 as.character(maxit), "is reached without convergence"))
    }
    # if this iteration is not converged 
    else {  
      # if H is not positive definite, perturb it to be so
      if ('try-error' %in% class(try(chol(H), silent = T))){
        # the absolute value of the minimum eigenvalue of H + 1e-6
        lambda <- abs(eigen(H)$values[-1])+1e-6
        # H+lambda*I = U*Lambda*t(U)+lambda*I = U*(Lambda+lambda*I)*t(U)
        H <- H + lambda*diag(nrow(H))
      }
      
      # calculate the delta using Cholesky decomposition of 
      # positive definite delta
      R <- chol(H)
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
          stop(paste('The step fails to reduce the objective despite trying', 
                     as.character(max.half), 'step halvings'))
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
      if (is.null(hess)){
        H <- null.hess(theta, grad, eps, ...)
      } else {H <- hess(theta, ...)}
      # Update the iteration number
      iter <- iter + 1
    }
  }
}




