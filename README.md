# spproj4_Newton_g13


Programming Task: in your work group of 3, to write an R function, newt, implementing Newton’s method for
minimization of functions. Only functions in base R and its packages (including recommended packages) can
be used and your code must not install any packages. Do not use functions attach or solve. Your submitted
code should only define functions, and not have any code outside functions except for loading permitted package
(it can of course contain functions other than newt).
Specification: Your newt optimization function should operate broadly in the same way as nlm and optim, but
the purpose is to have an independent implementation: you must code the optimization yourself, not simply call
optimization code written by someone else. The function arguments should be as follows:
newt(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6)
theta is a vector of initial values for the optimization parameters.
func is the objective function to minimize. Its first argument is the vector of optimization parameters. Remaining
arguments will be passed from newt using ‘...’.
grad is the gradient function. It has the same arguments as func but returns the gradient vector of the objective
w.r.t. the elements of parameter vector.
hess is the Hessian matrix function. It has the same arguments as func but returns the Hessian matrix of the
objective w.r.t. the elements of parameter vector. If not supplied then newt should obtain an approximation to the
Hessian by finite differencing of the gradient vector (hint: (t(A)+A)/2 is exactly symmetric for any matrix A).
... any arguments of func, grad and hess after the first (the parameter vector) are passed using this.
tol the convergence tolerance.
fscale a rough estimate of the magnitude of func near the optimum - used in convergence testing.
maxit the maximum number of Newton iterations to try before giving up.
max.half the maximum number of times a step should be halved before concluding that the step has failed to
improve the objective.
eps the finite difference intervals to use when a Hessian function is not provided.
newt should return a list containing:
f the value of the objective function at the minimum.
theta the value of the parameters at the minimum.
iter the number of iterations taken to reach the minimum.
g the gradient vector at the minimum (so the user can judge closeness to numerical zero).
Hi the inverse of the Hessian matrix at the minimum (useful if the objective is a negative log likelihood).
The function should issue errors or warnings (using stop or warning as appropriate) in at least the following
cases. 1. If the objective or derivatives are not finite at the initial theta; 2. If the step fails to reduce the objective
despite trying max.half step halvings; 3. If maxit is reached without convergence; 4. If the Hessian is not
positive definite at convergence.
Convergence should be judged by seeing whether all elements of the gradient vector have absolute value less
than tol times the absolute value of the objective function plus fscale. Don’t forget to deal with the case
in which the Hessian is not positive definite1
. Also make sure that your step halving procedure deals with the
possibility of a trial step leading to a non-finite objective value (not only an increased objective). One appropriate
set of test functions are these (minimum at 1,1). . .


You should also design your own tests (of different dimensions for theta for example). The tests are not to
be submitted, but simply to check that the function does what is intended in general, and will pass marking tests.
