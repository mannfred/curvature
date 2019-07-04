#i would like an R function that calculates total curvature...
#totcurv(f, a, b) 
#maybe: package(seewave, 'roughness')

#pracma - arclength

function (f, a, b, nmax = 20, tol = 1e-05, ...) #f is a parameterized function, a and b are endpoints, nmax = #iterations,  
{
  stopifnot(is.numeric(a), length(a) == 1, is.numeric(b), 
            length(b) == 1)
  fun <- match.fun(f)
  f <- function(x) fun(x, ...)
  if (abs(b - a) < tol) 
    return(list(length = 0, niter = 0, rel.err = tol))
  fa <- f(a)
  fb <- f(b)
  m <- length(fa)
  if (length(fa) < 2) 
    stop("Argument 'f' must be a parametrized function.")
  if (length(f(c(a, b))) != 2 * m) 
    stop("Argument 'f' must be a vectorized function.")
  h <- (b - a)
  A <- matrix(0, nmax, nmax)
  A[1, 1] <- sqrt(sum((fb - fa)^2))
  for (i in 1:(nmax - 1)) {
    h <- h/2
    x <- seq(a, b, by = (b - a)/2^i)
    y <- c(f(x))
    X <- matrix(y, ncol = m)
    dX <- diff(X)
    A[i + 1, 1] <- sum(sqrt(rowSums(dX^2)))
    for (j in 1:i) {
      A[i + 1, j + 1] <- (4^j * A[i + 1, j] - A[i, j])/(4^j - 
                                                          1)
    }
    if (abs(A[i + 1, i + 1] - A[i, i]) < tol && i > 3) 
      break
  }
  e <- abs(A[i + 1, i + 1] - A[i, i])
  list(length = A[i + 1, i + 1], niter = i + 1, rel.err = e)
}


#soilphysics - maxcurv

function (x.range, fun, method = c("general", "pd", "LRP", "spline"), 
          x0ini = NULL, graph = TRUE, ...) 
{
  stopifnot(is.atomic(x.range))
  if (!is.numeric(x.range)) 
    stop("'x.range' must be a numeric vector!")
  if (length(x.range) != 2) 
    stop("'x.range' must be a vector of length two!")
  if (diff(x.range) < 0) 
    stop("please, reorder 'x.range'.")
  if (!inherits(fun, "function")) 
    stop("'fun' must be a 'function' of x!")
  method <- match.arg(method) #indicates that the method parameter must match one of "general", "pd", etc.
  dfun <- deriv3(fun2form(fun), "x", func = TRUE)
  if (attr(dfun(x.range[1]), "gradient") == attr(dfun(x.range[2]), #make sure function is not a straight line
                                                 "gradient")) 
    stop("'fun' should not be a linar function of x!")
  x <- seq(x.range[1], x.range[2], length.out = 5000) #splits the x range into 5000 even pieces
  y <- fun(x) 
  if (method == "general") {
    gr <- attr(dfun(x), "gradient") #the tangents (first derv) of the 5000 x components
    he <- attr(dfun(x), "hessian")[, , "x"] 
    k <- abs(he)/(1 + gr^2)^(3/2) #there are possibly 5000 curvature points stored in k... could sum them?
    mcp <- x[which.max(k)]
  }
  else if (method == "pd") {
    b <- lm(range(fun(x.range)) ~ x.range)$coef
    if (fun(x.range[1]) > fun(x.range[2])) {
      si <- -1
    }
    else {
      si <- 1
    }
    b1 <- si * b[2]
    b0 <- mean(fun(x.range)) - b1 * mean(x.range)
    ang <- atan(b1)
    a1 <- -1/b1
    a0 <- y - a1 * x
    xjs <- (b0 - a0)/(a1 - b1)
    yjs <- b0 + b1 * xjs
    k <- sqrt((x - xjs)^2 + (b0 + b1 * xjs - y)^2)
    mcp <- x[which.max(k)]
  }
  else if (method == "LRP") {
    if (is.null(x0ini)) 
      stop("please, inform 'x0ini', a initial value for x0")
    ini <- coef(lm(y ~ x))
    fit <- try(nls(y ~ flrp(x, a0, a1, x0), start = list(a0 = ini[1], 
                                                         a1 = ini[2], x0 = x0ini)), silent = TRUE)
    if (class(fit) == "try-error") {
      fit <- try(nls(y ~ flrp(x, a0, a1, x0, left = FALSE), 
                     start = list(a0 = ini[1], a1 = ini[2], x0 = x0ini)), 
                 silent = TRUE)
    }
    if (class(fit) == "try-error") {
      stop("LRP could not get convergence!")
    }
    else {
      mcp <- coef(fit)[3]
      k <- rep(0, length(x))
    }
  }
  else {
    if (is.null(x0ini)) 
      stop("please, inform 'x0ini', a initial value for x0")
    lini <- coef(lm(y[x < x0ini] ~ x[x < x0ini]))
    rini <- coef(lm(y[x > x0ini] ~ x[x > x0ini]))
    fit <- try(nls(y ~ fs(x, a0, a1, b1, x0), start = list(a0 = lini[1], 
                                                           a1 = lini[2], b1 = rini[2], x0 = x0ini)), silent = TRUE)
    if (class(fit) == "try-error") {
      stop("spline could not get convergence!")
    }
    else {
      mcp <- coef(fit)[4]
      k <- rep(0, length(x))
    }
  }
  if (graph) {
    curve(fun, from = x.range[1], to = x.range[2], ...)
    lines(x = c(mcp, mcp, -9e+09), y = c(-9e+09, fun(mcp), 
                                         fun(mcp)), lty = 3)
    if (method == "pd") {
      abline(b0, b1, lty = 3)
      lines(x = c(mcp, xjs[which.max(k)]), y = c(fun(mcp), 
                                                 b0 + b1 * xjs[which.max(k)]), lty = 3)
    }
    else if (method == "LRP" || method == "spline") {
      lines(x, predict(fit), col = 4, lty = 2)
    }
    devAskNewPage(ask = TRUE)
    plot(x, k, type = "l", ...)
    lines(x = c(mcp, mcp, -9e+09), y = c(-9e+09, max(k), 
                                         max(k)), lty = 3)
  }
  out <- list(fun = fun, x0 = mcp, y0 = fun(mcp), method = method)
  class(out) <- "maxcurv"
  return(out)
}



#####total curvature?

#soilphysics - maxcurv


totalK<-function (x.range, fun) 
{
  stopifnot(is.atomic(x.range))
  if (!is.numeric(x.range)) 
    stop("'x.range' must be a numeric vector!")
  if (length(x.range) != 2) 
    stop("'x.range' must be a vector of length two!")
  if (diff(x.range) < 0) 
    stop("please, reorder 'x.range'.")
  if (!inherits(fun, "function")) 
    stop("'fun' must be a 'function' of x!")
  
  dfun <- deriv3(fun2form(fun), "x", func = TRUE)
  
  if (attr(dfun(x.range[1]), "gradient") == attr(dfun(x.range[2]), #make sure function is not a straight line
                                                 "gradient")) 
    stop("'fun' should not be a linear function of x!") #corrected spelling from "linar"
  x <- seq(x.range[1], x.range[2], length.out = 5000) #splits the x range into 5000 even pieces
  y <- fun(x) 
  gr <- attr(dfun(x), "gradient") #the tangents (first derv) of the 5000 x components, dfun() is defined 7 lines above. The gradient matrix has elements that are the first deriv of a function
  he <- attr(dfun(x), "hessian")[, , "x"] # x is in the third dimension of this object (df?). The hessian matrix has elements that are the second deriv of a function
  k <- abs(he)/(1 + gr^2)^(3/2) #there are possibly 5000 curvature points stored in k... could sum them? also, see: https://en.wikipedia.org/wiki/Curvature "curvature of the graph of a function"
}


##test
library(soilphysics)

fun=function(x) 0.5*x^2

curvature<-totalK(x.range=c(0,1), fun=fun)

total_curvature<-sum(curvature[c(1:5000),]) 

#4472.33 is the total curvature for y=x^2 [0,1]
#4851.77 is the total curvature for y=2*x^2 [0,1]
#3535.50 is the total curvature for y=0.5*^2 [0,1]

#to do:

#sanity check the totalcurvature computation

#divide total curvature by arc length
#write a pipe that can 
#1. fit polynomials to a list of semi-landmarks from many epimedium individuals
#2. calculate arclength and total curvature/arc length for a list of polynomials
#3. $$$profit

#translate maxcurv() output into arc length by
#1. parameterizing polynomials by arclength (pracma)
#2. taking (x,y) coords from maxcurv() and putting them into the arclength parameterized polynomial eq. 
#3. adding the arc length value to the maxcurv() graph?? maybe just in text in the bottom-right corner?

#fun2form
function (fun, y = NULL) 
{
  .doTrace({
  }, "on entry")
  {
    if (!inherits(fun, "function")) 
      stop("'fun' must be an object of class 'function'!")
    if (!is.null(y) & !inherits(y, "character")) 
      stop("'y' must be a 'character' which is going to define the left side of the formula!")
    form <- as.formula(paste(y, "~", deparse(fun)[2]))
    return(form)
  }
}