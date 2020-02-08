#run trace(fun2form, edit=TRUE) and change width.cutoff to 500L under deparse()

#total curvature function
totalK_fun<-function (x_range, fun, subdiv) 
{
  stopifnot(is.atomic(x_range)) # is.atomic checks that x.range cannot be a list or expression
  if (!is.numeric(x_range)) 
    stop("'x_range' must be a numeric vector!")
  if (length(x_range) != 2) 
    stop("'x_range' must be a vector of length two!")
  if (diff(x_range) < 0) #calculates the difference between x1 and x2 to ensure it's >0
    stop("please reorder 'x_range'.")
  if (!inherits(fun, "function")) 
    stop("'fun' must be a 'function' of x!")
  
  dfun <- deriv3(fun2form(fun), "x", func = TRUE)
  
  if (attr(dfun(x_range[1]), "gradient") == attr(dfun(x_range[2]), #make sure function is not a straight line
                                                 "gradient")) 
    stop("'fun' should not be a linear function of x!") 
  
  
  
  #gather coeffs to create expressions for deriv3() and functions for arclength()
  coeffs_lst <- 
    sapply(fun, function(v) v[1])
  
  #convert fun to param_fun
  as_param <- 
    function(x) x %>%
    as.numeric() %>%
    as.mpoly() %>% #creates a polynomial object of class 'mpoly'
    print() %>% 
    c("x", .) %>% #parametrize polynomial as t=(x, ax^3+bx^2+cx+d)
    mp() %>% #mpoly
    as.function()
  
  param_fun<-as_param(coeffs_lst)
  
  iter<- seq(0, 1, by=1/subdiv) #create vector of subdivisions to calculate arclength parameter
  arcfct_lst<- list() #empty bin
  b<- arclength(param_fun, x_range[1], x_range[2])$length #arc length of t-parameterized function
  
  for(i in seq_along(iter)){ 
    arcfct_lst[[i]] <- 
      local({
        b_sub<-iter[i]*b
        function(u) arclength(param_fun, x_range[1], u)$length - b_sub
      }) 
  }
  
  root_find<- function(x) uniroot(x, x_range)$root
  
  x <- sapply(arcfct_lst, root_find) #find roots for a list of 
  y <- fun(x) 
  gr <- attr(dfun(x), "gradient") #the tangents (first derv) of the x_n components, dfun() is defined above. The gradient matrix has elements that are the first deriv of a function
  he <- attr(dfun(x), "hessian")[, , "x"] # x is in the third dimension of this object (df?). The hessian matrix has elements that are the second deriv of a function
  k <- abs(he)/(1 + gr^2)^(3/2) #has n=subdiv measurements of k
  k_total<-(sum(k) /subdiv) *(180/pi) #add all measurements of k, rescale depending on #of subdivisions, and convert from rad to degrees
}




#add assign("dx_list", dX, envir = .GlobalEnv) to trace(arclength, edit=TRUE) to inspect objects created during arclength() calculations
# b<- arclength(param_fun, x_range[1], x_range[2])$length
