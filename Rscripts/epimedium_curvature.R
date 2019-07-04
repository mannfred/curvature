########################################
#measure curvature from epimedium photos
########################################
library(here)
library(Momit)
library(Momocs)
library(geomorph)
library(tidyverse)
library(mpoly)
library(pracma)
library(soilphysics)
library(polynom)

#import dorsal data####
dorsal<- here("data/epimedium_photos/koreanum/koreanum_dorsal_appended.TPS") %>% #mounts tps file
         readland.tps(specID='imageID') %>% #geomorph func, auto-applies the scaling factor to LMs
         a2l() %>% # Momit func, converts array to list
         Ldk() #adds "$shp{i}" as list element headers
        

#inspect raw LMs
str(dorsal)#LMs stored in $coo
dorsal %>% paper %>% draw_curve #draws all curves on one plot (not yet procrustes superimposed)
rapply(dorsal, coo_plot) #recursively plots all curves
coo_plot(dorsal[i]) #for i={1:31} to plot individually

#calculate polynomials
dorsal_curv <- npoly(dorsal$coo, degree=3) #31 polynomial curves estimated 
str(dorsal_curv[[1]]) #"Call:" = function+parameters  used to create the model

#draw polynomials estimates over raw LMs
dorsal_curv[[i]] %>% #for i={1:31} -- LMs must be previously plotted w coo_plot() 
npoly_i() %>%
coo_draw()

#save plot
dev.copy(figure,"Figure_1.png",width=8,height=6,units="in",res=100)
dev.off()


#calculate arc length####


#create polynomial functions from coeffs

coeffs<-sapply(dorsal_curv, function(v) v[1]) 
#extracts the first element ($coeff) from each element ($shp) in the list. 
#sapply stores results as a list of atomic vectors

aspolyfunc<-function(x) x %>%
  as.numeric() %>%
  as.mpoly() %>% #creates a polynomial object of class 'mpoly'
  print() %>% 
  c("x", .) %>% #parametrize polynomial as t=(x, ax^3+bx^2+cx+d)
  mp() %>% #mpoly
  as.function()

listpolyfunc<-lapply(coeffs, aspolyfunc) # a list of 31 parameterized polynomial functions


#extract the x-coords from the "baseline" entry for every $shp
#max/min baselines are stored in [4] and [5], [[1]][1] extracts only the x-coord
max_baselines<-(lapply(dorsal_curv, function(b) b[4][[1]][1]))
min_baselines<-(lapply(dorsal_curv, function(b) b[5][[1]][1]))

#calculates arclength for every polynomial bounded by baselines
lengths<-mapply(arclength, listpolyfunc, min_baselines, max_baselines) %>%
         as.tibble() %>% #row 1 contains arclengths 
         slice(., 1) %>% #keep only row 1
         as.list()






#calculate total curvature####


#extracts the lower xy-boundary from b[5], unlists it, and isolates the x-boundary by [1]
#then does the same for b[4] which is the upper xy-boundary stored in dorsal_curv$baseline1
baselines <- dorsal_curv %>% 
             lapply(., function(b) c( unlist(b[5])[1], unlist(b[4])[1]))


polylist<- coeffs %>%
           lapply(., polynomial) %>% #convert coeffs to a list of polynomials
           lapply(., as.character) #convert polynomials to character vectors
  

#convert polynomials stored as character strings to quoteless expressions
poly2func<- function (p) 
{
  f<- function(x) NULL
  body(f) <- parse(text=p)
  f
}


funclist<- polylist %>%
           lapply(., poly2func) #create a list of one-line polynomial functions


#run trace(fun2form, edit=TRUE) and change width.cutoff to 500L under deparse()

#total curvature function
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

curvaturelist <- totalK %>%
                 mapply(., baselines, funclist) %>%
                 as.tibble


totalcurvaturelist <- curvaturelist %>% 
                      rapply(., sum) %>%
                      as.list()





#procrustes superimposition####

dorsal_procrustes<-fgProcrustes(dorsal) 
dorsal_procrustes%>% paper %>% draw_curve #draws all procrustes-superimposed curves
#two are inverted?
rapply(dorsal_procrustes, coo_plot) 


#fit polynomials
curves<-opoly(dorsal_procrustes) #degree=5 by default
