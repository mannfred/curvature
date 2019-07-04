library(here)
library(Momit)
library(Momocs)

#####################
#import ventral data#
#####################

ventral<-from_tps("centropogon_ventral_appended.tps") 

#inspect 
str(ventral)
ventral$LM[1]

ventral_coo<-ventral$LM[[1]] #isolates coordinates as a matrix
coo_plot(ventral_coo)

ventral_op <- opoly(ventral_o, degree=4)
ventral_op

ventral_opi<-opoly_i(ventral_op)
coo_draw(ventral_opi)
lines(ventral_opi, col='red')

#####################
##import dorsal data#
#####################

#import dorsal data
dorsal<-from_tps("centropogon_dorsal_appended.tps") 

#inspect 
str(dorsal) #dorsal is an object containing 4 classes of subobjects: 'mom_df', 'tbl_df' (aka ), 'tbl', and 'data.frame'
dorsal$LM[1] #Centropogon_leucocarpus
dorsal$LM[2] #Centropogon_yungasensis

#procrustes superimposition on both species 
Opn(dorsal$LM) %>% 
  rescale(scaling_factor = 1/500) %>% 
  fgProcrustes() %>%
  coo_plot() 
  


coo_plot(dorsal_coo) 

dorsal_curve <- npoly(dorsal_coo, degree=4)
dorsal_curve

dorsal_shape<-npoly_i(dorsal_curve)
coo_draw(dorsal_shape)
lines(dorsal_shape, col='red')



#rescale

dorsal_rescaled<-rescale(dorsal_coo, scaling_factor=1/500) #1cm = 500pixels

#get equation for polynomial by calling "dorsal_curve" (= 0.08683742+0.47066380*x^1-1.67068169*x^2-0.47594076*x^3-3.77902817*x^4)
#next up: estimate curve parameters (max curvature, arclength, total curvature/arc length)

######################

library(pracma) #for calculating arc length

#parameterized wrt to t
f1_para<-function(t) c(t, 0.08683742+0.47066380*t^1-1.67068169*t^2-0.47594076*t^3-3.77902817*t^4)

dorsal_imp$shp2[,1] #find x coordinates for arc length boundary

arclength(f1_para, -0.305862254, 0.346909032) #0.8550248 'units' 

#################

#arc length parameterization

#r(t)
func<- function(t) c(t, t^2)

t1<-0
t2<-1

a<-0
b<-arclength(func, t1, t2)$length

#r(s)
#u and w are the bounds of the arc length function
#so fct is the length of func bw u and w 
func_para<- function(w) {
  fct<- function(u) arclength(func, a, u)$length - w 
  urt<- uniroot(fct, c(a, t2))
  urt$root
  
} 

#func_para(1) = 0.7639264 
#i.e. when s=1, t=0.7639264 which checks out with func(0.7639264)=1.00


ts <- linspace(0, 1, 250)
plot(matrix(func(ts), ncol=2), type='l', col="blue", 
     asp=1, xlab="", ylab = "",
     main = "poo", sub="whatever")

for (i in seq(0.05, 0.95, by=0.05)) {
  v <- func_para(i*b); fv <- func(v)
  points(fv[1], func(v)[2], col="darkred", pch=20)
} 

#################

library(soilphysics) #for estimating and plotting max curvature

#a small bug:

#maxcurv() uses fun2form() which uses deparse(). 
#deparse() has a default width.cutoff of 60 characters, which my polynomial curvature function exceeds.
#I ran  trace(fun2form, edit=TRUE), and add "width.cutoff=500L" as a parameter of deparse(), which allows maxcurv() to proceed.

#otherwise gives Error in parse(text = x, keep.source = FALSE) : 
#<text>:2:0: unexpected end of input (<text>:2:0 means "line 2, character 0")
#1:  ~ 0.08683742 + 0.4706638 * x^1 - 1.67068169 * x^2 - (0.47594076^3) - 
#  ^

f1<-function(x) 0.08683742+0.47066380*x^1-1.67068169*x^2-0.47594076*x^3-3.77902817*x^4

#bug fix
trace("fun2form", edit=TRUE) #set width.cutoff=500L under deparse()
maxcurv<-maxcurv(x.range = c(-1, 1), fun = f1)
maxcurv

#critical x:  0.1754351 
#critical y:  0.1118395

#maxcurv uses deriv3(), which computes derivatives of expressions, including second derivatives (hessian=TRUE).
#this is usefull because curvature (k) is a second derivative (?)

#calculate arc length from base to max K 

arclength(f1_para, -0.305862254, 0.1754351) #0.6300547 'units'

0.6300547/0.8550248 #point of max K divided by total arc length
#max K is found at 0.7368847 (73%) along the total length of the curve










