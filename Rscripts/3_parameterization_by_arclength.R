f <- function(t) c(t, 7.024797 + 1.840848*t - 0.09349572*t^2) #  7.024797 + 1.840848*t - 0.09349572*t^2 ,  baselines 6.591728  21.259368 

t1 <- 6.591728 
t2 <- 21.259368
a  <-  arclength(f, 0, t1)$length
b  <-  arclength(f, t1, t2)$length

fParam <- function(w) {
  fct <- function(u) arclength(f, t1, u)$length - w
  urt <- uniroot(fct,  c(6.591728, 21.259368))
  urt$root 
}

ts <- linspace(0, 25, 250)
plot(matrix(f(ts), ncol=2), type='l', col="blue", 
     asp=1, xlab="", ylab = "",
     main = "marmalade!", sub="20 subparts of equal length")

#plotting 

for (i in seq(0, 1, by=0.05)) {
  v <- fParam(i*b); fv <- f(v)
  points(fv[1], f(v)[2], col="darkred", pch=20)
} 



# get the xy coords of plotted points
q<-capture.output(
  for (i in seq(0, 1, by=0.05)) {
    v <- fParam(i*b)
    cat(v,"\n")
  } ) %>% as.numeric()


#corresponding y values for uniroot outputs
q_fq<-
  lapply(q, f) %>% 
  transpose() %>%  
  simplify2array()

points(q_fq, col="darkgreen", pch=20) #using different colour to check that q_fq is equivalent to v (above)

coo_draw(dorsal_lst[1]) #plot epimedium curve #1, t values align. 



#check arclength distances
r<-vector("list", length=nrow(xy_coords))

for (i in 1:(nrow(q_fq)-1)) {
  r[[i]]<-
  arclength(f, q_fq[[i,1]], q_fq[[i+1,1]])$length }

#all arc length distances are equivalent
