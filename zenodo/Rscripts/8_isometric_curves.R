# 15 random, increasing x coordinates bw 0 and 200 
a <- sort(sample.int(n=200, size=15))
  
# 7 random, increasing y coordinates bw 0 and 100
b <- sort(sample.int(n=100, size=7))
  
# 8 random, decreasing y coordinates
c <- c(sort(sample.int(n=100, size=8), decreasing=TRUE), b)
  
coords <- cbind(a,c)

# inspect curve
plot(coords)


# isometrically scale curve 
scaledcoords <- vector("list", 100)
for (i in 1:100){
  
  scaledcoords[[i]] <- coords*i
  
}


# extract baselines
baselines_list <- 
  scaledcoords %>% 
  lapply(., function(b) c(unlist(b[,1])[1], unlist(b[,1])[15]))


# fit splines and compute curvature
curvature_tbl <-
  mapply(curvature_spline, scaledcoords, baselines_list) %>% 
  enframe() %>% 
  mutate(total_K = abs(value)*(180/pi)) 


# estimate arclength as perimeter
perim_list <- Momocs::coo_perim(Ldk(scaledcoords))


# estimate chord length
chordlength <- Momocs::coo_length(Ldk(scaledcoords))


# estimate centroid size
csize <- coo_centsize(Ldk(scaledcoords))

alldata <- 
  perim_list %>% 
  unlist() %>% 
  enframe() %>% 
  mutate(totalK = curvature_tbl$total_K,
         chord = chordlength,
         size = csize) %>% 
  rename(perimeter = value) %>% 
  mutate(acratio = perimeter/chord)

plot(alldata$perimeter, alldata$chord)
plot(alldata$perimeter, alldata$totalK)
plot(alldata$chord, alldata$totalK)
plot(alldata$acratio, alldata$size) #all have the exact same ACR

