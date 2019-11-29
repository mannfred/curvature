#######################
Vector-valued functions
#######################

The curves found in nature are sometimes best described by multi-valued fuctions, i.e. functions that allow the curve to have several Y-values for a single X-position (e.g. spirals). For this reason, the standard y=f(x) notation will not suffice, and instead we introduce the concept of a vector-valued function. 

We first introduce the concept of a vector-valued function. Physically, a vector can be thought of as an arrow-like object with length and direction [Figure 4](Figures/Figure_4.jpg). A vector that exists in a coordinate plane "points" to an $(x, y)$ location and is called a _position vector_.  The $x$ and $y$ components of a position vector can be handled separately as scaled _unit vectors_ of length 1.  Therefore, the length of a position vector is computed as the Euclidian distance (hypotenuse) between its rectangular $x$ and $y$  unit vector components (Figure 4). A vector's direction is determined by its relative orientation to an arbitrarily chosen fixed line (_e.g._ the $x$-axis). Because vectors can "point" to a location on a coordinate plane, an $(x, y)$ coordinate can 
can be rewritten as a position vector <center> $\mathbf{r} = x\mathbf{i} + y\mathbf{j}$,  </center> 

Where $\mathbf{i}$  and $\mathbf{j}$ are orthogonal unit vectors fixed on the $x$ and $y$ axes, respectively, and $\mathbf{r}$ is the resultant position vector pointing to some point $(x, y)$ ([Figure 4](Figures/Figure_4.jpg)).  

Now, instead of the typical $f(x) = y$ notation for a function, consider a curve as the path traced by a moving particle that passes through a set of $(x, y)$ points, described by two time-dependent functions <center> $x = \hat{x}(t)$, <br> $y = \hat{y}(t)$, </center>

where $\hat{x}$ and $\hat{y}$ denote that $x$ and $y$ are functions of time $t$.  This parameterization is convenient because we would like to treat position vectors as composites of the unit vectors $\mathbf{i}$ and $\mathbf{j}$.  Consequently, we can combine the concepts of eqs 1 and 2 to produce a _vector-valued function_  <center> $\mathbf{\hat{r}}(t) = \hat{x}(t)\mathbf{i} + \hat{y}(t)\mathbf{j}$ ,  </center> 

where $\mathbf{\hat{r}}(t)$ is a function that produces a position vector $\mathbf{r} = x\mathbf{i} + y\mathbf{j}$ determined by $t$. 

$$
\begin{align}
\mathbf{\bar{r}}(s) = \bar{x}(s)\mathbf{i} + \bar{y}(s)\mathbf{j} && \text{(eq. 1)}
\end{align}
$$

Where $\mathbf{\bar{r}}(s)$ indicates that the position vector $\mathbf{r}$ is determined by its arc length $s$. 
When a curve is parameterized by its arc length, for every arc length increment ($ds$) that we move along the curve, we generate a position vector $\mathbf{r}$ indicating the $(x, y)$ location in the coordinate plane. 


This is because as ${\Delta s \to 0}$, the length of the tangent vector ($d\mathbf{r}$) and the length of the arc segment ($ds$) become equal. In this case, the tangent vector is called the _unit tangent vector_. 
#################
Measuring flowers
#################

In E. grandiflorum, the distance separating the outer sepals (sensu Stearn 2002) was measured using dial calipers (graduation = 0.1 mm), until the length of the inner sepals exceeded the length of the outer sepals (stage 4, Table 2 - describing the stages). From this point onwards, the inner sepal distance was measured. Because the aestivation was imbricate, we measured the sepals of the major axis (Figure demonstrating measurement technique). In _E. koreanum_, the inner sepals lack pigmentation and adhere closely to the petals, making them difficult to measure accurately in situ. For this reason, the outer sepals were measured until they abscised (stage 5, table 2). Flowers were sampled opportunistically and preserved in 70% ethanol. 

############
Landmarking
############

1. Rotate the photographs so that the opening of the corolla tube is parallel to the y-axis. 
2. Build tps file (a file listing all specimens) using tpsUtil. This tps file is used by tpsdig to add landmarks to.
3. Landmark specimens from tps file using tpsDig (steps 1 and 2 will soon be possible in MomX and could be done in geomorph). Landmarks used to measure the dorsal arc are 1) the farthest point on the apex of the spur before the spur diminishes to a tip (E. violaceum) or widens into a nectar bucket (E. koreanum), and 2) the dorsal point at which the spur widens to become an attachment point for the petal to the stem(?). 13 semi-landmarks are placed between them (15 points total). 
4. Curve points are drawn in tpsDig using the “pencil tool”, from landmark 1 (see above) to landmark 2. Following the placement of points, a curve is drawn that connects them. Right-click the curve and select “Resample Curve” and then space the points evenly “by length”. Manually adjust re-sampled points onto specimen and again “resample curve”. This does not usually need to be repeated more than twice.
5.1 Set scale by going to Options->image tools and typing in desired length and units. Press ‘set scale’ and then click on both ends of the scale bar in your image. Then go back to the image options box and select ‘OK’. 
5.2 Semi-landmarks need to be treated like landmarks for curve-fitting. To do this, use the ‘Append tps curve to landmarks’ function in tpsUtil. https://www.researchgate.net/post/Which_software_should_I_use_for_placing_sliding_semilandmarks 
6.  Import into R using  from_tps() function from Momit. 
7. Superimpose the shape using the fgProcrustes() function in Momocs.
8. Calculate polynomials and plot curves following the “Olea” example at: https://momx.github.io/Momocs/reference/opoly.html 