TODO:
Sit down with JS to review curvature function
Get results and write  about them
format citations
import figures

# _The inside curve: geometry and pollination ecology of curved flowers._

#### Abstract

The curvature of a flower along the lateral plane (‘floral curvature’) is a widespread, convergent trait with important ecological and evolutionary implications. This review summarizes the methods used to measure floral curvature, suggests a clarification of the definition of curvature, and its translation into a field-portable methodology. Intuitively, curvature is the degree to which a line is not straight. In plane geometry curvature is defined as the rate at which the unit derivative changes direction with respect to arc length. To apply this definition we suggest a protocol wherein a line is regressed against landmarks placed on a lateral image of an organism, then computing curvature at many points along the fitted line and taking the sum. The utility of this metric was tested by studying the development of nectar spur curvature in _Epimedium_ (Berberidaceae). Differences in the development of floral curvature were detected between _Epimedium koreanum_ and _Epimedium grandiflorum_ var. _violaceum_. Inflection points are found in wild-type _Epimedium grandiflorum_, but are lost in the cultivated variety ‘_violaceum_’, resulting in loss of total curvature. The suite of functions used to quantify floral curvature in this study are available as an open-source R package ‘curvy’.  The major advantages of this method are 1) precision of measurement is increased without introducing expensive field equipment or computing power, 2) precision of terminology within pollination ecology is improved by borrowing from the existing mathematical literature for studying line-curves, and 3)  the opportunity is opened for investigating the genetic basis of (lateral plane) curvature measured at the cellular scale.

#### 1. The ecology of floral curvature

> “We are beginning to understand why some hummingbird bills are long, whereas others are short, and why some hummingbird flowers are wide, whereas others are narrow. Now, why are bills of some hummingbirds and the tubes of the flowers they visit curved?” – _Temeles_ (1996).

At the center of plant-pollinator diversification is a remarkable variety of floral form.  The notion that flower shape is under selection to reduce interspecific mating ('floral isolation', Grant 1949), points to the importance of morphological diversity in initiating and reinforcing reproductive isolation (Muchhala and Armbruster 2009). For example, in the rapid radiation of Andean _Centropogon_ (Campanulaceae) , competition for pollination led to the divergence of floral traits associated with bat and hummingbird pollination (Lagomarsino and Muchhala 2019). Meanwhile, in South African _Lapeirousia_ (Iridaceae), geographic variation in floral tube length has initiated reproductive isolation between morphs with short and long flower tubes, despite sharing the same fly pollinator (Minnaar et al 2019). Thus, floral morphology is a crucial phenotypic feature underlying the diversification of plants and pollinators (Kay and Sargent 2009, Ollerton 2017). 

Flower-pollinator curvature as viewed from the side (lateral plane), has been a trait of special interest since the post-Darwin era of pollination ecology. In making pollinator observations of the Cape flora, Scott-Elliott (1890) noticed that the flowers of _Leonotis ocymifolia_ (Lamiaceae) visited by _Nectarinia_ sunbirds were "curved with the same curvature as that of the bird's beak." (p. 272). Robertson (1889) insightfully notes that the curved nectar spurs of _Viola_ spp. (Violaceae) "serves to limit the insect visits much more than the mere length of the spur." (p. 172). From these early observations curvature has been synonymous with specialization; we expect curvature to limit the range of functional taxa in a plant-pollinator mutualism (Scott-Elliot 1890) and strengthen interactions between the existing participants (Robertson 1889). And these expectations have largely been supported: Stiles (1975) first posited that neotropical _Heliconia_ partition hummingbird visitation by flower-bill curvature, and that specialization by curve-billed hummingbirds allow co-existence within the species-rich _Heliconia_ clade.  Subsequent research supports this hypothesis (Maglianesi et al 2014):  along the slopes of the Central Cordillera (Costa Rica), the degree of flower-hummingbird bill curvature is proportional to plant-pollinator interaction strength (Dehling et al 2014) and extent of specialization (_d'_, Bluthgen et al 2006). More recently, the spatial scope of plant-pollinator research has expanded to address the biogeography of curvature. As predicted by Stiles (2004), Maglianesi et al (2015) and Sonne et al (2019) find curvature to be most prevalent in the lowland environments of the neotropics. Explanations for this pattern range from heightened competition at lower elevations to environmental filtering in the Andean highlands (Stiles 2004, Graham et al 2009). Because the neotropical hummingbird subfamily Phaethornithinae comprises the majority of species with curved bills, we might expect plant-hummingbird curvature to have a predictable global distribution. 

Pollinator specialization has profound effects on macroevolutionary and biogeographic patterns (Kay and Sargent 2009, Armbruster and Muchhala 2009, Vamosi et al 2018), and curvature is a component, but widespread feature of specialist systems. Therefore,  to synthesize our knowledge of curved plant-pollinator systems, curvature is a concept that needs an exact definition and method of measurement. In this review we summarize the approaches to measuring curvature within the field of plant-hummingbird pollination, identify strengths and shortcomings, and offer a solution with the aim of standardizing how curvature is studied within the field of pollination ecology. Although this review is motivated by the problem of measuring curvature in plant-hummingbird systems, the solution is general to any biological form modelled as a line curve: this case is hopefully made in the demonstration to follow. 

#### 2. Summary of the literature 

We searched the scientific literature for studies focusing on or considering floral curvature - a common measurement used in pollination ecology as a proxy for specialization. We make the distinction between measuring curvature in the lateral plane versus measuring curvature of surfaces (e.g. petals). The methods for handling surface (Gaussian) curvature are relatively developed and well-defined (Nath et al 2003, Coen and Rebocho 2016). 

The literature was sourced by querying Google Scholar (Also do second, up-to-date sweep on WOS and litsearchr) for the terms "corolla curvature", "flower curvature", and "hummingbird curvature".  Results are summarized in Table 1. 

25 studies of hummingbird-plant systems were found using some form of curvature metric (Table 1).  An additional 10 studies of other plant, insect, and passerine bird systems were found (Table S1). There were numerous studies on plant-pollinator shape without a component on curvature that were omitted. 

The discussion of lateral curvature in plant-hummingbird interactions begins with Hainsworth (1973) and is first empirically studied by Stiles (1975), though methods for measuring curvature of bird bills outside of a pollination context can be found much earlier (Baldwin 1931). We identified six common approaches to measuring curvature. First, there are qualitative descriptions, e.g. "very curved", "less curved", but these have generally been out of use since the 1970s. Second, the _arc:chord_ method wherein curvature is a ratio of two lines: a straight line (chord) from tip to base (of the flower or bill) and a line that traverses a path along the arc of the flower/bill (Figure 3). Third, the _orthogonal length_ method which defines curvature as a ratio of two lines: a straight line from base to tip and a perpendicular line that measures the width of the flower/bill. This method is another form of the _arc:chord_ method because for a given chord length, the length of the perpendicular line will be proportional to the arc length.  Fourth, the _angle of deflection_ method which considers curvature as the angle between the base of the flower/bill and its tip. This is another form of the _inverse radius_ method which approximates the entire length of the flower/bill as a segment of a circle. These methods are interchangeable because the radius of a circle can be calculated from the length and angle of a line that passes through it (Bell 1956, Temeles 2009, see: Figure S2).  Sixth, geometric morphometrics, which quantifies shape as a configuration of homologous points (landmarks) existing on a coordinate plane. 

The strength of methods 2 and 3 are their portability and accessibility. These measurements can be taken in the field, or soon after from photographs. The methods are intuitive and in the simplest case, require only a ruler, string, and protractor. However, even if the measurements are made using imaging software these methods have some shortcomings (discussed in Berns and Adams 2010) - principally, that there are many shapes that could produce the same curvature value. For the _inverse radius_ method, a curve is approximated with the segment of a circle. This method is insufficient for any flower and bill shapes that deviate from having constant curvature (e.g. nectar spurs of _Delphinium_ ). Similarly, the _angle of deflection_ is not sensitive to local features along the length of the flower/bill - only the start and end points are considered in the calculation.  

Starting in 2010, geometric morphometrics (GM) emerges in the pollination literature. GM has since steadily gained in popularity in this field due to its statistical rigour, reproducibility, and the intuitive visual results diagramed as transformations of shape (e.g. illustrations flower shape variation between samples). We briefly outline the reasoning of a GM protocol to introduce relevant concepts, but recommend the authoritative and lucid introduction to the topic written by Webster and Sheets (2010). A GM protocol for a 2-D object begins by placing the specimens on an _xy_ grid and assigning landmarks to locations on the specimen that are topologically or biologically homologous - a landmark is defined so that its location can be reproduced within and between samples. The set of landmarks representing the shape of an organism is a 'landmark configuration'.  In a comparative study, the samples must first be optimally aligned and their shape information isolated from their orientation, location, and size. This is done using a least-squares type protocol, most commonly the Generalized Procrustes Analysis (GPA). GPA-adjusted landmark configurations hereafter exist in a multidimensional shape space defined by the number of landmarks used and whether the landmarks were assigned to 2-D or 3-D specimens. Each landmark configuration contains unique information about the specimen's shape, and as such, occupies a unique position in the corresponding shape space.  These configurations are then "projected" onto a simpler Euclidian space, similar to the reduction of a spherical Earth onto a two-dimensional map (Webster and Sheets 2010). From here, familiar statistical procedures (e.g. PCA) can be performed to quantify variation in landmark configurations (shape) between samples. 

This is giant leap forward for morphological studies because GM is a complete protocol for measuring, quantifying, and comparing shapes with high precision and accuracy, as well as the covariation of these shapes with ecological variables of interest. The limitation of GM in explicitly analyzing curvature is that this method is concerned with quantifying configurations of landmarks, the entirety of a shape summarized as a set of landmarks.  Once the specimen has been reduced to a landmark configuration, it exists as a point in shape space. Parsing segments of landmark configurations for separate analyses (e.g. for curvature) is not currently part of the geometric morphometrics toolkit. Therefore, studies that have used this technique to analyse biological forms are able to compare shapes in their entirety, but are ultimately limited to making descriptive statements about how segments of shapes appear to have different curvatures (e.g. Berns and Adams 2013, p. 251). 

#### 3. What is curvature?

Reviewing the literature leads us to ask, “what is curvature?”. In fact there exists numerous, but related mathematical definitions.  The formalization of this concept is born from several independent threads through history (summarized in @bardini_2016)

Fortunately, this question has been addressed by geometers, and we paraphrase the better explanations of Coolidge (1952), Casey (1996), Rutter (2000), and Jia (2018). First, we narrow the scope to plane geometry (Playfair 1846, Heath 1956). In this space, a straight line is traced when a particle (point) moves without changing direction.  A curve is traced when a particle moves from its origin and does not necessarily maintain its original heading as it moves. At every step along its path, our point can be said to have a measure of direction and speed.  The degree to which the particle’s direction changes from one instant to the next is its curvature (Figure 6). Intuitively, a particle that traces a straight line will not change its direction from one position to the next and will have zero curvature. A particle that traces a turning line will have relatively greater curvature.  

-Plant physiologists talk about curvature in the context of growth in relation to gravity (graviception). In this field, curvature is discussed as the rate of change of the angle of a material element (e.g. a cell) with respect to its position (s). See: @bastien_2014. 

To compute the curvature of a line a(t)=(x(t), y(t)) at a given moment t, it is useful to parametrize the curve in terms of its arc length s. If the relationship between t and s is known (u=t(s)), then the curve can be re-parameterized as a(u)= (x(u), y(u)). This ensures that the particle will trace out the line in equal intervals of s (ds). At every point along s there is a unit tangent vector T=T/|T|indicating the direction of the particle traveling the curve. Curvature (K) is the change in the value of T, dT, from s to ds . K = dT/ds. Jia (2018) defines curvature as ‘the change in the direction of the tangent with respect to arc length’. This definition has utility in an eco-morphological context, and we suggest that the term 'curvature' be reserved for this meaning. 

Total curvature (Rutter 2000)

#### 4. A proposed protocol for measuring curvature

In order to apply the above definition of curvature, a biological organ or tissue needs to be reduced to a continuous function.  We propose a workflow as illustrated in Figure 7. Cosgrove (1990) uses an analogous approach to study the development of cucumber hypocotyls. By fitting cubic splines to hand-marked seedlings, curvature was computed using the same definition as above.  Since Cosgrove (1990), the entire field of landmark-based geometric morphometrics has unfolded  (Bookstein 1990, Adams Rohlf Slice 2013). This rigorous, reproducible toolkit has been used extensively in pollination ecology, but has not yet been leveraged to calculate curvature (Table 1).  Terral et al (2004) use these tools to digitally landmark olive stones and fit polynomials to the landmarks: synthesizing the concepts of Cosgrove (1990) and Terral et al (2004) produces a modernized method for fitting curves and computing curvature from biological forms (Figure 7). 

Takahashi (2006) uses polynomials to estimate pelvis curvature in humans.

Comment on why using the simplest case (polynomials) and not splines, fourier, etc 

#### 5. Proof of concept: A study of the development of curvature in _Epimedium_

We tested the utility of this curvature metric by studying floral development in _Epimedium grandiflorum_ C.Morren and _Epimedium koreanum_ Nakai (Berberidaceae, Table 2 - sample sizes). Flower size was measured daily from April 9 to May 2, 2019 at the UBC Botanical Garden (Supp Mat 1). By correlating changes in flower size to developmental landmarks (Supp Fig 1), we were able to define 7 discrete stages of flower development (Table 2, Figure 5 - photographs of the stages). Flowers were sampled haphazardly and preserved in 70% ethanol. Preserved flowers were later transferred to a glass slide and imaged in the lateral view using a stereo microscope at 0.63x (Zeiss Stemi 508 with Axiocam 301). Specimens that did not fit within the field of view were imaged 2x or 3x and the images joined using the Stitching Plugin in the Fiji distribution of ImageJ2 (Preibisch et al 2009, Rueden et al 2017).

Photographed specimens were landmarked digitally using tpsDig [@rolhf_2015]. We placed landmarks along the edge of dorsal petals (in lateral view) as an approximation of the flower's total shape (see discussion of geometric morphometrics above).  Landmarks used to measure the dorsal arc were 1) the farthest point on the apex of the spur before the inflection point where either the spur diminishes to a tip (_E. violaceum_) or widens into a saccate reservoir (_E. koreanum_), and 2) the inflection point at which the spur widens to become an attachment for the petal to the stem(?). 13 semi-landmarks were placed between landmarks 1 and 2 (illustrated in Figure S3).

Landmark files (.tps) were imported into R using the _Momocs_ package v.1.3.0 [@bonhomme_2014]. Polynomial functions were regressed to the landmark coordinates for each specimen using _Momocs_- we chose polynomials of the third degree based on the recommendations of @rohlf_1990. Arc length was calculated from bounded polynomial functions using the _pracma_ package v.2.2.5 [@borchers_2019]. Curvature, as defined in the previous section, was computed using a custom function modified from the _soilphysics_ package v.3.1 [@silva_2017]. All custom functions used in this analysis are available as an R package _curvy_ - install by devtools::install_github("mannfred/curvy"). R scripts used in this analysis are hosted at github.com/mannfred/epidmedium.




Results: Curvature in 


Supp Mat 1:
#### Measuring flowers
In E. grandiflorum, the distance separating the outer sepals (sensu Stearn 2002) was measured using dial calipers (graduation = 0.1 mm), until the length of the inner sepals exceeded the length of the outer sepals (stage 4, Table 2 - describing the stages). From this point onwards, the inner sepal distance was measured. Because the aestivation was imbricate, we measured the sepals of the major axis (Figure demonstrating measurement technique). In _E. koreanum_, the inner sepals lack pigmentation and adhere closely to the petals, making them difficult to measure accurately in situ. For this reason, the outer sepals were measured until they abscised (stage 5, table 2). Flowers were sampled opportunistically and preserved in 70% ethanol. 

Supp Mat 2:
#### Landmarking
1. Rotate the photographs so that the opening of the corolla tube is parallel to the y-axis. 
2. Build tps file (a file listing all specimens) using tpsUtil. This tps file is used by tpsdig to add landmarks to.
3. Landmark specimens from tps file using tpsDig (steps 1 and 2 will soon be possible in MomX and could be done in geomorph). Landmarks used to measure the dorsal arc are 1) the farthest point on the apex of the spur before the spur diminishes to a tip (E. violaceum) or widens into a nectar bucket (E. koreanum), and 2) the dorsal point at which the spur widens to become an attachment point for the petal to the stem(?). 13 semi-landmarks are placed between them (15 points total). 
4. Curve points are drawn in tpsDig using the “pencil tool”, from landmark 1 (see above) to landmark 2. Following the placement of points, a curve is drawn that connects them. Right-click the curve and select “Resample Curve” and then space the points evenly “by length”. Manually adjust re-sampled points onto specimen and again “resample curve”. This does not usually need to be repeated more than twice.
5.1 Set scale by going to Options->image tools and typing in desired length and units. Press ‘set scale’ and then click on both ends of the scale bar in your image. Then go back to the image options box and select ‘OK’. 
5.2 Semi-landmarks need to be treated like landmarks for curve-fitting. To do this, use the ‘Append tps curve to landmarks’ function in tpsUtil. https://www.researchgate.net/post/Which_software_should_I_use_for_placing_sliding_semilandmarks 
6.  Import into R using  from_tps() function from Momit. 
7. Superimpose the shape using the fgProcrustes() function in Momocs.
8. Calculate polynomials and plot curves following the “Olea” example at: https://momx.github.io/Momocs/reference/opoly.html 
How do I decide what degree polynomial to use? Consult: Rohlf 1990…





<br>

```{r include=FALSE} 
knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE)
```

<br>

 You can also embed plots, for example:

```{r curve_example, results='hide', fig.cap="Figure 1. Very cool figure that I done made."}
<insert r code here>
```

<br>

```{r tables_lit_review}
table_1<-read.csv('Table_1_review_lit_curvature.csv', header=TRUE)
knitr::kable(table_1, caption = 'Table 1. Summary of literature reviewed for the role of curvature in plant-pollinator systems.')

table_2<-read.csv('Table_2_epimedium_flower_stages.csv', header=TRUE)
knitr::kable(table_2, caption = 'Table 2. Description of flower stages for Epimedium.')

```