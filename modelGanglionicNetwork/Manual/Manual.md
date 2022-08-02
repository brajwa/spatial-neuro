Generative Modeling of the Enteric Nervous System Employing Point
Pattern Analysis and Graph Construction
================
Abida Sanjana Shemonti
August 2, 2022

# Introduction

We construct a generative model employing spatial point pattern analysis
and graph theoratic approaches to describe the architecture of the
enteric nervous system (ENS) in the colon that works on the data from
images of human and mouse tissue samples obtained through confocal
microscopy.

# Data

The directory **spatial-neuro/modelGanglionicNetwork/Data/** contains
the following subdirectories:

-   **TIF Images:** This subdirectory contains ENS images (*.tif*) of
    the tissue samples collected from human colon with 3 color channels.

-   **Branch Information (in um):** There are *.csv* files in this
    subdirectory. We use an open source image processing package Fiji
    (Schindelin et al. 2012) based on ImageJ2 for processing the ENS
    images and extracting data from them. The extracted information
    contain the 2D coordinates of the end vertices of the ENS network,
    the network branch (edge) lengths (actual and Euclidean) etc. All
    distances are measured in microns.

    We describe the image processing steps done in Fiji for segmenting
    the ENS images here. Note that, the parameter values mentioned here
    for different filters are subjected to change. We tried to tune them
    to the values that worked for our images mostly.

    -   We load the *.tif* image into Fiji and flatten it. *Images \>
        Overlay \> Flatten*.

    -   We set the appropriate pixel to micron scale in *Analyze \> Set
        Scale*.

    -   If the entire image is not ideal for analysis, we crop out the
        region of interest (ROI) and make it an 8-bit image. *Images \>
        Type \> 8-bit*.

    -   Despeckle. *Process \> Noise \> Despeckle*. Subtract Background
        with rolling ball radius 50 pixels. *Process \> Subtract
        Background…*

    -   Enhance local contrast with block size 199. *Process \> Enhance
        Local Contrast (CLAHE)*. This step can take a few moments.

    -   We apply the fast Fourier transformation (FFT) based bandpass
        filter with large structures down to 150 pixels and small
        structures up to 50 pixels. *Process \> FFT \> Bandpass Filter*.

    -   At this point, we need to include an additional plugin to Fiji
        called MorphoLibJ. You may need to download its *.jar* to the
        Fiji plugins folder. Check out the source on git
        (<https://github.com/ijpb/MorphoLibJ>). *Plugins \> MorphoLibJ
        \> Segmentation \> Morphological Segmentation*. We run the
        process with Input Image *Border Image* and Watershed
        Segmentation tolerance 20. This step can take a few moments. We
        create image with Watershed Lines. This gives us the segmented
        image of the ENS network.

    -   We set the settings to Black Background in *Process \> Binary \>
        Options…* and invert the image. *Edit \> Invert*.

    -   *Analyze \> Skeleton \> Analyze Skeleton* gives us a table named
        *Branch information* with all the extracted numerical
        information about the ENS network that we need. We save this
        table in *.csv* format.

    **Note:** Fiji processes the images considering the top-left corner
    as the origin. So, before using the extracted coordinates in R
    codes, you MUST negate the y-coordinates, otherwise you will get a
    vertically flipped image.

-   **Neuron Intensity Profiles:** This subdirectory contains images
    (*.jpg*) of the intensity profile of the ganlionic neurons. These
    images are for visualization purpose and are used at the end of the
    model generation to create an overall picture of the simulated
    ganglionic network. We create a binary image of the simulated
    ganglionic network and process it in Fiji to construct an intensity
    profile for the ganglionic neurons. We invert the binary image
    (white network on black background) and may need to crop as
    required. We apply 3 image processing actions.

    -   Dilation. *Process \> Binary \> Dilate*.

    -   Euclidean distance mapping. *Process \> Binary \> Distance Map*.

    -   Binarization. *Process \> Binary \> Make Binary*.

    The crucial thing is that we need to tune to a reasonable
    combination of these 3 actions to get our desired intensity profile.
    This process can be quite subjective to individuals and prior
    biological knowledge. In our experiments we used a large number of
    dilation (10 to 15 times), followed by distance mapping and
    binarization twice or thrice, and we finished with a contrast
    enhancement.

-   **Supporting Files:** Some additional information of the ENS
    samples, such as sample ID, sample location and condition, age group
    and sex, and the scaling of the *.tif* images are saved in this
    subdirectory.

# Codes

The directory **spatial-neuro/modelGanglionicNetwork/Codes/** contains
the R code files. We utilize the R packages *spatstat* and *igraph* for
this purpose, along with other necessary libraries. In our experiments
we used *Fiji v1.53h* and *R version 4.0.2* in RStudio. The R files can
be run using RStudio or Linux terminal.

## Analyzing ENS characteristics

**Code: AnalyzeGanglionicNetwork.R**

### Constructing suitable data structures

We construct a few suitable data structures to store the ENS network
information extracted from the images, that will be useful in the
subsequent stages of the code. The function **constructDataStruct**
takes as input the name of the ENS sample under consideration
(*sample_id*), the path to the parent directory of the codes (*parent*),
the path to the corresponding sample’s output directory
(*output_folder_path*) and the path to the *.csv* file that has the ENS
network information (*branch_info_path*), and returns a dataframe object
(*branch_all*), a point pattern object (*branch_ppp*), a point pattern
on a linear network object (*branch_lpp*) and a graph object (*g1*).

**Note:** The ENS network information are stores in *.csv* format and
the file names are preceded with the sample ID. Any other format or
naming convension will require modification in the codes.

### Summary statistics

To understand the spatial organization of the ganglionic point patterns,
summary statistics such as - G, F, J, K, L functions are quite useful.
The following function **summaryStat** demonstrates the summary
statistics analysis for a given ganglionic point pattern. The function
**summaryStat** takes as input the name of the ENS sample under
consideration (*sample_id*), the path to the parent directory of the
codes (*parent*), the path to the corresponding sample’s output
directory (*output_folder_path*) and the path to the *.csv* file that
has the ENS network information (*branch_info_path*), and plots the
corresponding statistics as output.

### Analyzing ENS ganglia

We fit the ganglionic point pattern to an inhomogeneous hardcore-Strauss
process. The function **analyzeGanglia** takes as input the name of the
ENS sample under consideration (*sample_id*), the path to the parent
directory of the codes (*parent*), the path to the corresponding
sample’s output directory (*output_folder_path*) and the path to the
*.csv* file that has the ENS network information (*branch_info_path*),
and returns the fitted model (*ganglia_model*) along with its
parameters: intensity (*Beta*), interaction parameter (*Gamma*),
interaction distance (*R*), hardcore distance (*H*) and point pattern
window (*window*).

The library function **spatstat::ppm** has an *interaction* parameter.
In our case that parameter is the desired model object
**StraussHard()**, and we need to provide it with an interaction
distance of our preference. In the case of hardcore-Strauss process,
interaction distance *R* and interaction parameter *Gamma* are irregular
parameters, which can be determined by biological knowledge or by
analyzing the profile Akaike information criterion (AIC) of the point
pattern for a range of irregular parameters (See (Baddeley, Rubak, and
Turner 2015) Chapter 13.6.3 for details).

### Analyzing ENS network

We analyze the ENS network to figure out its structural characteristics,
which is very crucial for the next step in modeling - generating ENS
network on top of a simulated point pattern instance. We focus on the
degree of the vertices, the angle (in degree) and the length of the
edges. The function **analyzeBranch** computes these characteristics
values for the given ENS network and constructs their kernel density
estimates. It takes as input the name of the ENS sample under
consideration (*sample_id*), the path to the parent directory of the
codes (*parent*), the path to the corresponding sample’s output
directory (*output_folder_path*) and the path to the *.csv* file that
has the ENS network information (*branch_info_path*), and returns a
dataframe object (*branch_all*) containing all the feature values for
each of the edges of the ENS network, three kernel density estimates
(*orgKDE_angle*, *orgKDE_length* and *orgKDE_both*) containing the
kernel density estimates of the edge angle, edge length and both of the
given ENS network respectively, meshedness, compactness and density of
the given ENS network (*meshedness*, *compactness*, *network_density*).

The function **analyzeBranch** uses a few helping functions that are
also included below.

-   The function **range01** simply 0-1 normalizes a given vector.

-   The function **calcAngle** computes the angle of an edge with the
    horizontal x axis in degree. The input parameter *x* has four
    components - the coordinates of the end points of the edge *x1*,
    *y1*, *x2*, *y2*.

-   The function **calcDist** computes the Euclidean distance between
    the two nodes of an edge. The input parameter *x* has four
    components - the coordinates of the end points of the edge *x1*,
    *y1*, *x2*, *y2*.

-   The function **computeKDE** computes the distribution of angle
    and/or length of the ENS network. It takes as input a matrix *data*
    with column 1 that has the angle and column 2 that has the length of
    the edges, and an integer *mode* that is either 1/2/3 (1 for angle,
    2 for length, 3 for both).

-   The function **computeEdgeWeight** assigns a weight to each edge of
    the ENS network based on the degree of its two end vertices.

## Generating ganglia centers

**Code:GenerateGangliaCenters.R**

While analyzing the ENS ganglia, we fitted the ganglia point point
pattern to a hardcore-Strauss process. Given that fitted model or the
parameters of the model, the following function
**generateGangliaCenters** can generate simulated realizations of ENS
ganglia centers. If the parameter *with_model=1*, the function generates
ganglia centers directly from the model. If *with_model=0*, it depends
on the value of the parameter *process_type* in generating ganglia
centers from either of basic hardcore, basic Strauss or hardcore-Strauss
process with their required parameters. Finally, the function returns
the simulated realization as a point pattern object.

The coordinates of the generated ganglia centers are saved as *.csv*
file in the corresponding sample’s output subdirectory in
**spatial-neuro/modelGanglionicNetwork/Outputs/**.

## Generating ENS network

**Code:GenerateNetwork.R**

We generate simulated realizations of the ENS network on top of the
realization of ENS ganglia centers that is simulated above maintaining
statistical similarity with the given ENS network properties. We choose
Dealunay triangulation as the initial deterministic connection model.
Then rejection sampling of network edges is incorporated as a random
connection model, based on the kernel density estimations of the ENS
network characteristics.

### Deterministic connection model

The function **deterministicEdges** takes as input the simulated
realization of the ENS ganglia centers (*ganglia_ppp*) and the dataframe
object containing the given ENS network characteristics (*branch_all*).
It constructs the Delaunay triangluation on *ganglia_ppp* and suitable
data structures for the next edge sampling step (*network_extra*). It
also computes the kernel density estimations of the edge angle, length
and both for the triangulation. The function returns all these computed
objects.

### Random edge sampling process

The function **rejectionSampling** takes as input the simulated
realization of the ENS ganglia centers (*ganglia_ppp*), the kernel
density estimates of the edge angle, length and both of the given ENS
network (*orgKDE_angle*, *orgKDE_length*, *orgKDE_both*), meshedness,
network_density, compactness of the given ENS network and the structural
properties of the initial triangulation (*network_extra*, *g2_degree*,
*triKDE_angle*, *triKDE_length*, *triKDE_both*).

The rejection sampling process picks an edge from the triangulation
based on the weights assigned to them related to the degree of the end
vertices. It keeps the edge if it violates connectivity constraints or
the kernel density estimates are in good terms, and rejects it
otherwise. The process also considers some additional checkings based on
meshedness, network_density, compactness and dissection artifact. A
helping function **predictKDE** computes the probability for the chosen
edge given the kernel density estimation of the real ENS network.

The function **predictKDE** computes the probability for the chosen edge
during the rejection sampling process given the kernel density
estimation of the real ENS network. The parameter *kde* is the given
kernel density estimate of the ENS network, *data* is the data point
(edge) for which prediction to be done and *mode* is an integer; either
1/2/3; 1 for angle, 2 for length; 3 for both.

### ENS network

The functions **deterministicEdges** and **rejectionSampling** described
above are called from this function **generateNetworkEdges** to execute
the ENS network generation process. **generateNetworkEdges** creates the
final network object (*g2_lin*), plots the corresponding kernel density
estimates of edge angle and length to visualize the goodness of the edge
sampling algorithm, and computes the earth mover’s distance (EMD)
between the given ENS and the simulated networks.

The coordinates of the end vertices of the edges of the generated
network are saved as *.csv* file in the corresponding sample’s output
subdirectory in **spatial-neuro/modelGanglionicNetwork/Outputs/**.

## Putting everything together

**Code: Main.R**

We have described all the R functions implemented for the generative
modeling. We write a **main** function to call these functions
sequentially with proper parameters, maintain correct folder and file
paths and construct the intensity profile of the neurons in Fiji
**before** moving to the neuron generation part.

## Generating neuron centers

**Code:GenerateNeuronCenters.R**

At this point, we have a simulated realization of the ENS ganglia
centers and a simulated ENS network. We also have an intensity profile
(a greyscale image) for the location of the ganglionic neurons, centered
around the ganglia and influenced by the network. Now, we genearte the
neurons and complete the generative modeling process. The function
**generateNeuronCenters** chooses the intensity profile image with file
chooser.

Initially, the location of the neurons are generated with a homogeneous
hardcore-Strauss process of given intensity (*b_val*), interaction
parameter (*g_val*), interaction distance (*r_val*) and hardcore
distance (*h_val*). Then inhomogeneity is incorporated by thinning the
generated point pattern with the intensity profile.

The coordinates of the generated neuron centers are saved as *.csv* file
in the corresponding sample’s output subdirectory in
**spatial-neuro/modelGanglionicNetwork/Outputs/**.

# Visualization

The saved generated ganglia centers, the ENS network and the generated
neuron centers can plotted together (using Fiji or R functionalities) to
construct a complete visualization of the simulated ganglionic network.

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-baddeley2015spatial" class="csl-entry">

Baddeley, Adrian, Ege Rubak, and Rolf Turner. 2015. *Spatial Point
Patterns: Methodology and Applications with R*. London: Chapman;
Hall/CRC. <http://www.crcpress.com/books/details/9781482210200/>.

</div>

<div id="ref-schindelin2012fiji" class="csl-entry">

Schindelin, Johannes, Ignacio Arganda-Carreras, Erwin Frise, Verena
Kaynig, Mark Longair, Tobias Pietzsch, Stephan Preibisch, et al. 2012.
“Fiji: An Open-Source Platform for Biological-Image Analysis.” *Nature
Methods* 9 (7): 676–82.

</div>

</div>
