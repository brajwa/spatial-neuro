Quantify spatial inhomogeneity and anisotropy in the peripheral nerve
cross-sections
================
Abida Sanjana Shemonti
July 27, 2022

# What are we doing?

We are demonstrating a pipeline that represents the segmented
unmyelinated axons in several vagus and pelvic nerve cross-sections as
spatial point patterns, and computes Sinkhorn distance between certain
spatial features of the point patterns and their corresponding images,
to quantify spatial inhomogeneity and anisotropy in the peripheral nerve
cross-sections.

# Data

The directory **spatial-neuro/quantifyVagusNerve/Data/** contains the
following subdirectories:

-   **U-Net predictions:** This subdirectory contains images (*.png*) of
    the segmented unmyelinated axons in the vagus and pelvic nerve
    cross-sections under consideration. Plebani et al. (2022) describes
    the segmentation process.

-   **Inputs:** There are *.csv* and *.rds* files in this subdirectory.
    The *.csv* files contain several morphometric characteristics of the
    segmented unmyelinated axons in each nerve cross-section, along with
    their centroid locations, which we use in particular to construct
    the spatial point patterns. The morphometric characteristics are
    extracted using an open source image processing package Fiji
    (Schindelin et al. (2012)). The *.rds* files, drawn using R graphics
    tools, store the outer boundary and inner holes of the
    cross-sections that are used for spatial point pattern construction
    also. These are the main input files for the analysis.

-   **Demo Inputs:** This subdirectory contains a small subset of the
    inputs for quick demonstration purpose.

-   **Supporting Files:** Some additional information of the nerve
    samples, such as sample ID, nerve location and sex, are saved in
    this subdirectory. The pair-wise orientation of the cross-sections
    are computed at a certain stage of the pipeline. Those orientations
    are also stored inside this subdirectory.

# Codes

The directory **spatial-neuro/quantifyVagusNerve/Codes/** contains the R
code files. We utilize the R packages *spatstat* for spatial point
pattern analysis, *T4transport* and *Barycenter* for optimal
transportation problems extensively, along with the others. The R files
can be run using RStudio terminal or Linux terminal. Here we demonstrate
with the RStudio terminal commands. The changes required (if any) to run
the codes from Linux terminal are mentioned inside the code files.

## Computation of the scaling factor

**Code: ScalingFactorComputation.R**

The nerve cross-sections are scaled with the largest range along the
y-axis.

``` r
> Rscript --no-save ScalingFactorComputation.R
```

## Configuration of the interaction distance

**Code: LfunctAnova2.R**

The nerve cross-sections under consideration can vary a lot in size,
resulting into quite different ranges of interaction in the constructed
point patterns even after scaling. We find out a range of interaction
distance which is common in all the point patterns by investigating
their inhomogeneous L-function characteristics. We also figure out a
particular interaction distance, to be used for computations related to
*local* inhomogeneity and anisotropy, using analysis of variance that
can moderately separate the inhomogeneous L-function characteristics
between the vagus and pelvic samples.

This R file requires no command-line argument, the last line of the
outputs shows the configured interaction distance.

``` r
> Rscript --no-save LfunctAnova2.R
```

**Note:** This pre-processing step is completely dependent on the set of
nerve cross-sections under consideration. In case of addition or removal
of any cross-section, the interaction distance needs to be recomputed.

## Modification for local inhomogeneity and anisotropy

**Code: Kinhomsector.R, LocalKSector.R**

The R *spatstat* package provides several separate functions to compute
the inhomogeneous, anisotropic and local versions of the K- and
L-functions. We combine several of them to meet our requirement of
analysis. These R files are sourced during subsequent computations.

## Computation of Sinkhorn distance

The R code files for the computation of the Sinkhorn distance between
nerve cross-sections require the command line arguments to be passed in
the following order:

``` r
> Rscript --no-save code_file_name analysis_type scaling lambda version intr_dist
```

-   **analysis_type:** We analyze the spatial point patterns in terms of
    basic spatial intensity (*“basic_density”*) , local inhomogeneous
    L-function (*“local_linhom”*), local inhomogeneous L-function with
    horizontal sector (*“local_linhom_sector_horizontal”*) and local
    inhomogeneous L-function with vertical sector
    (*“local_linhom_sector_vertical”*). The horizontal and vertical
    sectors covers 15 degree around the horizontal and vertical axes
    respectively.

-   **scaling:** *0* - if the spatial point patterns are not to be
    scaled, *1* - if the spatial point patterns are to be scaled. We use
    the scaled point patterns.

-   **lambda:** The tunable entropic regularizer parameter for Sinkhorn
    distance computation.

-   **version:** *“complete”* - to use the complete input set,
    *“demo”* - to use the small subset of input for quick demonstration.

-   **intr_dist:** The interaction distance computed in LfunctAnova2.R

### Computation based on spatial point patterns

**Code: T4transportSpatialIntensity.R, T4transportSpatialFeature.R**

The R package *T4transport* is used for computing Sinkhorn distance
based on the spatial point patterns directly. At first, the Sinkhorn
distance between every pair of point patterns for basic spatial
intensity is computed, assuming each point in the point pattern carry an
uniform weight. While computing distance between a pair of point
patterns, one pattern is kept fixed and the other one is rotated by
every 45 degree. The orientation that gives the smallest Sinkhorn
distance is recorded as the orientation to be used for computations of
other spatial features.

``` r
> Rscript --no-save T4transportSpatialIntensity.R "basic_density" 1 0.01 "demo"  0.0046
```

The Sinkhorn distance between every pair of point pattern for every
orientation, their minimum, maximum, mean and variance and the Sinkhorn
distance matrix - everything is saved as *.xlsx* file, along with the
plots of the point patterns and their corresponding density maps, inside
**spatial-neuro/quantifyVagusNerve/Plots/** directory. The optimal
orientations are saved in
**spatial-neuro/quantifyVagusNerve/Data/Supporting Files/Angle of
Rotation/** directory as *.xlsx* file.

Then, the Sinkhorn distance for other spatial features can be computed
using the pre-computed orientations, where each point in the point
pattern carries weight equal to the value of the spatial feature. The
results are recorded in **spatial-neuro/quantifyVagusNerve/Plots/**
directory in a similar way.

``` r
> Rscript --no-save T4transportSpatialFeature.R "local_linhom" 1 0.01 "demo"  0.0046
> Rscript --no-save T4transportSpatialFeature.R "local_linhom_sector_horizontal" 1 0.01 "demo"  0.0046
> Rscript --no-save T4transportSpatialFeature.R "local_linhom_sector_vertical" 1 0.01 "demo"  0.0046
```

### Computation based on images of the spatial features

**Code: BarycenterSpatialIntensity.R, BarycenterSpatialFeature.R**

The R package *Barycenter* is used for computing Sinkhorn distance based
on the images of the spatial features of the point patterns, rather than
the point patterns themselves. The codes are run and the results are
stored as before.

``` r
> Rscript --no-save BarycenterSpatialIntensity.R "basic_density" 1 0.01 "demo"  0.0046

> Rscript --no-save BarycenterSpatialFeature.R "local_linhom" 1 0.01 "demo"  0.0046
> Rscript --no-save BarycenterSpatialFeature.R "local_linhom_sector_horizontal" 1 0.01 "demo"  0.0046
> Rscript --no-save BarycenterSpatialFeature.R "local_linhom_sector_vertical" 1 0.01 "demo"  0.0046
```

There is another code file in the directory named
*BarycenterSpatialIntensity_2.R* that helps in running the codes with
input files in batches.

# Output and visualization

**Code: VisualizeMDS.R**

The directory **spatial-neuro/quantifyVagusNerve/Plots/** will contain
the output plots and files in separate subdirectories. The *.xlsx* files
saved here contains the Sinkhorn distance matrices and we can visualize
the nerve samples in the Sinkhorn space of the corresponding spatial
features by applying multi-dimensional scaling on them.

``` r
> Rscript --no-save VisualizeMDS.R
```

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-cuturi2013sinkhorn" class="csl-entry">

Cuturi, Marco. 2013. “Sinkhorn Distances: Lightspeed Computation of
Optimal Transport.” *Advances in Neural Information Processing Systems*
26.

</div>

<div id="ref-plebani2022high" class="csl-entry">

Plebani, Emanuele, Natalia P Biscola, Leif A Havton, Bartek Rajwa, Abida
Sanjana Shemonti, Deborah Jaffey, Terry Powley, Janet R Keast, Kun-Han
Lu, and M Murat Dundar. 2022. “High-Throughput Segmentation of
Unmyelinated Axons by Deep Learning.” *Scientific Reports* 12 (1): 1–16.

</div>

<div id="ref-schindelin2012fiji" class="csl-entry">

Schindelin, Johannes, Ignacio Arganda-Carreras, Erwin Frise, Verena
Kaynig, Mark Longair, Tobias Pietzsch, Stephan Preibisch, et al. 2012.
“Fiji: An Open-Source Platform for Biological-Image Analysis.” *Nature
Methods* 9 (7): 676–82.

</div>

</div>
