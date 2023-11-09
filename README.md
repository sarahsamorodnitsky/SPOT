# SPOT 

## Description

This R package contains the code to perform the SPatial Omnibus Test (SPOT). SPOT takes a range of radii, starting with 0, calculates a spatial summary of the point locations 
across images, and calculates an omnibus p-value describing the strength of association between the spatial summary and the outcome across radii. 

## Installation and Dependencies

To install SPOT, run the following:
```{r}
devtools::install_github("sarahsamorodnitsky/SPOT")
```
This package requires the following dependencies: spatstat, survival, ACAT, dplyr, tidyselect, svMisc. 
