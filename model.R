#!/usr/bin/Rscript

# The dynamic tidal system modeling.
# Read detailed information below:
# 1) Intensity (I) is associated with Square projection (S) on the plate plane
# 2) Orbit is simple circle (e = 1)
# 3) An inclination angle of the orbit is placed on eye-beam
# 4) Reference point is placed in the geometrical center of bigger planet

# Astronomical parameters
pi <- 3.1415926

# Model parameters 
# radius of bigger planet
rb <- 1

# radius of smaller planet
rs <- 1

# real distance between planets (radius of circle orbit)
dist <- 10

# orbit inclination (in radians)
incl <- pi / 4

# position angle
pa <- 0

# visible distance
vdist <- function(dist,incl, pa)
{
  vdist <- coord * cos(incl) / sin(pa)
}

# flux
flux <- function()
{
  if (a > rb + rs)
    flux <- pi*rb**2 + pi*rs**2 
}
