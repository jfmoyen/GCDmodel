# GCDmodel
GCDkit routines for TE modelling

The recommended way to install these is by using 

devtools::install_github("jfmoyen/GCDmodel")

This will ensure that you get the most recent code. Alternately, go to the "release" section of this page and download the version appropriate for your system.

NB: Starting from version 0.6, the saturation code used for zircon stopped being the GCDkit code. GCDmodel is using its own code. GCDkit's millications() is still used, if GCDkit is loaded. Otherwise we use a local clone.
