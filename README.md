# GCDmodel

#### GCDkit routines for TE modelling

The recommended way to install these is by using

`devtools::install_github("jfmoyen/GCDmodel")`

This will ensure that you get the most recent code. Alternately, go to the "release" section of this page and download the version appropriate for your system.

# Version 0.7.0: 

-   Not many visible (user-facing) changes, but under the hood, a lot of code restructuring.
-   The GCDmodel version of `millications()` stops trying to emulate GCDkit and has been rewritten from scratch
-   The `Zrsat_XXX` functions, as well as `GCDmodel::millications()` should now be able to digest WR or milli (respectively wt. pct oxides or millications) in all(?) sensible formats: data frame, matrix or (named) vector.
-   Two global variables are now visible, `periodicTable` and `petroOxides`. `periodicTable`, from <https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee>, includes the atom mass and quite a few other properties
for all atoms in the peridoic table. `petroOxides` includes oxide information (mass, nb of cations and anions...) 
for commonly used oxides in petrology. They are exported and are visible with `GCDmodel::periodicTable` and `GCDmodel::petroOxides`

# Version 0.6.0: 

-   The saturation code for zircon has been modified. It does not use GCDkit code anymore, but internal versions - in file ZrSaturation.R

-   GCDkit's millications() is still used, if GCDkit is loaded. Otherwise we use a local clone.

-   This allows to select more easily the model to use. Now, the functions that use zircon saturation (i.e. `correctZrnSat` and `correctZrnMnzSat`) accept two arguments, `SatModel` and `SatArgs`. `SatModel` is the name of a saturation model (Boehnke, WH or Crisp). `SatArgs` are additional arguments passed to the model: TK and Zr are always passed, but other models require more parameters. Crisp, for instance, needs xH2O and P. In this case they should be passed with `SatArgs=list(H2O= ... ,Pbar= ...)`. This should allow to define relatively simply new zircon saturation functions, using any equation of choice. The header should look like the following: `ZrSat_XXX<-function(WR,milli=NULL,TK=NULL,Zr=NULL,...)` where XXX is the name by which the model will be called. Before the ..., include any other potentially useful argument, for instance `ZrSat_Crisp<-function(WR,milli=NULL,TK=NULL,Pbar=NULL,Zr=NULL,H2O=NULL,...)`
