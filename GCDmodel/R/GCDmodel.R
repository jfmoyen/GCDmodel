#### TODO
# Clean a bit Rcrust dependencies
# Missing GCDkit functions (see where they end up ??)
# Matt's upper cases (well, duh).

# Avoid redefining zrSat and mzSat - move to a library outside of main GCDkit ??? Allow usage without GCDkit + data loaded
# Better : source the plugin

# Require ?
# GCDkit, stats

#### Utilities for modelling
.onLoad<-function(lib, pkg){
    on.exit(options("show.error.messages"=TRUE))
    options("warn"=-1)

    packageStartupMessage("GCDmodel module - 0.5.01", appendLF = TRUE)

    # If we do not have access to GCDkit functions, redefine them here
    if(!exists("millications")){

      packageStartupMessage("GCDkit not available, redefining key functions - expect trouble", appendLF = TRUE)

      ##Definitions, variables
      lree.nm<<-c("Th","La","Ce","Pr","Nd","Sm","Gd")

      mw<<-c(232.04,138.91,140.12,140.91,144.24,150.3,157.2,30.974,15.999)
      names(mw)<<-c(lree.nm,"P","O")


      MW<<-c(60.07,79.865,101.961,159.677,71.839,71.839,70.937,40.299,56.077,61.979,94.195,18.015,44.008,141.943,18.998,32.000)
      names(MW)<<-c("SiO2","TiO2","Al2O3","Fe2O3","FeO","FeOt","MnO","MgO","CaO","Na2O","K2O","H2O.PLUS","CO2","P2O5","F","S")
      # SiO2     TiO2    Al2O3    Fe2O3      FeO     FeOt      MnO      MgO      CaO     Na2O
      # 60.078   79.865  101.961  159.677   71.839   71.839   70.937   40.299   56.077   61.979
      # K2O H2O.PLUS      CO2     P2O5        F        S
      # 94.195   18.015   44.008  141.943   18.998   32.000

      x.atoms<<-c(1,1,2,2,1,1,1,1,1,2,2,2,1,2,1,1)
      names(x.atoms)<<-names(MW)

      millications<<-function (x, print = FALSE, save = FALSE)
        # This is the GCDkit millications()
      {
        x <- as.matrix(x)
        if (ncol(x) == 1) {
          x <- t(x)
        }
        oxides <- c("SiO2", "TiO2", "Al2O3", "Fe2O3", "FeO", "FeOt",
                    "MnO", "MgO", "CaO", "Na2O", "K2O", "H2O.PLUS", "CO2",
                    "P2O5", "F", "S")
        oxides <- oxides[oxides %in% colnames(x)]
        ee <- sapply(1:nrow(x), function(f) {
          z <- x[f, oxides]/MW[oxides] * x.atoms[oxides] * 1000
          return(z)
        })
        milli <- t(ee)
        milli[is.na(milli)] <- 0
        rownames(milli) <- rownames(x)
        .atoms.from.formula <- function(ox) {
          z <- gsub("[0-9]", "", ox)
          z <- sapply((strsplit(z, "O")), paste, collapse = "")
          return(z)
        }
        if (save) {
          oxides <- oxides[c(1:11, 14)]
        }
        results <- milli[, oxides, drop = FALSE]
        ee <- paste(.atoms.from.formula(oxides), ".m", sep = "")
        i <- grep("^Fe.m$", ee)
        if (length(i) == 2)
          ee[i] <- c("Fe3.m", "Fe2.m")
        ee <- gsub("^Fet.m$", "Fe.m", ee)
        colnames(results) <- ee
        if (print)
          print(results)
        if (save) {
          assign("milli", milli, .GlobalEnv)
          assign("results", results, .GlobalEnv)
          invisible()
        }
        return(milli)
      }

      zrSaturation<<-function(cats=milli,T=0,Zr=filterOut(WR,"Zr",1)){
        on.exit(options("show.error.messages"=TRUE))
        # if(T==0){
        #   T<-winDialogString("Temperature (degrees C)","750")
        #   if(is.null(T)){cat("Cancelled.\n");options("show.error.messages"=FALSE);stop()}
        # }

        T<-as.numeric(T)+273.15
        oxides<-c("SiO2","TiO2","Al2O3","FeOt","MnO","MgO","CaO","Na2O","K2O","P2O5")
        oxides<-oxides[oxides%in%colnames(cats)]
        cats<-cats[,oxides,drop=FALSE]
        if(length(Zr)>1) cats[names(Zr),oxides]

        # NEW
        sums<-apply(cats,1,sum,na.rm=TRUE)
        ee<-sapply(1:nrow(cats),function(i){
          z<-cats[i,]/sums[i]
          return(z)
        })
        cats<-t(ee)

        #cats<-normalize2total(cats)/100

        M<-(cats[,"Na2O"]+cats[,"K2O"]+2*cats[,"CaO"])/(cats[,"Al2O3"]*cats[,"SiO2"])

        # Watson and Harrison (1983)
        DZr<-exp(-3.8-0.85*(M-1)+12900/T)
        Zr.sat<-497644/DZr

        DZr<-497644/Zr
        DZr<-as.vector(DZr)
        TZr.sat.C<-12900/(log(DZr)+3.8+0.85*(M-1))-273.15

        # Boehnke et al. (2013
        DZrB<-exp(10108/T-1.16*(M-1)-1.48)
        Zr.satB<-497644/DZrB

        DZrB<-497644/Zr
        DZrB<-as.vector(DZrB)
        TZr.satB.C<-10108/(log(DZrB)+1.16*(M-1)+1.48)-273.15


        y<-cbind(M,Zr,round(Zr.sat,1),round(TZr.sat.C,1),round(Zr.satB,1),round(TZr.satB.C,1))
        colnames(y)<-c("M","Zr.obs","Zr.sat","TZr.sat.C","Zr.sat (Boehnke)","TZr.sat.C (Boehnke)")
        # if(nrow(y)>1) y<-formatResults(y) else rownames(y)<-rownames(cats)
        # if(!getOption("gcd.shut.up"))print(y)
        # assign("results",y,.GlobalEnv)
        invisible(y[,-2])
      }

    }else{packageStartupMessage("Nothing to redefine", appendLF = TRUE)}


    packageStartupMessage("...Done!", appendLF = TRUE)
    #winMenuAdd("Plugins")
    invisible()
}

############### INTERFACE ############
# Melting functions:
#function(kd,c0,pm,
#          min.props,
#          melt.arg)  # melt.arg$mjrs, melt.arg$trc, melt.arg$PP, melt.arg$TT, melt.arg$H2O
# invisible(
# list(cL=.sanitize(cL), # A vector
#      cmins=cmins, # A matrix with proper names of phases/elements
#      kd=kd,
#      DD=.sanitize(ee$DD), # A vector
#      min.props=m.pr, # A vector
#      dont=dont)

# Correction functions:
# function(kd,c0,pm,
#          cmins,min.props,
#          melt.arg,  # melt.arg$mjrs, melt.arg$trc, melt.arg$PP, melt.arg$TT, melt.arg$H2O
#          dont)
# invisible(
#   list(cL=.sanitize(cL), # A vector
#        cmins=cmins,
#        kd=kd,
#        DD=.sanitize(ee$DD),# A vector
#        min.props=m.pr, # A vector
#        dont=dont)
################################################

######## Kd table cleaning
ppxCleanKdTable<-function(kd,ppxPhases,interactive=T){
  #' Prepare a Kd table for use with perplex models:
  #' Create multiple rows for duplicate phases (min_1 min_2 etc)
  #' Add rows with zeros for missing phases
  #' Replace NA's by zeros
  #' NB keep all Kds including phases not in perplex file
  #' @param kd: a kd table (minerals in rows, elements in columns)
  #' @param ppxPhases: "proper" phases in a perpleX model cube
  #' @param interactive: Give info messages in pop-up windows
  #' @return a matrix of Kds, with names matching ppxPhases

  #' @export

  used.mins<-union(ppxPhases,rownames(kd))

  kd<-as.matrix(kd)
  kd.new<-matrix(NA,nrow=length(used.mins),ncol=ncol(kd))
  rownames(kd.new)<-used.mins
  colnames(kd.new)<-colnames(kd)

  .replace.Kd<-function(old.name,new.name){
    on.exit(options("show.error.messages"=TRUE))
    if(!old.name%in%rownames(kd)){
      message<-paste("Phase ",old.name," not found in the kd file!\nReplacing by zeros.....",sep="\n")
      if(interactive){
        winDialog(type="ok",message)
      }
      cat(message,"\n")
      options(show.error.messages=FALSE)
      kd.new[new.name,]<<-rep(0,ncol(kd.new)) # NEW
      return()
    }
    # Should never happen ?
    # if(!new.name%in%allPhases){
    #   message<-paste("Phase ",new.name," not found in the Perplex file!",sep="\n")
    #   winDialog(type="ok",message)
    #   options(show.error.messages=FALSE)
    #   return()
    # }

    kd.new[new.name,]<<-kd[old.name,]
  }

  ee<-lapply(used.mins,function(i){
    old<-gsub("_[0-9]{1}$","",i)
    new<-i
    .replace.Kd(old,new)
  })
## if there is Ph_1 and Ph-2 in perplex, and Ph in Kd, this routine keeps the 3
  # So we end up with Ph, Ph_1, Ph_2
  # Probably not a problem ?
  kd.new[is.na(kd.new)]<-0
  return(kd.new)

}

##################################

massBal<-function(model,elt,tol=1e-6){
  #' Check whether a model honours mass balance constrains
  #' (for debugging purposes...)
  #' @param model. A model, generated by any of the melting/correction functions, i.e
  #' a list with the appropriate components
  #' @param elt Character. The (trace) element to check mass balance for
  #' @param tol real. Tolerance value, how close must the two values be (round of small numerical errors)
  #'
  #' @export
  ss<-model$c0[elt]
  cat("C(",elt,") in source: ",ss,"\n",sep="")
  ee<-sum(model$cmins[,elt]*model$min.props)*(1-model$FF/100)+model$cL[elt]*model$FF/100
  cat("C(",elt,") partitionned: ",ee,"\n",sep="")
  if(abs(ee-ss)<tol){
    cat("Ok\n")
  }else{
    cat("Error! Difference of ",ee-ss," ppm, i.e. ",(ee-ss)/ss*100," %\n",sep="")

  }
}


DummyPM<-function(kd,c0=NULL,pm=NULL,cmins=matrix(),min.props=NULL,melt.arg=list(),dont=character(0)){

  #' Test function

  #' @export

  cL<-c(1,2)
  names(cL)<-c("dummy1","dummy2")

  cmins<-matrix(NA,
             nrow=length(.sanitize(min.props)),
             ncol=ncol(kd)
             )
  rownames(cmins)<-names(.sanitize(min.props))
  colnames(cmins)<-colnames(kd)

  DD<-numeric(0)
  cS<-numeric(0)

  return(list(cL=cL,cmins=cmins,min.props=min.props,cS=cS,kd=kd,DD=DD))
}

BatchPM<-function(kd,c0,pm,cmins=matrix(),min.props,melt.arg=list(),dont=character(0)){

  #' Trace elements concentrations for a batch process
  #'
  #' More commonly used for batch melting, but does actually also work with crystalization.
  #' Make sure to use the same mineral/element names in all the variables
  #' @param kd Named matrix. Partition coefficients
  #' @param c0 Named vector. The composition of the source
  #' @param pm Scalar. Degree of partial melting (0-100)
  #' @param cmins Matrix. Composition of minerals (not used, for consistency)
  #' @param min.props Named vector. Mineral proportions (normalized to 1)
  #' @param melt.arg List. Further useful arguments (mostly used for correction modules, kept here for consistency)
  #' @param dont Character vectors. Do not calculate for these elements. Not used here, kept for consistency (mostly used in correction modules)
  #' @return list containing the following slots:
  #' \itemize{
  #'   \item \code{c0} Named vector. Composition of the source, as passed to the function
  #'   \item \code{cL} Named vector. Composition of the liquid
  #'   \item \code{cmins} Named matrix. Composition of all individual solid mineral phases
  #'   \item \code{cS} Named vector. Composition of the whole solid
  #'   \item \code{min.props} Named vector. Mineral proportions in the restite (and expressed as (0-1) of the restite)
  #'   \item \code{FF} Real. Melt proportion, 0-100, as passed to the function
  #'   \item \code{kd} Named matrix. The Kd table, here similar to the input (included for consistency with other modelling routines)
  #'   \item \code{DD} Named vector. The D (bulk distribution) values.
  #' }

  #' @export

  c0<-.sanitize(c0)
  min.props<-.sanitize(min.props,normalize.TO=1)

  which.elems<-intersect(names(c0),colnames(kd))
  which.mins<-names(min.props)

  # Test whether Kds for all minerals are available
  missing<-setdiff(names(min.props),rownames(kd))
  if(length(missing)!=0){
    stop(paste("Missing Kd value for mineral(s)",missing,"\n"))
  }

  # Keep only the required elements/minerals
  # drop=F is important to preserve the objects as matrices even if only one line
  c0<-c0[which.elems] # c0 has been sanitized to a vector by then
  kd<-kd[which.mins,which.elems,drop=F]

  # Do the real job...
  DD<-min.props%*%kd
  FF<-pm/100
  cL<-c0/(DD+FF*(1-DD))
  ee<-sapply(names(min.props),function(i){
    z<-cL[1,which.elems]*kd[i,which.elems]
    return(z)
  },simplify=TRUE)
  cmins<-t(ee)
  cmins<-cmins[which.mins,which.elems,drop=F]

  cS<-cL[1,which.elems]*DD[1,which.elems]

  invisible(list(c0=c0,cL=.sanitize(cL),cmins=cmins,cS=cS,min.props=min.props,FF=pm,kd=kd,DD=.sanitize(DD)))
}

correctZrnSat<-function(kd,c0,pm,cmins,min.props,melt.arg=list(),dont=character(0),SatModel="Boehnke"){

  #' Correction for zircon saturation, Watson & Harrison 1983
  #' @param kd Named matrix. Partition coefficients. Must include Zircon (Zrn)
  #' @param pm Scalar. Degree of partial melting (0-100)
  #' @param cmins Matrix. Composition of minerals (not used, for consistency)
  #' @param min.props Named vector. Mineral proportions (normalized to 1)
  #' @param melt.arg List. Further useful arguments.
  #' Includes here the melt major elements (melt.arg$mjrs),trace (melt.arg$trc), temperature in KELVIN (melt.arg$TT)
  #' @param dont Character vector. Do not calculate for these elements. Use to avoid overptinting previous saturation calc (i.e. if you use a saturation routine /after/ this one, )
  #' @param SatModel model to use Currently "Boehnke" (the default; Boehnke et al. 2013) or "WH" (Watson & Harrison 1983)
  #' @return list containing the following slots:
  #' \itemize{
  #'   \item \code{c0} Named vector. Composition of the source, as passed to the function
  #'   \item \code{cL} Named vector. Composition of the liquid
  #'   \item \code{cmins} Named matrix. Composition of all individual solid mineral phases
  #'   \item \code{cS} Named vector. Composition of the whole solid
  #'   \item \code{min.props} Named vector. Mineral proportions in the restite (and expressed as (0-1) of the restite)
  #'   \item \code{FF} Real. Melt proportion, 0-100, as passed to the function
  #'   \item \code{kd} Named matrix. The Kd table, here similar to the input (included for consistency with other modelling routines)
  #'   \item \code{DD} Named vector. The D (bulk distribution) values.
  #'   \item \code{dont} Character vector. Elements that should not be touched by further correction routines
  #'   \item \code{tot.mass} Real. The real mass of the system (it increased since we added Zr, so it is now 100% + [Zr])
  #'   \item \code{corr.masses} Named vector. Corrected masses for each mineral, sum = tot.mass (NOT 100)
  #'   \item \code{si.used} Real. Amount of SiO2 used to form ZrSiO4, and thus technically not available to major phases (wt %)
  #'   \item \code{o.used} Real. Amount of O2 used to form ZrSiO4, and thus technically not available to major phases (wt %)
  #'   \item \code{sat} Named vectors. Saturation parameters
  #' }
  #' @export

  c0<-.sanitize(c0)
  m.pr<-.sanitize(min.props,normalize.TO=1)

  # Protection against Matt converting oxides to upper case !!
  names(melt.arg$mjrs)<-.TrueNames(names(melt.arg$mjrs))

  milcats <- millications(melt.arg$mjrs)
  # Prevent against various forms of Fe
  .cleanMillicats(milcats)

  # Model selector
  if(SatModel=="Boehnke"){
    zrSatName <- "Zr.sat (Boehnke)"
    ZrTName <- "TZr.sat.C (Boehnke)"
  }else{
    zrSatName <- "Zr.sat"
    ZrTName <- "TZr.sat.C"
  }
  # call Saturation plugin, Zr just an arbitrary value as it is not needed
  sat.data <- zrSaturation(milcats,T=melt.arg$TT-273.15,Zr=melt.arg$trc["Zr"])
  Zr.sat <- sat.data[zrSatName]
  Zr.T <- sat.data[ZrTName]
  M <- sat.data["M"]

  # Housekeeping
  all.but.zrn <- rownames(cmins)

  # Calculate masses of everything, for 100 g of system
  mass<-c((1-pm/100)*m.pr*100,pm)
  names(mass)[length(mass)]<-"Melt"

  # .. however, the system is now bigger by [Zr]
  tot.mass <- 100 + c0["Zr"]/1e4 # 100 g + the Zr we added

  if(melt.arg$trc["Zr"]<=Zr.sat){
    ## The melt is not saturated in zircon, do nothing more!
    # min.props are inherited, don't change
    cm<-cmins # cmins are inherited, don't change
    cL<-melt.arg$trc # inherited also
    zrc.prop<-0

    # We update the mass of each phase by the amount of Zr it takes
    zr.mass<-c(cm[,"Zr"],cL["Zr"])*mass/1e6
    mass<-mass+zr.mass

    # Various masses & outputs
    si.used<-0
    o.used<-0


  }else{
    ## The melt is saturated
    # We know that cL = Zr.sat
    # First we recalc the Zr amount in the other minerals based on cL
    cZr <- kd[all.but.zrn,"Zr"] * Zr.sat

    # So we now now the mass of Zr that goes to each phase
    zr.mass<-c(cZr,Zr.sat)*mass/1e6

    # Zr mass balance
    zr.excess <- c0["Zr"]/1e4 - sum(zr.mass)

    # Excess Zr is affected to Zrn
    zrc.mass <- zr.excess/0.497644 # relative to the whole system (100 g)

    # Forming Zrn requires other elements !
    si.used<- zrc.mass*0.327765 # in g per 100 g of system
    o.used <- zrc.mass*0.174570

    # We take a bit of that SiO2 and O2 from everywhere - no partitionning here
    si.mass<-mass/sum(mass)*si.used
    o.mass<-mass/sum(mass)*o.used

    # Final mass of the system
    mass<-c(mass+zr.mass-si.mass-o.mass,zrc.mass)
    names(mass)[length(mass)]<-"Zrn"

    # Expressed as proportions
    foo<-mass[c(rownames(cmins),"Zrn")]
    m.pr<-foo/sum(foo)

    # # Zr mass balance
    # used.Zr <- Zr.sat * pm/100 + sum(cZr*m.pr) * (1-pm/100)
    #
    # # Excess Zr is affected to Zrn
    # zr.excess <- c0["Zr"] - used.Zr
    # zrc.prop <- zr.excess/497644 # relative to the whole system (0-1)
    # m.pr["Zrn"] <- zrc.prop / (1 - pm/100) # relative to the solid

    # m.pr <- m.pr/sum(m.pr) ##### ???? changes the major elem mass prop and alters Zr mass balance...
    # But not normalizing will alter th other traces mass balance....

    #... and we update the Kd of Zr in zircon
    kd["Zrn","Zr"]<-497644/Zr.sat

    # Now we have three sort of elements:
    # 1. "ordinary" elements
    # 2. "dont" elements, not recalculated, for which we inherit the original values
    # 3. Zr, that gets a special treatment

    # 1. ordinary elements - we do partitioning with the new min.props
    ee<-BatchPM(kd=kd,c0=c0,pm=pm,
                min.props=m.pr)
    cL<-ee$cL
    cm<-ee$cmins

    # 2. "dont" elements - we inherit the values
    cL[dont] <- melt.arg$trc[dont]
    #all.but.zrn<-setdiff(rownames(cm),"Zrn")
    cm[all.but.zrn,dont]<-cmins[,dont]

    # 3. Zr - special case
    cL["Zr"] <- Zr.sat
    cm["Zrn","Zr"]<-497644

    # We must adjust the Zr conc to honour mass balance
    alpha<-(c0["Zr"]-497644*m.pr["Zrn"]*(1-pm/100)-Zr.sat*pm/100)/sum(cZr*m.pr[all.but.zrn]*(1-pm/100))
    cm[all.but.zrn,"Zr"]<-cZr*alpha
  }

  # D is an "observed" D that we calculate from cS/cL, in all cases
  cS<-.sanitize(m.pr%*%cm)
  DD<-cS/cL

  # Prepare the various outputs
  dont<-c(dont,"Zr")

  # Proportions relative to the whole system
  #all.but.zrn.masses<-c(m.pr[all.but.zrn]*100*(1-pm/100),pm)
  #av.mass<-tot.mass - zrc.prop*100

  # Normalize to the mass available
  #corr.masses<-c(all.but.zrn.masses*av.mass/100,zrc.prop*100)

  #nm<-names(corr.masses)
 # nm[nm=="FF"]<-"Melt"
  #names(corr.masses)<-nm

  # SiO2 used to make zircon, relative to whole system
 # si.used<-zrc.prop*327765/1e4 # in wt% of whole system
  # ZrSiO4 = 597644 ppm of Zr, 327765 ppm of SiO2, 174570 ppm of O2

  invisible(list(c0=c0,cL=.sanitize(cL),
                 cmins=cm,
                 cS=cS,
                 min.props=m.pr,
                 FF=pm,
                 kd=kd,
                 DD=.sanitize(DD),
                   dont=dont,
                 tot.mass=tot.mass,
                 corr.masses=mass,
                 si.used=si.used,
                 o.used=o.used,
                 sat=c(M,Zr.sat,Zr.T)))

}



correctMnzSat<-function(kd,c0,pm,cmins,min.props,melt.arg=list(),dont=character(0)){

  #' Correction for monazite saturation, Montel 1992 w. elements partitionning
  #' @param kd Named matrix. Partition coefficients. Must include Monazite (Mnz)
  #' @param pm Scalar. Degree of partial melting (0-100)
  #' @param cmins Matrix. Composition of minerals (not used, for consistency)
  #' @param min.props Named vector. Mineral proportions (normalized to 1)
  #' @param melt.arg List. Further useful arguments.
  #' Includes here the melt major elements (melt.arg$mjrs),trace (melt.arg$trc), water (melt.arg$H2O), temperature in KELVIN (melt.arg$TT)
  #' @param dont Character vector. Do not calculate for these elements. Use to avoid overptinting previous saturation calc (i.e. if you use a saturation routine /after/ this one, )
  #' @return list containing the following slots:
  #' \itemize{
  #'   \item \code{c0} Named vector. Composition of the source, as passed to the function
  #'   \item \code{cL} Named vector. Composition of the liquid
  #'   \item \code{cmins} Named matrix. Composition of all individual solid mineral phases
  #'   \item \code{cS} Named vector. Composition of the whole solid
  #'   \item \code{min.props} Named vector. Mineral proportions in the restite (and expressed as (0-1) of the restite)
  #'   \item \code{FF} Real. Melt proportion, 0-100, as passed to the function
  #'   \item \code{kd} Named matrix. The Kd table, here similar to the input (included for consistency with other modelling routines)
  #'   \item \code{DD} Named vector. The D (bulk distribution) values.
  #'   \item \code{dont} Character vector. Elements that should not be touched by further correction routines
  #'   \item \code{tot.mass} Real. The real mass of the system (it increased since we added Zr, so it is now 100% + [Zr])
  #'   \item \code{corr.masses} Named vector. Corrected masses for each mineral, sum = tot.mass (NOT 100)
  #' }
  #' @export

  c0<-.sanitize(c0)
  m.pr<-.sanitize(min.props,normalize.TO=1)

  ##Definitions, variables
  lree.nm<-c("Th","La","Ce","Pr","Nd","Sm","Gd")

  # Protection against Matt converting oxides to upper case !!
  names(melt.arg$mjrs)<-.TrueNames(names(melt.arg$mjrs))

  milcats <- millications(melt.arg$mjrs)
  # Prevent against various forms of Fe
  .cleanMillicats(milcats)

  W.R<-c(-520,-260,-116,31,177,460,750)
  names(W.R)<-lree.nm
  MW.REE<-mw[lree.nm]
  REE<-melt.arg$trc[lree.nm]
  M<-rep(0,length(lree.nm))

 # Assume 1% of U in monazite
  Xmz<-0.99

  ## Actual amount of LREE in the melt, atomic (uncorrected)
  S<-REE[lree.nm]/MW.REE[lree.nm]
  ST<-sum(S)

  # Call the Saturation plugin and get the saturation value of LREE + Th
  ee <- mzSaturationWithTh(cats=milcats,REE=REE,H2O=melt.arg$H2O,Xmz=Xmz,Temp=melt.arg$TT-273.15) # Temperature is in C
  AiT<-ee[1,"Sum.A.REE.sat"]

  # Housekeeping
  all.but.mnz <- rownames(cmins)

  # Calculate masses of everything, for 100 g of system
  mass<-c((1-pm/100)*m.pr*100,pm)
  names(mass)[length(mass)]<-"Melt"

  # .. however, the system is now bigger by [LREE]
  tot.mass <- 100 + sum(c0[lree.nm])/1e4 # 100 g + the LREE we added

  if(ST<=AiT){
    # The melt is not saturated in monazite, do nothing more!
    cm<-cmins # cmins are inherited, don't change
    cL<-melt.arg$trc

    # We update the mass of each phase by the amount of LREE it takes
    for(e in lree.nm){
      lree.mass<-c(cm[,e],cL[e])*mass/1e6
      mass<-mass+lree.mass
    }

    # Various masses & outputs
    o.used<-0
    p.used<-0

  }else{
    # The melt is saturated, we calculate the mnz composition

    # Partitioning equation: gamma=(Xmelt/XMnz)
    gamma<-exp(W.R/(melt.arg$TT))

    # We do not know MiT and A, will optimize the equation for MiT
    .compute<-function(MiT){
      A<-gamma*AiT*S/(MiT+gamma*AiT)
      ee<-(sum(A)-AiT)^2
      return(ee)
    }

    m.guess<-0
    out<-optim(m.guess,.compute,method="BFGS",lower=0,upper=1)
    R2<-out$value
    MiT<-out$par

    # Concentrations in the liquid, in ppm
    A<-gamma*AiT*S/(MiT+gamma*AiT)
    lree.c<-A*MW.REE

    # Concentration in all the other phases, from Kd
    cphases<-sapply(lree.nm,function(e){
             return(kd[all.but.mnz,e] * lree.c[e])
               })

    # Extra mass
    min.extra.mass<-sapply(rownames(cphases),function(p){
      ee<-sum(cphases[p,]*mass[p])/1e6
      return(ee)
    })
    liq.extra.mass<-sum(mass["Melt"]*lree.c)/1e6
    # Corrected masses
    mass<-mass+c(min.extra.mass,liq.extra.mass)

    # REE : Atomic REE concentration in the system (vec)
    # c0[lree.nm] : Weight REE concentration in the system (vec)
    # lree.c : Weight LREE concentration in the liquid (vec)
    # cphases : Weight LREE concentration in each phase (mat)
    # mass : mass of each phase, melt and Mnz included (vec)

    cphases.melt<-rbind(cphases,lree.c)
    rownames(cphases.melt)<-names(mass)

    # Mass of LREE in each component of the system (mat)
    lree.mass<-sapply(lree.nm,function(e){
      return(cphases.melt[,e]*mass/1e6)
    })

    # Total mass of each LREE already afected
   lree.affected<-apply(lree.mass,2,sum)

    # LREE mass balance
    lree.excess <- c0[lree.nm]/1e4 - lree.affected

    # Stoechiometric proportion of each REE in the relevant end-member
    mol.mass.mnz<-MW.REE+mw["P"]+4*mw["O"]
    st<-MW.REE/mol.mass.mnz

    # Excess LREE are affected to Mnz
    mnz.mass<-lree.excess/st # relative to whole system (100g)

    # Forming Mnz requires other elements !
    p.used<- mnz.mass*mw["P"]/mol.mass.mnz # in g per 100 g of system
    o.used <- 4*mnz.mass*mw["O"]/mol.mass.mnz

    # We take a bit of that O2 from everywhere - no partitionning here
    # We do nothing to P
    o.mass<-mass/sum(mass)*sum(o.used)
    p.mass<-mass/sum(mass)*sum(p.used)

    # Final mass of the system
    # It is heavier by LREE and P that we added out of nowhere...
    mass<-c(mass-o.mass,sum(mnz.mass))
    names(mass)[length(mass)]<-"Mnz"

    # Expressed as proportions
    foo<-mass[c(rownames(cmins),"Mnz")]
    m.pr<-foo/sum(foo)

    # Mnz composition, in ppm
    cmnz<-lree.excess/mass["Mnz"]*1e6

    # (apparent) Kd of monazite
    kd["Mnz",lree.nm]<-cmnz/lree.c ###

    # As for Zrn (see comments), three types of elements
    # 1. ordinary elements - we do partitioning with the new min.props
    ee<-BatchPM(kd=kd,c0=c0,pm=pm,
                min.props=m.pr)

    cL<-ee$cL
    cm<-ee$cmins

    # 2. "dont" elements - we inherit the values (except for Mnz)
    cL[dont] <- melt.arg$trc[dont]
    cm[all.but.mnz,dont]<-cmins[,dont]

    # 3. Mnz forming elements - special case
    cL[lree.nm] <- lree.c
    cm["Mnz",lree.nm]<-cmnz
    cm[all.but.mnz,lree.nm]<-cphases

  }
  # D is an "observed" D that we calculate from cS/cL, in all cases
  cS<-.sanitize(m.pr%*%cm)
  DD<-cS/cL

  dont<-c(dont,lree.nm)

  invisible(list(c0=c0,cL=.sanitize(cL),cmins=cm,cS=cS,min.props=m.pr,FF=pm,kd=kd,DD=.sanitize(DD),
                 dont=dont,tot.mass=tot.mass,corr.masses=mass,o.used=sum(o.used),p.used=sum(p.used)))

}


correctZrnMnzSat<-function(kd,c0,pm,cmins,min.props,melt.arg=list(),dont=character(0)){
  #' @export

  zz<-correctZrnSat(kd=kd,
                    c0=c0,
                    pm=pm,
                    min.props=min.props,
                    melt.arg=melt.arg,
                    cmins=cmins,dont=dont)
  mm<-correctMnzSat(kd=kd,
                    c0=c0,
                    pm=pm,
                    min.props=zz$min.props,
                    melt.arg=list(TT=melt.arg$TT,
                                  mjrs=melt.arg$mjrs,
                                  H2O=melt.arg$H2O,
                                  trc=zz$cL
                    ),
                    cmins=zz$cmins,dont=zz$dont)

   invisible(list(cL=mm$cL,
                  cmins=mm$cmins,
                  kd=mm$kd,
                  DD=.sanitize(mm$DD),
                  min.props=mm$min.props,
                  dont=mm$dont))

  # ee<-correctZrnSat(c0,cL=cL,cmins=cmins,kd=kd,pm=pm,min.props=min.props,melt.arg=melt.arg,dont=dont)
  # ee<-.correctMnzSat(c0=c0,cL=ee$cL,cmins=ee$cmins,kd=ee$kd,pm=pm,min.props=ee$min.props,melt.arg=melt.arg,dont=ee$dont)
  # #ee<-.correctMnzSat.Rc(c0,cL=cL,cmins=cmins,kd=kd,pm=pm,min.props=min.props,melt.arg=melt.arg,dont=dont)
  #invisible(list(cL=ee$cL,cmins=ee$cmins,kd=ee$kd,min.props=ee$min.props))
}



############ PRIVATE ######################

###############################
# Clean a one-line matrix or data.frame, turn it into a vector
# And normalize, if needed
################################

.sanitize<-function(obj,normalize.TO=0){
  #' Sanitize objects that should be one (numeric) line
  #' typically mineral proportions or C0
  #' Make sure it's a vector
  #' Remove NA's
  #' Normalize to specified value
  #' @param obj the one-line thing to be trated
  #' @param normalize.TO new total, if =0 then no normalization occurs
  #' @return a VECTOR with no empty element, and (possibly) normalized to the desired value

  # Convert one-line matrix to vector
  if(class(obj)=="matrix"){
    nn<-colnames(obj)
    obj<-as.vector(obj)
    names(obj)<-nn
  }

  # Idem for data.frame
  if(class(obj)=="data.frame"){
    nn<-colnames(obj)
    obj<-as.numeric(obj)
    names(obj)<-nn
  }

  # Protect against empty values
  mis<-is.na(obj)
  obj<-obj[!mis]

  # Normalize
  if(normalize.TO>0)
  obj<-obj/sum(obj)*normalize.TO

  return(obj)
}

###############################
# Revert all-caps oxides to their true names
################################
#Convert upper case to mixed case for major elements
.TrueNames<-function(x){

  mix_case<-c("SiO2","K2O","TiO2","FeO","Fe2O3","FeOt","Fe2O3t","O2","CaO","MgO","Na2O","Al2O3","H2O","MnO","NiO","ZrO2","Cl2","CO2","CuO","Cr2O3","S2","F2")
  converted_case<-NULL
  for(i in 1:length(x)){
    #browser()
    converted_case<-c(converted_case,mix_case[which(toupper(mix_case)==toupper(x[i]))])
  }
  return(converted_case)
}


.TrueNamesOld<-function(ee){
  nm<-names(ee)
  tn<-c("SiO2","TiO2","Al2O3","Fe2O3","FeO","FeOt","MnO","MgO","CaO","Na2O","K2O","H2O.PLUS","CO2","P2O5","F","S")
  ucn<-toupper(tn)

  needed.names<-ucn %in% nm

  corr.names<-tn
  corr.names[needed.names]<-

  names(ee)<-corr.names
  return(ee)
}

###############################
# Clean millications
# GCDkit's function does not gurantee that all variants of Fr will be returned
# (because GCDkit knows that WR is clean - we don't)
################################

.cleanMillicats<-function(milli){
  # Make sure we have a matrix all the time
  if(is.vector(milli)){
    milli<-t(as.matrix(milli))
    vec<-T
  }else{
    vec<-F
  }

  oxd<-colnames(milli)
  if(!("FeOt" %in% oxd)){
    # FeOt is missing, we must do something !
    if(any(colnames(milli)=="FeO")){fe2<-milli[,"FeO"]}else{fe2<-0}
    if(any(colnames(milli)=="Fe2O3")){fe3<-milli[,"Fe2O3"]}else{fe3<-0}
    FeOt<-fe2+fe3
    milli<-cbind(milli,FeOt)
  }

  # Return what we have been given...
  if(vec){
    milli<-milli[1,]
  }

  return(milli)
}

###############################
# A slightly modified version of Monazite sat, taking Th into account (and returning more info).
# Also added some drop=F to make it work with one line
################################

mzSaturationWithTh<-function(cats=milli,REE=filterOut(WR,c("Th","La","Ce","Pr","Nd","Sm","Gd"),1),
                             H2O=3,Xmz=0,Temp=750){
  #' @export
  on.exit(options("show.error.messages"=TRUE))
  if(Temp==0){
    Temp<-winDialogString("Temperature (degrees C)","750")
    if(is.null(Temp)){
      cat("Cancelled.\n")
      options("show.error.messages"=FALSE)
      stop()
    }
  }

  Temp<-as.numeric(Temp)+273.15
  if(is.na(H2O))H2O<-0.0000001

  if(all(H2O==0)){
    H2O<-winDialogString("Water contents in the melt (%)","3")
    if(is.null(H2O)){cat("Cancelled.\n");options("show.error.messages"=FALSE);stop()}
    H2O<-as.numeric(H2O)
  }

  if(Xmz==0){
    Xmz<-winDialogString("X mz","0.99") # Default value assumes 1% U in monazite
    if(is.null(Xmz)){cat("Cancelled.\n");options("show.error.messages"=FALSE);stop()}
  }
  Xmz<-as.numeric(Xmz)
  if(nrow(cats)>1) cats<-cats[rownames(REE),,drop=F]

  lree.nm<-c("Th","La","Ce","Pr","Nd","Sm","Gd")
  # Quick fix in the case of missing GCDkit
  mw<-c(232.04,138.91,140.12,140.91,144.24,150.3,157.2,30.974,15.999)
  names(mw)<-c(lree.nm,"P","O")

  MW.REE<-mw[lree.nm]
  ree<-t(t(REE)/MW.REE)
  ree<-apply(ree,1,sum,na.rm=TRUE)
  reex<-ree/Xmz

  sums<-apply(cats,1,sum,na.rm=TRUE)
  ee<-sapply(1:nrow(cats),function(i){
    z<-cats[i,]/sums[i]
    return(z)
  })
  x<-t(ee)

  D.M<-(x[,"Na2O"]+x[,"K2O"]+2*x[,"CaO"])/x[,"Al2O3"]*1/(x[,"Al2O3"]+x[,"SiO2"])*1

  T.calc<-13318/(9.5+2.34*D.M+0.3879*sqrt(H2O)-log(reex))-273.15

  # Difference with mzSaturation: here we also return the sat value of REE
  sum.ree.sat<-exp(9.5+2.34*D.M+0.3879*sqrt(H2O)-13318/Temp)

  y<-cbind(D.M,round(T.calc,3),sum.ree.sat)

  colnames(y)<-c("Dmz","Tmz.sat.C","Sum.A.REE.sat")
 # if(!getOption("gcd.shut.up"))print(y)
  assign("results",y,.GlobalEnv)
  invisible(y)
}


############ NOT USED ? ###################

kd.table<-function(file,data=WR,na.iszero=F){
  #' Read a tab-separated file containing Kd values and convert it into a properly
  #' formatted table
  #' @param file Fully qualified (i.e. with path) name of the file to read (use
  #'   \code{setwd()} as appropriate).
  #' @param data Matrix or data frame containing the "reference" data: only
  #'   elements presnet both in \code{data} and \code{file} will be kept.
  #' @param na.iszero If \code{TRUE}, replace missing values by zeros. Alters the
  #'   original data, of course, but allows calculations even when some values are
  #'   missing.
  #' @return a data.frame containing the Kd values

  cat("WARNING - if you see this message you are using a deprecated function (kd.table)\n")


  # Read Kd file, should be "properly" formatted
  ee<-read.table(file,sep="\t")

  # Test whether the data table exists
  # If it does, filter the kd table to keep only usable elements

  cn<-try(colnames(data),silent=T)
  if(!class(cn)=="try-error"){
    cn<-intersect(cn,colnames(ee))
    ee<-ee[,cn,drop=F]
  }

  # If required, replace NA's by zeros
  if(na.iszero){
    ee[is.na(ee)]<-0
  }

  return(as.matrix(ee))
}

############ ORIGINAL FUNCTIONS ###########
# DO NOT CHANGE

# Loading Kd file
.loadKdTableOriginal<-function(){
  #' @export

  cat("This is .loadKdTableOriginal\n")
    on.exit(options("show.error.messages"=TRUE))
    filename<-choose.files(caption = "Select Kd file",multi = FALSE)
    kd<-read.table(filename,sep="\t")
    kd<-as.matrix(kd)
    assign("kd.table.bak",kd,.GlobalEnv) # Keep a backup of complete Kd table
    kd.new<-matrix(NA,nrow=length(properPhases),ncol=ncol(kd))
    rownames(kd.new)<-properPhases
    colnames(kd.new)<-colnames(kd)
    replace.Kd<-function(old.name,new.name){
        on.exit(options("show.error.messages"=TRUE))
        if(!old.name%in%rownames(kd)){
            message<-paste("Phase ",old.name," not found in the kd file!\nReplacing by zeros.....",sep="\n")
            winDialog(type="ok",message)
            options(show.error.messages=FALSE)
            kd.new[new.name,]<<-rep(0,ncol(kd.new)) # NEW
            return()
        }
        if(!new.name%in%allPhases){
            message<-paste("Phase ",new.name," not found in the Perplex file!",sep="\n")
            winDialog(type="ok",message)
            options(show.error.messages=FALSE)
            return()
        }
        kd.new[new.name,]<<-kd[old.name,]
    }
    ee<-lapply(properPhases,function(i){
        old<-gsub("_[0-9]{1}$","",i)
        new<-i
        replace.Kd(old,new)
    })
    return(kd.new)
}

# Loading Kd file
.loadKdTable<-function(){
  #' @export

  cat("This is .loadKdTable\n")

  on.exit(options("show.error.messages"=TRUE))
  filename<-choose.files(caption = "Select Kd file",multi = FALSE)
  kd<-read.table(filename,sep="\t")
  kd<-as.matrix(kd)
  assign("kd.table.bak",kd,.GlobalEnv) # Keep a backup of complete Kd table
  kd.new<-matrix(NA,nrow=length(properPhases),ncol=ncol(kd))
  rownames(kd.new)<-properPhases
  colnames(kd.new)<-colnames(kd)
  replace.Kd<-function(old.name,new.name){
    on.exit(options("show.error.messages"=TRUE))
    if(!old.name%in%rownames(kd)){
      message<-paste("Phase ",old.name," not found in the kd file!\nReplacing by zeros.....",sep="\n")
      winDialog(type="ok",message)
      options(show.error.messages=FALSE)
      kd.new[new.name,]<<-rep(0,ncol(kd.new)) # NEW
      return()
    }
    if(!new.name%in%allPhases){
      message<-paste("Phase ",new.name," not found in the Perplex file!",sep="\n")
      winDialog(type="ok",message)
      options(show.error.messages=FALSE)
      return()
    }
    kd.new[new.name,]<<-kd[old.name,]
  }
  ee<-lapply(properPhases,function(i){
    old<-gsub("_[0-9]{1}$","",i)
    new<-i
    replace.Kd(old,new)
  })
  return(kd.new)
}

# Specify melting module, can be replaced by another one under different name
.BatchPMOriginal<-function(kd,c0,pm,min.props,melt.arg){
  #' @export

        which.elems<-names(c0)
        DD<-min.props%*%kd
        FF<-pm/100
        cL<-c0/(DD+FF*(1-DD))
        ee<-sapply(names(min.props),function(i){
            z<-cL[1,which.elems]*kd[i,which.elems]
            return(z)
        },simplify=TRUE)
        cmins<-t(ee)
        cS<-cL[1,which.elems]*DD[1,which.elems]
        return(list(cL=cL,cmins=cmins,cS=cS,kd=kd))
}

.BatchPM<-function(kd,c0,pm,min.props,melt.arg){
  #' @export


  which.elems<-names(c0)
  DD<-min.props%*%kd
  FF<-pm/100
  cL<-c0/(DD+FF*(1-DD))
  ee<-sapply(names(min.props),function(i){
    z<-cL[1,which.elems]*kd[i,which.elems]
    return(z)
  },simplify=TRUE)
  cmins<-t(ee)
  cS<-cL[1,which.elems]*DD[1,which.elems]
  return(list(cL=cL,cmins=cmins,cS=cS,kd=kd))
}

.correctZrnSat<-function(c0,cL,cmins,kd,pm,min.props,melt.arg,dont=character(0)){
  #' @export
   DD<-min.props%*%kd
        which.elems<-names(c0)
        which.elems<-setdiff(which.elems,dont)

        milcats <- millications(melt.arg$mjrs)
        # call Saturation plugin, Zr just an arbitrary valueas it is not needed
        sat.data <- zrSaturation(milcats,T=melt.arg$TT-273.15,Zr=100)
        Zr.sat <- sat.data["Zr.sat"]

    # Add Kds for zircon to the table
    kd.full<-rbind(kd,kd.table.bak["Zrn",colnames(kd)])
    rownames(kd.full)[nrow(kd.full)]<-"Zrn"

    m.pr <- min.props
    # Correction for Zrn
    if(cL[,"Zr"]<=Zr.sat){
        # The melt is not saturated in zircon, do nothing more!
        zrc.prop <- 0
        m.pr["Zrn"] <-0
        cS<-cL[1,which.elems]*DD[1,which.elems]
    }else{
        # The melt is saturated, we correct for saturation amount of zircon
        zr.excess <- cL[,"Zr"]-Zr.sat
        zrc.prop <- zr.excess/497644
        m.pr["Zrn"] <- zrc.prop
        m.pr <- m.pr/sum(m.pr)


        D.c <- m.pr%*%kd.full
        cL[,which.elems] <- c0[which.elems]/(D.c[,which.elems]+pm/100*(1-D.c[,which.elems]))
        cL[,"Zr"] <- Zr.sat

        #cS<-cL[1,which.elems]*D.c[1,which.elems]
        #cS["Zr"] <- (c0["Zr"]-pm/100*cL[,"Zr"])/(1-pm/100)
    }
        # All minerals compositions
        ee<-sapply(names(m.pr),function(i){
            z<-cL[1,which.elems]*kd.full[i,which.elems]
            return(z)
        },simplify=TRUE)
        prd<-cmins[,dont,drop=FALSE]
        prd<-rbind(prd,NA)
        rownames(prd)[length(prd)]<-"Zrn"
        cmins<-cbind(t(ee),prd)
        cmins["Zrn","Zr"]<-497644
        dont<-c(dont,"Zr")
        invisible(list(cL=cL,cmins=cmins,kd=kd.full,min.props=m.pr,dont=dont))
}

# Specify the correction module, this one accounts for monazite saturation
.correctMnzSat<-function(c0,cL,cmins,kd,pm,min.props,melt.arg,dont=character(0)){
  #' @export
   DD<-min.props%*%kd
        which.elems<-names(c0)
        which.elems<-setdiff(which.elems,dont)

        ##Definitions, variables
        lree.nm<-c("Th","La","Ce","Pr","Nd","Sm","Gd")
        W.R<-c(-520,-260,-116,31,177,460,750)
        names(W.R)<-lree.nm
        MW.REE<-mw[lree.nm]
        REE<-cL[1,lree.nm,drop=FALSE]
        M<-rep(0,length(lree.nm))

        #Xmz<-0.99 # Assume 1% of U in monazite
        Xmz<-1

        milcats <- millications(melt.arg$mjrs)

        ## Actual amount of LREE in the melt (uncorrected)
        S<-REE[lree.nm]/MW.REE[lree.nm]
        ST<-sum(S)

        # Call the Saturation plugin and get the saturation value of LREE + Th
        ee <- mzSaturationWithTh(cats=milcats,REE=REE,H2O=melt.arg$H2O,Xmz=Xmz,Temp=melt.arg$TT-273.15) # Temperature is in C
        AiT<-ee[1,"Sum.A.REE.sat"]

        # Add Kds for monazite to the table
        kd.full<-rbind(kd,kd.table.bak["Mnz",colnames(kd)])
        rownames(kd.full)[nrow(kd.full)]<-"Mnz"

        m.pr <- min.props
        # Correction for Mnz
        if(ST<=AiT){
            # The melt is not saturated in monazite, do nothing more!
            mnz.prop<-0
            m.pr["Mnz"] <- 0
        }else{
            # The melt is saturated, we calculate the mnz composition

            # Partitioning equation: gamma=(Xmelt/XMnz)
            gamma<-exp(W.R/(melt.arg$TT))

            # We do not know MiT and A, will optimize the equation for MiT
            .compute<-function(MiT){
                A<-gamma*AiT*S/(MiT+gamma*AiT)
                ee<-(sum(A)-AiT)^2
                return(ee)
            }

            m.guess<-0
            out<-optim(m.guess,.compute)
            R2<-out$value
            MiT<-out$par

            # Proportions in the liquid
            A<-gamma*AiT*S/(MiT+gamma*AiT)
            # Proportions in monazite
            M<-S-A

            mf.mnz<-sum((MW.REE+mw["P"]+4*mw["O"])*M)/1e6
            mnz.prop<-mf.mnz*pm/100

            # Concentrations of LREE + Th in the melt
            lree.c<-A*MW.REE

            # The melt is saturated, we correct for saturation amount of monazite
            m.pr <- min.props
            m.pr["Mnz"] <- mnz.prop
            m.pr <- m.pr/sum(m.pr)
            D.c <- m.pr%*%kd.full
            cL[,which.elems] <- c0[which.elems]/(D.c[,which.elems]+pm/100*(1-D.c[,which.elems]))
            cL[,lree.nm] <- lree.c
            #cS<-cL[1,which.elems]*D.c[1,which.elems]
            #cS[lree.nm] <- (c0[lree.nm]-pm/100*cL[,lree.nm])/(1-pm/100)

    }
        ee<-sapply(names(m.pr),function(i){
            z<-cL[1,which.elems]*kd.full[i,which.elems]
            return(z)
        },simplify=TRUE)

        # All minerals composition
        prd<-cmins[,dont,drop=FALSE]
        prd<-rbind(prd,NA)
        rownames(prd)[length(prd)]<-"Mnz"
        cmins<-cbind(t(ee),prd)

        # Mnz composition
        x<-M/sum(M)
        MW.mnz<-sum(x*MW.REE)+mw["P"]+4*mw["O"]
        ci<-x*MW.REE/MW.mnz*1e6
        cmins["Mnz",lree.nm]<-ci
        dont<-c(dont,lree.nm)
        invisible(list(cL=cL,cmins=cmins,kd=kd,min.props=m.pr,dont=dont))
}

# Combined, sequential, first Zrn then Mnz
.correctZrnMnzSat<-function(c0,cL,cmins,kd,pm,min.props,melt.arg,dont=character(0)){
  #' @export
  ee<-.correctZrnSat(c0,cL=cL,cmins=cmins,kd=kd,pm=pm,min.props=min.props,melt.arg=melt.arg,dont=dont)
    ee<-.correctMnzSat(c0=c0,cL=ee$cL,cmins=ee$cmins,kd=ee$kd,pm=pm,min.props=ee$min.props,melt.arg=melt.arg,dont=ee$dont)
    #ee<-.correctMnzSat.Rc(c0,cL=cL,cmins=cmins,kd=kd,pm=pm,min.props=min.props,melt.arg=melt.arg,dont=dont)
    invisible(list(cL=ee$cL,cmins=ee$cmins,kd=ee$kd,min.props=ee$min.props))
}


########## BACKUP COPIES
## FIXME return NA for cmins where a mineral is not present !
BatchPMOKZ<-function(kd,c0,pm,cmins=matrix(),min.props,melt.arg=list(),dont=character(0)){

  #' Trace elements concentrations for a batch process
  #'
  #' More commonly used for batch melting, but does actually also work with crystalization.
  #' Make sure to use the same mineral/element names in all the variables
  #' @param kd Named matrix. Partition coefficients
  #' @param c0 Named vector. The composition of the source
  #' @param pm Scalar. Degree of partial melting (0-100)
  #' @param cmins Matrix. Composition of minerals (not used, for consistency)
  #' @param min.props Named vector. Mineral proportions (normalized to 1)
  #' @param melt.arg List. Further useful arguments (mostly used for correction modules, kept here for consistency)
  #' @param dont Character vectors. Do not calculate for these elements. Not used here, kept for consistency (mostly used in correction modules)
  #' @return list containing the following slots:
  #' \itemize{
  #'   \item \code{cL} Named vector. Composition of the liquid
  #'   \item \code{cmins} Named matrix. Composition of all individual solid mineral phases
  #'   \item \code{cS} Named vector. Composition of the whole solid
  #'   \item \code{kd} Named matrix. The Kd table, here similar to the input (included for consistency with other modelling routines)
  #'   \item \code{DD} Named vector. The D (bulk distribution) values.
  #' }

  c0<-.sanitize(c0)
  min.props<-.sanitize(min.props,normalize.TO=1)

  ## TODO c0 or min.props are data.frame !!!
  # Protect against c0 and min.props being formatted as vector or matrix
  # if(class(c0)=="matrix"){
  #   nn<-colnames(c0)
  #   c0<-as.vector(c0)
  #   names(c0)<-nn
  # }
  #
  # # # Protect against c0 and min.props being formatted as vector or matrix
  # # if(class(min.props)=="matrix"){
  # #   nn<-colnames(min.props)
  # #   min.props<-as.vector(min.props)
  # #   names(min.props)<-nn
  # # }
  #
  # # Protect against empty C0 values
  # mis<-is.na(c0)
  # c0<-c0[!mis]
  #
  # # # Protect against empty min.props values
  # # mis<-is.na(min.props)
  # # min.props<-min.props[!mis]
  #
  # # # Normalize min props
  # # min.props<-min.props/sum(min.props)

  which.elems<-intersect(names(c0),colnames(kd))
  which.mins<-names(min.props)

  # Test whether Kds for all minerals are available
  missing<-setdiff(names(min.props),rownames(kd))
  if(length(missing)!=0){
    stop(paste("Missing Kd value for mineral(s)",missing,"\n"))
  }

  # Keep only the required elements/minerals
  # drop=F is important to preserve the objects as matrices even if only one line
  #c0<-c0[which.elems,drop=F]
  c0<-c0[which.elems] # c0 has been sanitized to a vector by then
  kd<-kd[which.mins,which.elems,drop=F]

  # Do the real job...
  DD<-min.props%*%kd
  FF<-pm/100
  cL<-c0/(DD+FF*(1-DD))
  ee<-sapply(names(min.props),function(i){
    z<-cL[1,which.elems]*kd[i,which.elems]
    return(z)
  },simplify=TRUE)
  cmins<-t(ee) # TODO subset to keep only the ones that exist
  cS<-cL[1,which.elems]*DD[1,which.elems]
  return(list(cL=.sanitize(cL),cmins=cmins,cS=cS,min.props=min.props,kd=kd,DD=.sanitize(DD)))
}

correctZrnSatOKZ<-function(kd,c0,pm,cmins=matrix(),min.props,melt.arg=list(),dont=character(0)){

  #' Correction for zircon saturation, Watson & Harrison 1983
  #' @param kd Named matrix. Partition coefficients. Must include Zircon (Zrn)
  #' @param pm Scalar. Degree of partial melting (0-100)
  #' @param cmins Matrix. Composition of minerals (not used, for consistency)
  #' @param min.props Named vector. Mineral proportions (normalized to 1)
  #' @param melt.arg List. Further useful arguments.
  #' Includes here the melt major elements (melt.arg$mjrs),trace (melt.arg$trc), temperature in KELVIN (melt.arg$TT)
  #' @param dont Character vector. Do not calculate for these elements. Use to avoid overptinting previous saturation calc (i.e. if you use a saturation routine /after/ this one, )

  c0<-.sanitize(c0)
  m.pr<-.sanitize(min.props,normalize.TO=1)

  # DD<-min.props%*%kd
  # which.elems<-names(c0) #### CAREFUL, also cross ref with kd
  # which.elems<-setdiff(which.elems,dont)

  # Need more protection, cf. BatchPM

  # GCDkit must be loaded & plugins sourced for the following to work
  # .. and loadData must have operated, to get access to MW etc.
  # Silent output for the plugin ??

  milcats <- millications(melt.arg$mjrs)
  # call Saturation plugin, Zr just an arbitrary value as it is not needed
  sat.data <- zrSaturation(milcats,T=melt.arg$TT-273.15,Zr=100)
  Zr.sat <- sat.data["Zr.sat"]

  # # Add Kds for zircon to the table
  # kd.full<-rbind(kd,kd.table.bak["Zrn",colnames(kd)]) ##### CAREFUL
  # rownames(kd.full)[nrow(kd.full)]<-"Zrn"


  # Protect against empty min.props values
  # mis<-is.na(min.props)
  # min.props<-min.props[!mis]


  # Correction for Zrn
  if(melt.arg$trc["Zr"]<=Zr.sat){
    # The melt is not saturated in zircon, do nothing more!
    zrc.prop <- 0
    m.pr["Zrn"] <-0
    #cS<-cL[1,which.elems]*DD[1,which.elems]
    cL<-melt.arg$trc
    # cmins are inherited
  }else{
    # The melt is saturated, we correct for saturation amount of zircon

    zr.excess <- melt.arg$trc["Zr"]-Zr.sat
    zrc.prop <- zr.excess/497644
    m.pr["Zrn"] <- zrc.prop
    m.pr <- m.pr/sum(m.pr)
    # We must now recalculate the melt compo taking the Zrn into account
    ee<-BatchPM(kd=kd,c0=c0,pm=pm,
                min.props=m.pr)
    cL<-ee$cL
    cL["Zr"] <- Zr.sat
    cmins<-ee$cmins
    cmins["Zrn","Zr"]<-497644
  }




  # D.c <- m.pr%*%kd.full
  # cL[,which.elems] <- c0[which.elems]/(D.c[,which.elems]+pm/100*(1-D.c[,which.elems]))
  # ... except for the elements that we explicitly do not touch !
  cL[dont] <- melt.arg$trc[dont] #### and cmins ?!!
  #cS<-cL[1,which.elems]*D.c[1,which.elems]
  #cS["Zr"] <- (c0["Zr"]-pm/100*cL[,"Zr"])/(1-pm/100)

  # All minerals compositions
  # ee<-sapply(names(m.pr),function(i){
  #   z<-cL[1,which.elems]*kd.full[i,which.elems]
  #   return(z)
  # },simplify=TRUE)
  # prd<-cmins[,dont,drop=FALSE]
  # prd<-rbind(prd,NA)
  # rownames(prd)[length(prd)]<-"Zrn"
  # cmins<-cbind(t(ee),prd)

  #cmins[,dont] <- melt.arg$trc[dont] ### TODO


  dont<-c(dont,"Zr")
  invisible(list(cL=.sanitize(cL),cmins=cmins,kd=kd,DD=numeric(0),min.props=m.pr,dont=dont))
  ### MAKE SURE WE RETURN VECTORS & not 1 LINE MATRIX
  ### DD ?
}

correctMnzSatOKZ<-function(kd,c0,pm,cmins=matrix(),min.props,melt.arg=list(),dont=character(0)){

  #' Correction for monazite saturation, Montel 1992 w. elements partitionning
  #' @param kd Named matrix. Partition coefficients. Must include Monazite (Mnz)
  #' @param pm Scalar. Degree of partial melting (0-100)
  #' @param cmins Matrix. Composition of minerals (not used, for consistency)
  #' @param min.props Named vector. Mineral proportions (normalized to 1)
  #' @param melt.arg List. Further useful arguments.
  #' Includes here the melt major elements (melt.arg$mjrs),trace (melt.arg$trc), water (melt.arg$H2O), temperature in KELVIN (melt.arg$TT)
  #' @param dont Character vector. Do not calculate for these elements. Use to avoid overptinting previous saturation calc (i.e. if you use a saturation routine /after/ this one, )


  c0<-.sanitize(c0)
  m.pr<-.sanitize(min.props,normalize.TO=1)

  #DD<-min.props%*%kd
  which.elems<-names(c0)
  which.elems<-setdiff(which.elems,dont)

  ##Definitions, variables
  lree.nm<-c("Th","La","Ce","Pr","Nd","Sm","Gd")
  W.R<-c(-520,-260,-116,31,177,460,750)
  names(W.R)<-lree.nm
  MW.REE<-mw[lree.nm]
  REE<-melt.arg$trc[lree.nm]
  M<-rep(0,length(lree.nm))

  #Xmz<-0.99 # Assume 1% of U in monazite
  Xmz<-1

  milcats <- millications(melt.arg$mjrs)

  ## Actual amount of LREE in the melt (uncorrected)
  S<-REE[lree.nm]/MW.REE[lree.nm]
  ST<-sum(S)

  # Call the Saturation plugin and get the saturation value of LREE + Th
  ee <- mzSaturationWithTh(cats=milcats,REE=REE,H2O=melt.arg$H2O,Xmz=Xmz,Temp=melt.arg$TT-273.15) # Temperature is in C
  AiT<-ee[1,"Sum.A.REE.sat"]

  # # Add Kds for monazite to the table
  # kd.full<-rbind(kd,kd.table.bak["Mnz",colnames(kd)])
  # rownames(kd.full)[nrow(kd.full)]<-"Mnz"

  # m.pr <- min.props
  # Correction for Mnz
  if(ST<=AiT){
    # The melt is not saturated in monazite, do nothing more!
    mnz.prop<-0
    m.pr["Mnz"] <- 0
    # The compo of the Mnz is not determined (!)
    ci<-rep(NA,length(lree.nm))
    names(ci)<-lree.nm
    lree.c<-rep(NA,length(lree.nm))
    names(lree.c)<-lree.nm
    cL<-melt.arg$trc
    # cmins<-cmins ## inherited

  }else{
    # The melt is saturated, we calculate the mnz composition

    # Partitioning equation: gamma=(Xmelt/XMnz)
    gamma<-exp(W.R/(melt.arg$TT))

    # We do not know MiT and A, will optimize the equation for MiT
    .compute<-function(MiT){
      A<-gamma*AiT*S/(MiT+gamma*AiT)
      ee<-(sum(A)-AiT)^2
      return(ee)
    }

    m.guess<-0
    out<-optim(m.guess,.compute)
    R2<-out$value
    MiT<-out$par

    # Proportions in the liquid
    A<-gamma*AiT*S/(MiT+gamma*AiT)
    # Proportions in monazite
    M<-S-A

    mf.mnz<-sum((MW.REE+mw["P"]+4*mw["O"])*M)/1e6
    mnz.prop<-mf.mnz*pm/100

    # Concentrations of LREE + Th in the melt
    lree.c<-A*MW.REE

    # The melt is saturated, we correct for saturation amount of monazite
    m.pr <- min.props
    m.pr["Mnz"] <- mnz.prop
    m.pr <- .sanitize(m.pr,1)

    # D.c <- m.pr%*%kd.full
    # cL[,which.elems] <- c0[which.elems]/(D.c[,which.elems]+pm/100*(1-D.c[,which.elems]))
    # cL[,lree.nm] <- lree.c
    # #cS<-cL[1,which.elems]*D.c[1,which.elems]
    #cS[lree.nm] <- (c0[lree.nm]-pm/100*cL[,lree.nm])/(1-pm/100)
    # Mnz composition
    x<-M/sum(M)
    MW.mnz<-sum(x*MW.REE)+mw["P"]+4*mw["O"]
    ci<-x*MW.REE/MW.mnz*1e6
    # We must now recalculate the melt & mins compo
    ee<-BatchPM(kd=kd,c0=c0,pm=pm,
                min.props=m.pr)
    cL<-ee$cL
    # Melt compo is the saturation value
    cL[lree.nm] <- lree.c

    cmins<-ee$cmins
    # Mnz compo is determined by the above eqn, for the LREE
    cmins["Mnz",lree.nm]<-ci
  }


  # Whatever comes in "dont" must be preserved
  cL[dont] <- melt.arg$trc[dont] #### and cmins ?!!

  # ee<-sapply(names(m.pr),function(i){
  #   z<-cL[1,which.elems]*kd.full[i,which.elems]
  #   return(z)
  # },simplify=TRUE)
  #
  # All minerals composition
  # prd<-cmins[,dont,drop=FALSE]
  # prd<-rbind(prd,NA)
  # rownames(prd)[length(prd)]<-"Mnz"
  # cmins<-cbind(t(ee),prd)


  dont<-c(dont,lree.nm)
  invisible(list(cL=.sanitize(cL),cmins=cmins,kd=kd,DD=numeric(0),min.props=m.pr,dont=dont))


}


correctMnzSatOKZ2<-function(kd,c0,pm,cmins,min.props,melt.arg=list(),dont=character(0)){

  #' Correction for monazite saturation, Montel 1992 w. elements partitionning
  #' @param kd Named matrix. Partition coefficients. Must include Monazite (Mnz)
  #' @param pm Scalar. Degree of partial melting (0-100)
  #' @param cmins Matrix. Composition of minerals (not used, for consistency)
  #' @param min.props Named vector. Mineral proportions (normalized to 1)
  #' @param melt.arg List. Further useful arguments.
  #' Includes here the melt major elements (melt.arg$mjrs),trace (melt.arg$trc), water (melt.arg$H2O), temperature in KELVIN (melt.arg$TT)
  #' @param dont Character vector. Do not calculate for these elements. Use to avoid overptinting previous saturation calc (i.e. if you use a saturation routine /after/ this one, )

  #' @export

  c0<-.sanitize(c0)
  m.pr<-.sanitize(min.props,normalize.TO=1)

  #DD<-min.props%*%kd
  # which.elems<-names(c0)
  # which.elems<-setdiff(which.elems,dont)

  ##Definitions, variables
  lree.nm<-c("Th","La","Ce","Pr","Nd","Sm","Gd")

  # Quick fix in the case of missing GCDkit
  mw<-c(232.04,138.91,140.12,140.91,144.24,150.3,157.2,30.974,15.999)
  names(mw)<-c(lree.nm,"P","O")

  # Protection against Matt converting oxides to upper case !!
  names(melt.arg$mjrs)<-.TrueNames(names(melt.arg$mjrs))

  W.R<-c(-520,-260,-116,31,177,460,750)
  names(W.R)<-lree.nm
  MW.REE<-mw[lree.nm]
  REE<-melt.arg$trc[lree.nm]
  M<-rep(0,length(lree.nm))

  #Xmz<-0.99 # Assume 1% of U in monazite
  Xmz<-0.99

  milcats <- millications(melt.arg$mjrs)

  ## Actual amount of LREE in the melt, atomic (uncorrected)
  S<-REE[lree.nm]/MW.REE[lree.nm]
  ST<-sum(S)

  # Call the Saturation plugin and get the saturation value of LREE + Th
  ee <- mzSaturationWithTh(cats=milcats,REE=REE,H2O=melt.arg$H2O,Xmz=Xmz,Temp=melt.arg$TT-273.15) # Temperature is in C
  AiT<-ee[1,"Sum.A.REE.sat"]

  if(ST<=AiT){
    # The melt is not saturated in monazite, do nothing more!
    # mnz.prop<-0
    # m.pr["Mnz"] <- 0
    cm<-cmins # cmins are inherited, don't change
    cL<-melt.arg$trc

    # The compo of the Mnz is not determined (!)
    # ci<-rep(NA,length(lree.nm))
    # names(ci)<-lree.nm
    # lree.c<-rep(NA,length(lree.nm))
    # names(lree.c)<-lree.nm
    # cL<-melt.arg$trc
    # cmins<-cmins ## inherited

  }else{
    # The melt is saturated, we calculate the mnz composition

    # Partitioning equation: gamma=(Xmelt/XMnz)
    gamma<-exp(W.R/(melt.arg$TT))

    # We do not know MiT and A, will optimize the equation for MiT
    .compute<-function(MiT){
      A<-gamma*AiT*S/(MiT+gamma*AiT)
      ee<-(sum(A)-AiT)^2
      return(ee)
    }

    m.guess<-0
    out<-optim(m.guess,.compute,method="Brent",lower=0,upper=1)
    R2<-out$value
    MiT<-out$par

    # Proportions in the liquid
    A<-gamma*AiT*S/(MiT+gamma*AiT)
    # Proportions in monazite
    M<-S-A

    mf.mnz<-sum((MW.REE+mw["P"]+4*mw["O"])*M)/1e6
    mnz.prop<-mf.mnz*pm/100

    # Concentrations of LREE + Th in the melt
    lree.c<-A*MW.REE

    # Mnz composition
    x<-M/sum(M)
    MW.mnz<-sum(x*MW.REE)+mw["P"]+4*mw["O"]
    cmnz<-x*MW.REE/MW.mnz*1e6

    # (apparent) Kd of monazite
    kd["Mnz",lree.nm]<-cmnz/lree.c

    # Now we correct for saturation amount of monazite
    m.pr["Mnz"] <- mnz.prop
    m.pr <- m.pr/sum(m.pr)

    # As for Zrn (see comments), three types of elements
    # 1. ordinary elements - we do partitioning with the new min.props


    ee<-BatchPM(kd=kd,c0=c0,pm=pm,
                min.props=m.pr)

    cL<-ee$cL
    cm<-ee$cmins

    # 2. "dont" elements - we inherit the values (except for Mnz)
    cL[dont] <- melt.arg$trc[dont]
    all.but.mnz<-setdiff(rownames(cm),"Mnz")
    cm[all.but.mnz,dont]<-cmins[,dont]

    # 3. Mnz forming elements - special case
    cL[lree.nm] <- lree.c
    cm["Mnz",lree.nm]<-cmnz
  }
  # D is an "observed" D that we calculate from cS/cL, in all cases
  cS<-.sanitize(m.pr%*%cm)
  DD<-cS/cL


  #
  #
  #   # D.c <- m.pr%*%kd.full
  #   # cL[,which.elems] <- c0[which.elems]/(D.c[,which.elems]+pm/100*(1-D.c[,which.elems]))
  #   # cL[,lree.nm] <- lree.c
  #   # #cS<-cL[1,which.elems]*D.c[1,which.elems]
  #   #cS[lree.nm] <- (c0[lree.nm]-pm/100*cL[,lree.nm])/(1-pm/100)
  #   # Mnz composition
  #   x<-M/sum(M)
  #   MW.mnz<-sum(x*MW.REE)+mw["P"]+4*mw["O"]
  #   ci<-x*MW.REE/MW.mnz*1e6
  #   # We must now recalculate the melt & mins compo
  #   ee<-BatchPM(kd=kd,c0=c0,pm=pm,
  #               min.props=m.pr)
  #   cL<-ee$cL
  #   # Melt compo is the saturation value
  #   cL[lree.nm] <- lree.c
  #
  #   cmins<-ee$cmins
  #   # Mnz compo is determined by the above eqn, for the LREE
  #   cmins["Mnz",lree.nm]<-ci
  # }
  #
  #
  # # Whatever comes in "dont" must be preserved
  # cL[dont] <- melt.arg$trc[dont] #### and cmins ?!!
  #
  # # ee<-sapply(names(m.pr),function(i){
  # #   z<-cL[1,which.elems]*kd.full[i,which.elems]
  # #   return(z)
  # # },simplify=TRUE)
  # #
  # # All minerals composition
  # # prd<-cmins[,dont,drop=FALSE]
  # # prd<-rbind(prd,NA)
  # # rownames(prd)[length(prd)]<-"Mnz"
  # # cmins<-cbind(t(ee),prd)


  dont<-c(dont,lree.nm)
  invisible(list(cL=.sanitize(cL),cmins=cm,kd=kd,DD=DD,min.props=m.pr,dont=dont))

}

