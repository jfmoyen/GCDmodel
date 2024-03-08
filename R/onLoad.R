####### OnLoad #######
# The main purpoise of this function is to check whether GCDkit is loaded.
# If not, we redefine some of the functions required
######################

.onLoad<-function(lib, pkg){
  on.exit(options("show.error.messages"=TRUE))
  options("warn"=-1)

  packageStartupMessage("GCDmodel module - 0.5.9", appendLF = TRUE)

  # If we do not have access to GCDkit functions, redefine them here
  if(!exists("millications")){

    packageStartupMessage("GCDkit not available, redefining key functions",
                          appendLF = TRUE)

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

  }else{packageStartupMessage("GCDkit present, using GCDkit version of millications",
                              appendLF = TRUE)}


  packageStartupMessage("...Done!", appendLF = TRUE)
  #winMenuAdd("Plugins")
  invisible()
}
