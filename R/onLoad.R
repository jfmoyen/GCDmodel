####### OnLoad #######
# The main purpoise of this function is to check whether GCDkit is loaded.
# If not, we redefine some of the functions required
######################

.onLoad<-function(lib, pkg){
  on.exit(options("show.error.messages"=TRUE))
  options("warn"=-1)

  packageStartupMessage("GCDmodel module - 0.5.9", appendLF = TRUE)

  # If we do not have access to GCDkit functions, redefine them here
  #if(!exists("millications")){
  if(! "GCDkit" %in% .packages() ){
    packageStartupMessage("GCDkit not loaded, using local variant of millications",
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

   }else{packageStartupMessage("GCDkit present, using GCDkit version of millications",
                              appendLF = TRUE)}


  packageStartupMessage("...Done!", appendLF = TRUE)
  #winMenuAdd("Plugins")
  invisible()
}
