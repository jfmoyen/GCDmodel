####### OnLoad #######
# The main purpose of this function is to check whether GCDkit is loaded.
# If so, then we overwrite our functions with GCDkit's variant
######################

.onLoad<-function(lib, pkg){
  on.exit(options("show.error.messages"=TRUE))
  options("warn"=-1)

  packageStartupMessage("GCDmodel module - 0.7.0", appendLF = TRUE)

  # If we do not have access to GCDkit functions, redefine them here
  #if(!exists("millications")){
  if(! "GCDkit" %in% .packages() ){
    packageStartupMessage("GCDkit not loaded, using local variant of millications",
                          appendLF = TRUE)

    ##Definitions, variables - these must go somewhere else
    # lree.nm<<-c("Th","La","Ce","Pr","Nd","Sm","Gd")
    # lree.nm.2 <- c(lree.nm,"P","O")
    #
    # mw <<- periodicTable[lree.nm.2,"AtomicMass",drop=T]
    # names(mw) <- lree.nm.2

   }else{
     millications <- GCDkit::millications
     zrSaturation <- GCDkit::zrSaturation

     packageStartupMessage("GCDkit present, using GCDkit version of millications",
                              appendLF = TRUE)}


  packageStartupMessage("...Done!", appendLF = TRUE)

  invisible()
}
