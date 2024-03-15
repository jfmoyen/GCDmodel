## code to prepare dataset
# This code should be run ONLY when the underlying csv, in data-raw, are updated

periodicTable <- read.csv("data-raw/PeriodicTable.csv",sep=";")
rownames(periodicTable) <- periodicTable$Symbol

usethis::use_data(periodicTable, overwrite = TRUE)

####

petroOxides <- read.csv("data-raw/PetroOxides.csv",sep=";")
rownames(petroOxides) <- petroOxides$OxideName

petroOxides$CationMass <- periodicTable[petroOxides$CationName,"AtomicMass"]
petroOxides$OxideMass <- petroOxides$CationMass * petroOxides$NbCations + periodicTable["O","AtomicMass"] * petroOxides$NbO

usethis::use_data(petroOxides, overwrite = TRUE)
