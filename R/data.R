#' Periodic table of elements
#'
#' Properties of chemical elements
#'
#' @format
#' A data frame with 118 rows and 28 columns.
#' See documentation in https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee
#'
#' @source <https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee>
"periodicTable"

#' Oxides of petrological interest
#'
#' Properties of chemical elements
#'
#' @format
#' A data frame with 22 rows and 6 columns.
#' \describe{
#'   \item{OxideName}{Full formula of the oxide}
#'   \item{CationName}{}
#'   \item{NbCations, NbO}{Number of cations and oxygen atoms in the oxide}
#'   \item{CationMass}{Mass of the cation, taken from GCDmodel::periodicTable}
#'   \item{OxideMass}{Mass of the oxide, calculated from the above data}
#' }
#'
#' @source J.F. Moyen
"petroOxides"

