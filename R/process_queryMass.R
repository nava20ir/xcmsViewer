#' Annotate Features Based on MS1 Masses
#'
#' Performs MS1-based annotation of LC–MS features by matching measured m/z values 
#' to reference compound masses, taking into account common adducts for positive 
#' or negative ion mode.
#'
#' @param x An object containing LC–MS features (e.g., a `prunedXcmsSet` or similar) 
#'   with feature metadata accessible via `fData(x@featureSet)`.
#' @param mode Character, either `"pos"` or `"neg"`. Specifies ionization mode 
#'   for selecting appropriate adducts. Default is `"pos"`.
#' @param ref Data frame of reference compounds with at least a mass column. 
#'   Used to match measured features.
#' @param ppmtol Numeric, mass tolerance in parts per million (ppm) for matching features to reference masses. Default is 10.
#' @param ips Numeric between 0 and 1, minimum IPS (Intensity/Presence Score) threshold for adducts to consider. Default is 0.75.
#' @param fun_parallel Function for parallel processing (e.g., `parallel::mclapply`). Default is `parallel::mclapply`.
#' @param ... Additional arguments passed to the parallel function.
#'
#' @return A data frame with MS1 annotation results:
#' \describe{
#'   \item{FeatureID}{ID of the LC–MS feature.}
#'   \item{MatchedCompound}{Name of the matched reference compound.}
#'   \item{Adduct}{Adduct type assigned to the feature.}
#'   \item{MassDiff}{Mass difference between observed and theoretical m/z.}
#'   \item{Score}{Optional score reflecting quality of the match.}
#' }
#' The returned object has an attribute `ppmtol` storing the mass tolerance used.
#'
#' @details
#' 1. Loads the MAIT package adduct tables (`posAdducts` or `negAdducts`) according to `mode`.  
#' 2. Filters adducts based on the IPS threshold.  
#' 3. For each feature, `massQuery()` is called to match the measured m/z to all possible
#'    reference compound masses with selected adducts within the ppm tolerance.  
#' 4. Parallelization is supported via `fun_parallel`.
#'
#' @examples
#' \dontrun{
#' # x is a prunedXcmsSet object
#' # ref is a reference compound table with at least a 'mass' column
#' annotations <- annotateMS1(x, mode = "pos", ref = ref, ppmtol = 10, ips = 0.8)
#' head(annotations)
#' }
#'
#' @export
annotateMS1 <- function(
  x, mode = c("pos", "neg")[1], ref, ppmtol = 10, ips = 0.75, fun_parallel = parallel::mclapply, ...
) {
  features <- fData(x@featureSet)
  maitEnv <- environment()
  data("MAITtables", package = "MAIT", envir = maitEnv)
  if (mode == "pos") {
    at <- maitEnv$posAdducts    
  } else if (mode == "neg") {
    at <- maitEnv$negAdducts
  } else 
      stop("'mode' should be either pos or neg!")
  at <- data.frame(
    Adduct = as.character(at$name),
    MassDiff = at$massdiff,
    Nmol = as.numeric(at$nmol),
    IPS = as.numeric(at$ips),
    stringsAsFactors = FALSE)
  at <- at[which(at$IPS >= ips), ]
  
  ll <- fun_parallel(1:nrow(features), function(i) {
    massQuery(m = features$mzmed[i], ppmtol = ppmtol, refTab = ref, 
              addTable = at, ID = features$ID[i]) 
  }, ...)
  v <- do.call(rbind, ll)
  attr(v, "ppmtol") <- ppmtol
  v
}


#' Query Possible Metabolites by m/z
#'
#' Given a measured mass (m/z), this function queries a reference metabolite 
#' table and an adduct table to find all potential metabolite matches within 
#' a specified mass tolerance (ppm).
#'
#' @param m Numeric, the measured mass (m/z) of the feature.
#' @param ppmtol Numeric, mass tolerance in parts per million (ppm) for matching.
#' @param refTab Data frame of reference metabolites. Must contain at least the following columns:
#'   \describe{
#'     \item{name}{Metabolite name}
#'     \item{monoMass}{Monoisotopic molecular weight}
#'   }
#' @param addTable Data frame of adducts. Must contain at least:
#'   \describe{
#'     \item{Adduct}{Name of the adduct, e.g., "[M+H]+"}
#'     \item{MassDiff}{Mass difference of the adduct relative to the neutral molecule}
#'     \item{Nmol}{Number of molecules in the adduct (e.g., 2 for [2M-H])}
#'     \item{IPS}{Intensity/Presence Score; possible values: 0.25, 0.5, 0.75, 1}
#'   }
#' @param ID Optional character or numeric ID to be included in the output for each matched row.
#'
#' @return A data frame of all reference metabolites that match the measured mass, including:
#' \describe{
#'   \item{ID}{Optional feature ID if provided}
#'   \item{InChIKey}{Metabolite InChIKey}
#'   \item{CID}{Compound ID from reference table}
#'   \item{cpdName}{Metabolite name}
#'   \item{formula}{Chemical formula}
#'   \item{monoMass}{Monoisotopic mass}
#'   \item{Adduct}{Matched adduct}
#'   \item{MassDiff}{Mass difference of the adduct}
#'   \item{Nmol}{Number of molecules in the adduct}
#'   \item{IPS}{Intensity/Presence Score for the adduct}
#'   \item{MassWithAdduct}{Calculated mass including adduct}
#'   \item{MassQueried}{Measured mass that was queried}
#'   \item{DeltaPPM}{Deviation in ppm between measured and expected mass}
#' }
#'
#' @details
#' 1. Calculates the neutral mass for each adduct: `(m - MassDiff)/Nmol`.
#' 2. Compares neutral mass to the reference metabolite monoisotopic mass.  
#' 3. Selects all matches within `ppmtol` ppm.  
#' 4. Returns a combined data frame with reference metabolite info, adduct info, 
#'    calculated mass with adduct, queried mass, and deviation in ppm.
#'
#' @examples
#' \dontrun{
#' # Load reference data from MAIT
#' library(MAIT)
#' data(MAITtables)
#' 
#' # Example: Query m/z = 337.1 for positive adducts
#' refTab <- get_annotation_mass("HMDB")  # reference metabolite table
#' q <- massQuery(
#'   m = 337.1,
#'   ppmtol = 10,
#'   refTab = refTab,
#'   addTable = posAdducts,
#'   ID = "Feature1"
#' )
#' head(q)
#' }
#'
#' @export

massQuery <- function(m, ppmtol, refTab, addTable, ID=NULL) {
  
  if (is.null(refTab$monoMass))
    stop("The 'monoMass' column doesn't exist in the refTab!")
  
  m_mass <- (m - addTable$MassDiff)/addTable$Nmol
  dd <- sapply(m_mass, function(v) {
    v - refTab$monoMass
  })
  
  i <- which(abs(dd) < ppmtol*m/1e6, arr.ind = TRUE)
  if (nrow(i) == 0)
    return(NULL)
  
  at <- data.frame(ID = rep(ID, nrow(i)), stringsAsFactors = FALSE)
  at <- cbind(at, refTab[i[, 1], c("InChIKey","CID","cpdName","formula", "monoMass")], addTable[i[, 2], ])
  at$MassWithAdduct <- at$monoMass + at$MassDiff
  at$MassQueried <- m
  at$DeltaPPM <- abs(apply(i, 1, function(x) dd[x[1], x[2]])/m*1e6)
  at <- cbind(at, refTab[i[, 1], -(1:4)])
  at
}
