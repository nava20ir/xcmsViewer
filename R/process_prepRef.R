#' Prepare Reference Table for MS Annotation
#' 
#' Prepares a reference metabolite table for MS1/MS2 annotation by selecting the
#' appropriate columns based on ionization mode and setting default values for
#' retention time (RT) and primary compound status.
#'
#' @param ref Data frame of reference metabolites. Must contain columns for 
#'   positive and negative mode MS2 data: `POS_mass`, `POS_intensity`, `POS_purity`, 
#'   `POS_sourceId`, `NEG_mass`, `NEG_intensity`, `NEG_purity`, `NEG_sourceId`.
#' @param mode Character, either `"positive"` or `"negative"`. Determines which set 
#'   of MS2 columns to use.
#' @param primaryInChI Optional character vector of InChIKeys to flag as primary compounds. 
#'   If provided, the `primary` column will be set TRUE for matching compounds.
#'
#' @return A modified reference table with the following columns added or updated:
#' \describe{
#'   \item{ms2_mass}{MS2 m/z values corresponding to the selected ion mode}
#'   \item{ms2_intensity}{MS2 intensities corresponding to the selected ion mode}
#'   \item{ms2_purity}{MS2 purity corresponding to the selected ion mode}
#'   \item{ms2_sourceId}{MS2 source ID corresponding to the selected ion mode}
#'   \item{RT}{Retention time; set to NA if not present in the original table}
#'   \item{primary}{Logical column indicating primary compound status; FALSE by default or TRUE for compounds in `primaryInChI`}
#' }
#'
#' @details
#' - The function standardizes the reference table so that downstream annotation 
#'   functions can work with uniform column names regardless of ion mode.
#' - It ensures that required columns (`RT` and `primary`) exist.
#' - If `primaryInChI` is provided, only compounds whose `InChIKey` matches will 
#'   be flagged as primary.
#'
#' @examples
#' \dontrun{
#' ref <- data.frame(
#'   InChIKey = c("ABC", "DEF"),
#'   POS_mass = c("100;101", "150;151"),
#'   POS_intensity = c("1000;500", "2000;800"),
#'   POS_purity = c(0.9, 0.85),
#'   POS_sourceId = c("S1", "S2"),
#'   NEG_mass = c("99;100", "149;150"),
#'   NEG_intensity = c("900;450", "1800;700"),
#'   NEG_purity = c(0.88, 0.82),
#'   NEG_sourceId = c("S1", "S2")
#' )
#' ref_prepared <- prepRef(ref, mode = "positive", primaryInChI = c("ABC"))
#' }
#'
#' @export
prepRef <- function(ref, mode, primaryInChI = NULL) {
  
  mode <- tolower(mode)
  mode <- match.arg(mode, c("positive", "negative"))
  if (mode == "positive") {
    ref$ms2_mass <- ref$POS_mass
    ref$ms2_intensity <- ref$POS_intensity
    ref$ms2_purity <- ref$POS_purity
    ref$ms2_sourceId <- ref$POS_sourceId
  } else if (mode == "negative") {
    mod <- "neg"
    ref$ms2_mass <- ref$NEG_mass
    ref$ms2_intensity <- ref$NEG_intensity
    ref$ms2_purity <- ref$NEG_purity
    ref$ms2_sourceId <- ref$NEG_sourceId
  } else {
    stop("Unknown mode!")
  }
  
  if (!"RT" %in% colnames(ref))
    ref$RT <- NA
  if (! "primary" %in% colnames(ref))
    ref$primary <- FALSE
  if (!is.null(primaryInChI))
    ref$primary <- ref$InChIKey %in% primaryInChI
  
  ref
}