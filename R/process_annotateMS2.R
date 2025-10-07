#' Annotate Features Using MS2 Cosine Similarity
#'
#' Computes cosine similarity between measured MS2 spectra of features
#' and reference MS2 spectra provided in \code{ms1Annot}. Only high-purity
#' MS2 spectra (ms2_purity > 0.5) are considered.
#'
#' @param x An object containing feature data (e.g., XCMS \code{featureSet})
#'   with MS2 spectra stored in \code{fData(x@featureSet)$ms2spectrum}.
#' @param ms1Annot A data frame containing reference MS2 spectra and feature IDs.
#'   Must include columns: \code{ID}, \code{ms2_mass}, \code{ms2_intensity},
#'   and \code{ms2_purity}.
#' @param ppmtol Numeric, mass tolerance in ppm for comparing measured and reference spectra.
#'   Default is 20 ppm.
#'
#' @return The input \code{ms1Annot} data frame with an added column:
#'   \code{ms2_cos} containing cosine similarity scores for each feature.
#'
#' @details
#' For each feature with high MS2 purity:
#' \enumerate{
#'   \item Extract measured MS2 spectrum from \code{x@featureSet}.
#'   \item Parse the reference MS2 spectrum from \code{ms1Annot}.
#'   \item Compute the cosine similarity using \code{cospec()}.
#'   \item Store the score in the new column \code{ms2_cos}.
#' }
#' If no MS2 spectrum is present, the score is set to NA.
#'
#' @examples
#' \dontrun{
#' res <- annotateMS2(x = featureSetObj, ms1Annot = referenceDF, ppmtol = 20)
#' head(res$ms2_cos)
#' }
#'
#' @export

annotateMS2 <- function(x, ms1Annot, ppmtol = 20) {
  
  fdata <- fData(x@featureSet)
  im <- which(ms1Annot$ms2_purity > 0.5)
  
  v <- sapply(im, function(i) {
    cos <- NA
    fid <- ms1Annot$ID[i]
    measured_cs <- fdata$ms2spectrum[fdata$ID == fid]
    if (nchar(measured_cs) > 0) {
      measured_cs <- str2spectra(measured_cs)
      if (nchar(ms1Annot$ms2_mass[i]) > 0) {
        ref_cs <- data.frame(
          mz = as.numeric(strsplit(ms1Annot$ms2_mass[i], ";")[[1]]),
          intensity = as.numeric(strsplit(ms1Annot$ms2_intensity[i], ";")[[1]])
        )
        cos <- cospec(measured = measured_cs, standard = ref_cs, ppmtol = ppmtol)
      }
    }
    attr(cos, "ppmtol") <- ppmtol
    cos
  })
  
  ms1Annot$ms2_cos <- NA
  ms1Annot$ms2_cos[im] <- v
  ms1Annot
}

#' Compute Annotation Scores for Features
#'
#' Calculates a composite score for LC–MS features based on MS1 accuracy, 
#' MS2 spectral similarity, and retention time (RT) agreement. The function 
#' integrates individual scores into a final weighted score for annotation ranking.
#'
#' @param x An object containing feature data (e.g., XCMS \code{featureSet}) 
#'   with MS1 and MS2 information. MS2 spectra should be stored in 
#'   \code{fData(x@featureSet)$ms2spectrum}, and retention times in \code{fData(x@featureSet)$rtmed}.
#' @param an2 A data frame containing reference annotation information, including:
#'   \itemize{
#'     \item \code{ID}: Feature IDs corresponding to rows in \code{x}.
#'     \item \code{ms2_mass} and \code{ms2_intensity}: Reference MS2 spectra.
#'     \item \code{ms2_cos}: Cosine similarity scores from measured vs reference MS2.
#'     \item \code{internalStd}: Indicator for internal standard (used in scoring).
#'     \item \code{DeltaPPM}: Mass error for MS1 comparison.
#'     \item \code{RT}: Reference retention times.
#'     \item \code{IPS}: Intensity/presence score weight.
#'     \item \code{primary}: Indicator for primary annotation.
#'   }
#'
#' @return A copy of \code{an2} with additional columns:
#'   \describe{
#'     \item{score_ms1}{MS1-based score for mass accuracy.}
#'     \item{score_ms2}{MS2-based score reflecting spectral similarity.}
#'     \item{score_rt}{Retention time agreement score.}
#'     \item{Score}{Final integrated annotation score combining MS1, MS2, RT, and weighting factors.}
#'   }
#'
#' @details
#' The scoring combines three components:
#' \enumerate{
#'   \item \strong{MS2 score}: Uses \code{ms2Score()} to compare the number 
#'         of measured MS2 peaks against reference peaks and their cosine similarity.
#'   \item \strong{MS1 score}: Uses \code{ms1Score()} to assess mass accuracy 
#'         (deltaPPM) for each feature.
#'   \item \strong{RT score}: Uses \code{rtScore()} to measure retention time 
#'         agreement between observed and reference RT.
#' }
#' The final score is calculated as:
#' \deqn{Score = (score\_ms1 + score\_ms2 + score\_rt) * IPS + primary/2}
#' which integrates all three components weighted by intensity/presence and primary annotation status.
#'
#' @examples
#' \dontrun{
#' # Assuming 'featureSetObj' contains detected features
#' # and 'annotationDF' contains reference annotations
#' scored <- scoreAnnot(featureSetObj, annotationDF)
#' head(scored$Score)
#' }
#'
#' @export
scoreAnnot <- function(x, an2) {
  fdata <- fData(x@featureSet)
  ### scoring annotation
  # MS2 score
  nmea <- stringr::str_count(fdata$ms2spectrum, pattern = ";")/2 + as.integer(nchar(fdata$ms2spectrum) > 0)
  names(nmea) <- rownames(fdata)
  nref <- stringr::str_count(an2$ms2_mass, pattern = ";") + as.integer(nchar(an2$ms2_mass) > 0)
  nref[is.na(an2$ms2_intensity)] <- NA
  scr_ms2 <- ms2Score(n_measured=nmea[an2$ID], n_ref = nref, cos = an2$ms2_cos, std = an2$internalStd)
  
  # MS1 score
  scr_ms1 <- ms1Score(x, deltaPPM = an2$DeltaPPM)
  
  # RT score
  rt <- structure(fdata$rtmed, names = rownames(fdata))
  rtstd <- an2$RT
  scr_rt <- rtScore(rt[an2$ID], an2$RT)
  
  # Final score
  an2$score_ms1 <- scr_ms1
  an2$score_ms2 <- scr_ms2
  an2$score_rt <- scr_rt
  an2$Score <- (scr_ms1 + scr_ms2 + scr_rt) * an2$IPS + an2$primary/2
  an2
}

#' Compute Retention Time (RT) Score
#'
#' Calculates a score for a feature based on the agreement between observed
#' retention time (\code{rt}) and reference retention time (\code{rtStd}).
#' The score is scaled according to deviations and capped by thresholds.
#'
#' @param rt Numeric vector of observed retention times (in seconds or minutes)
#'   for each feature.
#' @param rtStd Numeric vector of reference retention times corresponding to
#'   the same features. Must be the same length as \code{rt}.
#' @param maxScore Numeric scalar, maximum deviation used for scaling. Default
#'   is \code{60*15} (i.e., 15 minutes if RT is in seconds).
#'
#' @return Numeric vector of retention time scores, with the following logic:
#' \describe{
#'   \item{1}{If normalized deviation < 0.05 → score = 1 (excellent match).}
#'   \item{0.75}{If normalized deviation < 0.10 → score = 0.75 (good match).}
#'   \item{-0.25}{If normalized deviation > 0.25 → score = -0.25 (poor match).}
#'   \item{0}{All other cases.}
#' }
#'
#' @details
#' The score is based on the **absolute difference** between observed and
#' reference RT, normalized by \code{maxScore}:
#' \deqn{nd = |rt - rtStd| / maxScore}
#' Thresholds:
#' \itemize{
#'   \item nd < 0.05 → perfect match
#'   \item nd < 0.10 → good match
#'   \item nd > 0.25 → poor match
#' }
#' Values outside these thresholds are scored as 0.
#'
#' @examples
#' rt_obs <- c(300, 305, 400)
#' rt_ref <- c(302, 310, 395)
#' scores <- rtScore(rt_obs, rt_ref)
#' print(scores)
#'
#' @export
rtScore <- function(rt, rtStd, maxScore = 60*15) {
  scr <- rep(0, length(rtStd))
  nd <- abs(rt - rtStd)/maxScore
  scr[which(nd > 0.25)] <- -0.25
  scr[which(nd < 0.1)] <- 0.75
  scr[which(nd < 0.05)] <- 1
  scr
}



#' Compute MS2 Similarity Score for Feature Annotation
#'
#' Calculates an MS2-based annotation score for features by combining cosine 
#' similarity between measured and reference MS2 spectra with additional adjustments 
#' for internal standards and extreme mismatches.
#'
#' @param n_measured Numeric vector of counts of measured MS2 peaks per feature.
#' @param n_ref Numeric vector of counts of reference MS2 peaks per feature.
#' @param cos Numeric vector of cosine similarity scores comparing measured and reference MS2 spectra.
#' @param std Logical vector indicating whether each feature is an internal standard.
#'
#' @return Numeric vector of MS2 scores. Scoring rules:
#' \describe{
#'   \item{Base score}{Set to cosine similarity (\code{cos}), NAs replaced with 0.}
#'   \item{Internal standard, no peaks}{If \code{std = TRUE} and \code{n_measured = n_ref = 0}, score = 0.25.}
#'   \item{Internal standard, small mismatch}{If \code{std = TRUE} and either \code{n_measured >= 2 & n_ref = 0} or \code{n_measured = 0 & n_ref >= 2}, subtract 0.25 from base score.}
#'   \item{Internal standard, large mismatch}{If \code{std = TRUE} and either \code{n_measured >= 6 & n_ref = 0} or \code{n_measured = 0 & n_ref >= 6}, subtract 0.5 from base score.}
#' }
#'
#' @details
#' This function is intended for feature annotation pipelines where MS2 spectral 
#' comparison is used alongside internal standard information. It adjusts the 
#' cosine similarity score to penalize extreme mismatches in internal standards 
#' while giving a small positive score if no peaks are expected.
#'
#' @examples
#' n_measured <- c(5, 0, 7)
#' n_ref <- c(5, 0, 0)
#' cos <- c(0.85, NA, 0.90)
#' std <- c(FALSE, TRUE, TRUE)
#' ms2_scores <- ms2Score(n_measured, n_ref, cos, std)
#' print(ms2_scores)
#'
#' @export
ms2Score <- function(n_measured, n_ref, cos, std) {
  scr_ms2 <- cos
  scr_ms2[is.na(scr_ms2)] <- 0
  
  i <- which(std & n_measured == 0 & n_ref == 0)
  scr_ms2[i] <- 0.25
  i <- which(std & ((n_measured >= 2 & n_ref == 0) | (n_measured == 0 & n_ref >= 2)))
  scr_ms2[i] <- scr_ms2[i] - 0.25
  i <- which(std & ((n_measured >= 6 & n_ref == 0) | (n_measured == 0 & n_ref >= 6)))
  scr_ms2[i] <- scr_ms2[i] - 0.5
  scr_ms2
}







#' Compute MS1 Mass Accuracy Score
#'
#' Calculates a score for LC–MS features based on the deviation of measured
#' m/z from the theoretical/reference mass (deltaPPM). The score is derived
#' from a normal distribution model of the mass error across features.
#'
#' @param x An object containing feature data (e.g., XCMS \code{featureSet}) 
#'   with m/z values stored in \code{x@peak@table$mz} and peak indices in 
#'   \code{fData(x@featureSet)$peakidx}.
#' @param deltaPPM Numeric vector of measured mass deviations in parts per million (ppm) for each feature.
#' @param lower Numeric scalar, minimum standard deviation for scaling. Default is 5 ppm.
#' @param upper Numeric scalar, maximum standard deviation for scaling. Default is 10 ppm.
#'
#' @return Numeric vector of MS1 scores (between 0 and 2) representing mass accuracy:
#' \describe{
#'   \item{Higher scores}{Indicate closer agreement between measured and reference m/z.}
#'   \item{Lower scores}{Indicate larger deviations from the reference m/z.}
#' }
#'
#' @details
#' The function computes the standard deviation of relative m/z values for each feature using 
#' the median absolute deviation (MAD). The SD is then scaled to be within \code{lower} and 
#' \code{upper} bounds. The probability of observing the given \code{deltaPPM} under a 
#' normal distribution is calculated using \code{pnorm}, and the score is set to 
#' \eqn{2 * P(deltaPPM)} to produce a final MS1 accuracy score.
#'
#' @examples
#' # Assuming 'featureSetObj' is an XCMS featureSet object
#' deltaPPM <- c(2.5, 5.1, 8.0)
#' scores <- ms1Score(featureSetObj, deltaPPM)
#' print(scores)
#'
#' @export
ms1Score <- function(x, deltaPPM, lower = 5, upper = 10) {
  fdata <- fData(x@featureSet)
  sds <- sapply(fdata$peakidx, function(i) {
    mz <- x@peak@table$mz[i]
    mad(1e6*mz/median(mz))
  })
  sds[sds==0] <- NA
  sds <- quantile(sds, na.rm = TRUE, probs = 0.75)*2
  sds <- min(upper, sds)
  sds <- max(lower, sds)
  pp <- pnorm(q = deltaPPM, sd = sds, mean = 0, lower.tail = FALSE)
  2*pp
}
