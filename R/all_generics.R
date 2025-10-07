################ dumpPrunedXcmsSet #############
#' @rdname prunedXcmsSet-class
#' @aliases dumpPrunedXcmsSet
#'
#' @title Export a prunedXcmsSet Object to SQLite Database
#'
#' @description
#' This method serializes all relevant components of a \code{\linkS4class{prunedXcmsSet}}
#' object into a SQLite database file. It stores feature tables, expression data,
#' annotations, chromatographic peaks, EIC traces, and associated metadata into
#' separate database tables for downstream visualization or data sharing.
#'
#' @param object A \code{\linkS4class{prunedXcmsSet}} object containing the
#' processed metabolomics data and associated metadata.
#' @param db.file Character string giving the path to the SQLite database file
#' to be created.
#' @param overwrite Logical (default \code{FALSE}). If \code{TRUE}, existing tables
#' in the target database file are replaced. If \code{FALSE}, existing tables
#' are preserved.
#'
#' @details
#' The function exports multiple components of the \code{prunedXcmsSet} object:
#' \itemize{
#' \item Feature-level metadata (\code{xcmsFeatureSet_fData})
#' \item Sample-level metadata (\code{xcmsFeatureSet_pData})
#' \item Feature expression matrix (\code{xcmsFeatureSet_exprs})
#' \item Chromatographic peaks (\code{peak})
#' \item Raw scan metadata and intensities (\code{scanMeta}, \code{scanIntensity})
#' \item Compound annotations (\code{annotation})
#' \item Extracted Ion Chromatograms (\code{EIC})
#' \item Additional visualization parameters and internal attributes
#' }
#'
#' Several indices are automatically created to improve query performance, including
#' unique indices on feature and scan identifiers. Optional normalized expression
#' values (if available) are also exported as a separate table
#' (\code{xcmsFeatureSet_exprs_normalized}).
#'
#' @return
#' No value is returned. The function is called for its side effect of writing
#' tables and indices to the SQLite database specified by \code{db.file}.
#'
#' @section Database Structure:
#' The following main tables are written:
#' \describe{
#' \item{\code{xcmsFeatureSet_fData}}{Feature metadata}
#' \item{\code{xcmsFeatureSet_pData}}{Sample metadata}
#' \item{\code{xcmsFeatureSet_exprs}}{Intensity matrix}
#' \item{\code{annotation}}{Compound annotations}
#' \item{\code{EIC}}{Extracted Ion Chromatograms}
#' }
#'
#' @seealso
#' \code{\linkS4class{prunedXcmsSet}}, \code{\link[DBI]{dbWriteTable}},
#' \code{\link[RSQLite]{SQLite}}
#'
#' @examples
#' \dontrun{
#' library(RSQLite)
#' data(exampleXcmsObject) # hypothetical example
#' dumpPrunedXcmsSet(exampleXcmsObject, db.file = "output.sqlite", overwrite = TRUE)
#' }
#'
#' @export 
setGeneric("dumpPrunedXcmsSet", function(object, ...) {
  standardGeneric("dumpPrunedXcmsSet")
})

#' @rdname prunedXcmsSet-class
setMethod(
  "dumpPrunedXcmsSet", 
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(object, db.file, overwrite = FALSE) {
    
    expr <- exprs(object@featureSet)
    cc <- pData(object@featureSet)$label
    if (is.null(cc)) {
      cc <- pData(object@featureSet)$file
      if (is.null(cc))
        cc <- rownames(cc)
    }
    colnames(expr) <- cc
    expr <- as.data.frame(expr)
    
    fd <- fData(object@featureSet)
    if (is.list(fd$peakidx))
      fd$peakidx <- sapply(fd$peakidx, paste, collapse = ";")
    
    # axis
    ax <- data.frame(
      axis = c("fx", "fy", "sx", "sy"),
      header = c(
        c(attr(object, "fx"), NA)[1],
        c(attr(object, "fy"), NA)[1],
        c(attr(object, "sx"), NA)[1],
        c(attr(object, "sy"), NA)[1]
      ),
      stringsAsFactors = FALSE
    )

    eictab <- object@EIC
    eictab <- lapply(names(eictab), function(x) {
      tab <- eictab[[x]]
      tab$featureID <- x
      tab
    })
    eictab <- do.call(rbind, eictab)

    obj <- list()
    obj$xcmsFeatureSet_fData <- fd
    obj$xcmsFeatureSet_pData <- pData(object@featureSet)
    obj$peak <- object@peak@table
    obj$scanMeta <- object@scan@meta
    obj$scanIntensity <- object@scan@intensity
    obj$annotation <- object@annot
    obj$xcmsFeatureSet_exprs <- expr
    obj$PPMTol <- attr(object@annot, "PPMTol")
    obj$defaultVis <- ax
    obj$misc <- data.frame(
      param = c("keepMS1", "IonMode"),
      value = c(object@scan@keepMS1, attr(object@annot, "IonMode"))
      )
    obj$EIC <- eictab
    # exprs.normlized
    en <- getFeatureExprs(object, normalized = TRUE)
    if (!is.null(en))
      obj$xcmsFeatureSet_exprs_normalized <- data.frame(en)
    
    ## make list as char
    obj$xcmsFeatureSet_fData$"General|Extended|peakidx" <- sapply(      
      obj$xcmsFeatureSet_fData$"General|Extended|peakidx", paste, collapse = ";"
      )
    
    mydb <- DBI::dbConnect(RSQLite::SQLite(), db.file)
    on.exit( DBI::dbDisconnect(mydb) )
    
    for (i in names(obj)) {
      cat(sprintf("Writing table %s ...\n", i))
      if (nrow(obj[[i]]) > 0) {
        DBI::dbWriteTable(mydb, name = i, value = obj[[i]], overwrite = overwrite)
      } else {
        if (overwrite) {
          rs <- RSQLite::dbSendStatement(mydb, sprintf("DROP TABLE IF EXISTS %s;", i))
          RSQLite::dbClearResult(rs)
        }
        ctc <- sprintf("CREATE TABLE %s(%s);", i, paste0("'", paste(colnames(obj[[i]]), collapse="','"), "'"))
        rs <- RSQLite::dbSendStatement(mydb, ctc, overwrite = overwrite)
        RSQLite::dbClearResult(rs)
      }
    }
    
    cat("Creating indices ...\n")
    indices <- c(
      "CREATE UNIQUE INDEX index_fData_ID ON xcmsFeatureSet_fData ('General|All|ID');",
      "CREATE UNIQUE INDEX index_scanMeta_ID ON scanMeta (ID);",
      "CREATE UNIQUE INDEX index_scanMeta_rt_ID_file_msLevel ON scanMeta (rt, ID, fromFile, msLevel);",
      "CREATE INDEX index_scanIntensity_ID ON scanIntensity (ID);",
      "CREATE INDEX index_scanIntensity_mz ON scanIntensity (mz);",
      "CREATE INDEX index_eic_fid ON EIC (featureID);"
    )
    for (stat in indices) {
      res <- RSQLite::dbSendStatement(mydb, stat)
      RSQLite::dbClearResult(res)
    }
  })







######################## annotateMetabolite ########################
#' Annotate Metabolites in a prunedXcmsSet Object
#'
#' @description
#' This method performs both MS1 and MS2 level annotation for features detected
#' in a \code{prunedXcmsSet} object. It integrates reference compound data,
#' computes similarity scores, and categorizes annotation confidence levels.
#'
#' @param object A \code{\link{prunedXcmsSet}} object containing LC–MS feature data.
#' @param mode Character string specifying the ionization mode:
#' either \code{"pos"} (positive) or \code{"neg"} (negative).
#' @param ref A reference compound table containing at least monoisotopic masses
#' and spectral data (e.g., from GNPS, MSDIAL, or in-house libraries).
#' @param ppmtol Numeric. The mass tolerance in parts per million (PPM) used for
#' both MS1 and MS2 annotation.
#' @param ips Numeric between 0 and 1. The minimum ionization probability score
#' used to filter adduct rules from the MAIT package.
#' @param fun_parallel A parallelization function, defaulting to
#' \code{parallel::mclapply}. Used for parallelized feature annotation.
#' @param ... Additional arguments passed to downstream functions.
#'
#' @details
#' The \code{annotateMetabolite()} method performs the following steps:
#' \enumerate{
#' \item Calls \code{\link{annotateMS1}} to match measured m/z values against reference compounds.
#' \item Calls \code{\link{annotateMS2}} to perform MS2 spectral matching.
#' \item Scores annotations using \code{\link{scoreAnnot}}.
#' \item Categorizes annotation levels using \code{\link{categorizeAnnotation}}.
#' \item Updates feature metadata with summarized annotation names (\code{annot_ms1}, \code{annot_ms2}).
#' }
#'
#' The resulting annotations and scores are stored in \code{object@annot}, and
#' corresponding feature-level annotations are written into
#' \code{fData(object@featureSet)}.
#'
#' @return A \code{\link{prunedXcmsSet}} object with added metabolite
#' annotation results.
#'
#' @seealso
#' \code{\link{annotateMS1}}, \code{\link{annotateMS2}}, \code{\link{scoreAnnot}},
#' \code{\link{categorizeAnnotation}}
#'
#' @examples
#' # Assuming 'pxs' is a prunedXcmsSet and 'ref' is a reference table:
#' # annotated_pxs <- annotateMetabolite(pxs, mode = "pos", ref = ref, ppmtol = 20)
#'
#' @export
#' @rdname annotateMetabolite
setGeneric("annotateMetabolite", function(object, ...) {
  standardGeneric("annotateMetabolite")
})

#' @rdname prunedXcmsSet-class
setMethod(
  "annotateMetabolite", 
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(
    object, mode = c('pos', "neg")[1], ref, ppmtol = 25, ips = 0.75, fun_parallel= parallel::mclapply, ...) {
    
    
    x <- object    
    fd <- fData(x@featureSet)    
    fd$annot_ms2 <- NA
    ptol <- data.frame(
      MSLevel = c("MS1", "MS2"),
      value = c(ppmtol, ppmtol),
      stringsAsFactors = FALSE)
    
    an1 <- annotateMS1(
      x, mode = mode, ref = ref, ppmtol = ppmtol, ips = ips, fun_parallel = fun_parallel, ...
    )    
    an2 <- annotateMS2(x = x, ms1Annot = an1, ppmtol = ppmtol)    
    an <- scoreAnnot(x = x, an2 = an2)
    an <- an[order(an$Score, decreasing = TRUE),]    
    an <- an[.xcmsViewerInternalObjects()$xcmsAnnot_column]        
    r <- categorizeAnnotation(fd = fd, an = an, ppmtol = ppmtol)
    an <- r$an
    attr(an, "PPMTol") <- ptol    
    attr(an, "IonMode") <- mode    

    fd <- r$fd
    ms1 <- sapply(unique(an$ID), function(x) {
      nam <- an$cpdName[an$ID == x]
      nam <- nam[1:min(3, length(nam))]
      substr(paste(nam, collapse = ";"), 1, 60)
    })    
    ms1 <- ms1[rownames(fd)]    
    ms1[is.na(ms1)] <- ""    
    fd$annot_ms1 <- ms1    
    
    fData(x@featureSet) <- fd    
    x@annot <- an    
    x
  })




#' Validate Restricted Data Frame Structure
#'
#' @description
#' Internal helper function used to validate that a given data frame has the
#' required column names and classes. Optionally, it can also check whether the
#' data frame contains at least a specified set of columns instead of requiring
#' an exact match.
#'
#' @param df A \code{data.frame} to be validated.
#' @param name Character vector of expected column names.
#' @param class Optional character vector of expected column classes, with the
#' same length and order as \code{name}.
#' @param str Character string used in diagnostic messages to indicate the name
#' or role of the data frame being checked (e.g., "feature table").
#' @param contain Logical (default \code{FALSE}). If \code{TRUE}, the function only
#' checks that all required columns listed in \code{name} are present (additional
#' columns are allowed). If \code{FALSE}, the column names must exactly match
#' those specified in \code{name}.
#'
#' @details
#' The function performs two main checks:
#' \enumerate{
#' \item It verifies that column names match (or contain) the expected names.
#' \item It verifies that each column inherits from the expected class.
#' }
#' If any validation step fails, a diagnostic string is returned describing the
#' issue. Otherwise, the function returns \code{NULL}, indicating success.
#'
#' This function is primarily intended for internal use in S4 object validation
#' (e.g., inside a \code{validObject()} method).
#'
#' @return
#' A character string describing the first problem found, or \code{NULL} if all
#' checks pass.
#'
#' @examples
#' df <- data.frame(a = 1:3, b = letters[1:3])
#' .validRestrictedDataFrame(df, name = c("a", "b"), class = c("integer", "character"), str = "example")
#'
#' # Returns an error message if a column is missing:
#' .validRestrictedDataFrame(df, name = c("a", "b", "c"), class = c("integer", "character", "numeric"), str = "example")
#'
#' @keywords internal
.validRestrictedDataFrame <- function(df, name, class, str, contain = FALSE) {
  if (ncol(df) != 0)  {
    if (contain) {
      i <- setdiff(name, colnames(df))
      if (length(i) > 0) 
        return(sprintf("Missing required columns in %s: %s", str, paste(i, collapse = ", ")))
      df <- df[, name]
    } else {
      if (!identical(colnames(df), name))
        return(sprintf("Problem in %s column name!", str))
    }
    if (!missing(class)) {
      ii <- sapply(seq_along(name), function(x) {
        inherits(df[[x]], c(class[x], "AsIs"))
        })
      # if (!all( sapply(df, class) == class))
      if (!all( ii ))      
        return(sprintf("Problem in %s column class!", str))
    }
  }
}




######################## exportTables ########################
#' Export metabolite and pheno table
#'
#' @param file xlsx file
#' @export
#' @import openxlsx
#' @importFrom Biobase exprs pData fData
#' @rdname prunedXcmsSet-class

setGeneric("exportTables", function(object, file) {
  standardGeneric("exportTables")
})

#' @rdname prunedXcmsSet-class
setMethod(
  "exportTables", 
  signature = signature(object = "prunedXcmsSet"), 
  definition = function(object, file) { 
    pd <- Biobase::pData(object@featureSet)    
    fd <- Biobase::fData(object@featureSet)
    fd$peakidx <- sapply(fd$peakidx, paste, collapse = ",")
    expr <- Biobase::exprs(object@featureSet)
    if ( "label" %in% colnames(pd) ) 
      coln <- pd$label else
        coln <- pd[, 1]
    colnames(expr) <- coln
    .writeMetaboliteTables(pd, fd, expr, file)
  })

#' @rdname prunedXcmsSet-class
setMethod(
  "exportTables", 
  signature = signature(object = "SQLiteConnection"), 
  definition = function(object, file) { 
    pd <- dbGetQuery(object, "SELECT * FROM xcmsFeatureSet_pData")
    fd <- dbGetQuery(object, "SELECT * FROM xcmsFeatureSet_fData")
    expr <- dbGetQuery(object, "SELECT * FROM xcmsFeatureSet_exprs_normalized")
    if ( "label" %in% colnames(pd) ) 
      coln <- pd$label else
        coln <- pd[, 1]
    colnames(expr) <- coln
    .writeMetaboliteTables(pd, fd, expr, file)
  })

.writeMetaboliteTables <- function(pd, fd, expr, file) {
  wb <- createWorkbook("BayBioMS")
  addWorksheet(wb, sheetName = "SampleInfo")
  addWorksheet(wb, sheetName = "MetaboliteInfo")
  addWorksheet(wb, sheetName = "MetaboliteIntensity")
  writeData(wb, sheet = "SampleInfo", x = pd)
  writeData(wb, sheet = "MetaboliteInfo", x = fd)
  writeData(wb, sheet = "MetaboliteIntensity", x = cbind(ID = fd$ID, expr))
  saveWorkbook(wb, file = file, overwrite = TRUE)
}


