####Example of mgf file.
#Example of mgf spectra.
# BEGIN IONS
# PEPMASS=413.26611887841
# CHARGE=1+
# 	TITLE=MS/MS scan at 1.535 min with Intensity: 604.0 
# 
# 189.48956 1.9
# 283.62076 3.4
# 301.22977 66.3
# 311.08008 1.3
# 399.99106 2.3
# 
# END IONS

setGeneric("toMgf", function(object, ...)
	standardGeneric("toMgf"))
#' Extract a spectra with the selected values filtered for a MSMS spectrum object.
#' 
#' Extract a spectra with the selected values filtered for a MSMS spectrum object.
#' @param object An MSMSspectrum object.
#' @param multiplicity a float giving the fraction of peak to be found or an integer.
#' @param The MS level spectrum to be extracted.
#' @param thresh A float giving the limit in relative intensity.
#' @param filename The file in which writing the files.
#' @param mzdigit The number of precision in mz.
#' @export
#' @return A chracter vector giving all the possible output.
#' @examples 
#' print("examples to be put here")
setMethod("toMgf", "MSMSspectrum", function(object,filename,level=c("MS2","MS1"),multiplicity=0.5, thresh=0.005, title = NULL, mzdigit = 4,intdigit = 0) {
	###Getting the correct spectra.
	level <- match.arg(level)
	if(!object@fused){
		stop("only fused spectra may be written, use the combineMSMSSpectra() function.")
	}
	reducemf <- function(x){sprintf(paste0('%0.',mzdigit,'f'),x)}
	reduceif <- function(x){sprintf(paste0('%0.',intdigit,'f'),x)}	
	
	if(is.null(title)) title = paste(level,'spectrum for m/z',reducemf(object@precursorMz))
	tPeaks <- getSpec(object,multiplicity= NULL,level=level,thresh=thresh))
	vlines <- character(7+nrow(tPeaks))
	vlines[1] <- 'BEGIN IONS'
	vlines[2] <- paste0('PEPMASS=',reducemf(object@precursorMz))
	vlines[3] <- paste0('CHARGE=',charge,ifelse(object@polarity=='pos','+','-'))
	vlines[4] <- title
	vlines[5:(5+nrow(tPeaks))] <- mapply(reducemf(tPeaks[,1]),reduceif(tPeaks[,2]),FUN=paste)
	vlines[(7+nrow(tPeaks))] <- 'END IONS'
    vlines
)