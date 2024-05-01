############
# Various tools to generate formula for a MSMS spectra.
############


### Interface with the sirius program using R.

####Siruis option and a data.frame
# Usage: [options] ARGUMENTS...
# [--annotate -a] : if set, a csv file is  created additional to the trees. It contains all annotated peaks together with their explanation
# [--auto-charge -Z] : Use this option if the charge of your compounds is unknown and you do not want to assume [M+H]+ as default. With the auto charge option SIRIUS will not care about charges and allow arbitrary adducts for the precursor peak.
# [--cite]
# [--elements -e value] : The allowed elements. Write CHNOPSCl to allow the elements C, H, N, O, P, S and Cl. Add numbers in brackets to restrict the maximal allowed occurence of these elements: CHNOP[5]S[8]Cl[1]
# [--format -O value] : file format of the tree output. Available are 'dot' and 'json'
# [--formula --formulas -f value...] : specify the neutral molecular formula of the measured compound to compute its tree or a list of candidate formulas the method should discriminate. Omit this option if you want to consider all possible molecular formulas
# [--help -h]
# [--ion -i value] : the ionization/adduct of the MS/MS data. Example: [M+H]+, [M-H]-, [M+Cl]-, [M+Na]+, [M]+.
# [--isotope -s value] : how to handle isotope pattern data. Use 'score' to use them for ranking or 'filter' if you just want to remove candidates with bad isotope pattern. Use 'omit' to ignore isotope pattern.
# [--noise value] : median intensity of noise peaks
# [--ms1 -1 value...] : MS1 spectrum file name
# [--ms2 -2 value...] : MS2 spectra file names
# [--no-html] : only for DOT/graphviz output: Do not use html for node labels
# [--no-ion] : only for DOT/graphviz output: Print node labels as neutral formulas instead of ions
# [--candidates -c value] : number of candidates in the output
# [--output -o value] : target directory/filename for the output
# [--ppm-max value] : allowed ppm for decomposing masses
# [--parentmass --precursor --mz -z value] : the mass of the parent ion
# [--profile -p value] : name of the configuration profile. Some of the default profiles are: 'qtof', 'orbitrap', 'fticr'.
# [--version]


###Examples of command line
#sirius -f C20H19NO5 -2 demo/demo-data/txt/chelidonine_msms2.txt demo/demo-data/txt/chelidonine_msms2.txt -O json




SIR.OPT.ANNOTATE <- 'annotate'
SIR.OPT.AUTOCHARGE <- 'auto-charge'
SIR.OPT.ELEMENTS <- 'elements'
SIR.OPT.FORMAT <- 'format'
SIR.OPT.ION <- 'ion'
SIR.OPT.ISOTOPE <- 'isotope'
SIR.OPT.NOISE <- 'noise'
SIR.OPT.CANDIDATES <- 'candidates'
SIR.OPT.PPM_MAX <- 'ppm-max'
SIR.OPT.PARENTMASS <- 'parentmass'
SIR.OPT.PROFILE <- 'profile'
SIR.OPT.FORMULA <- 'formula'
SIR.OPT.MSMS <- '2'
SIR.OPT.MS <- '1'
SIR.OPT.OUTPUT <- 'output'

sirius_opts <- list()
sirius_opts[[SIR.OPT.ANNOTATE]] <- 'character'
sirius_opts[[SIR.OPT.AUTOCHARGE]] <- 'logical'
sirius_opts[SIR.OPT.ELEMENTS[]] <- 'character'
sirius_opts[[SIR.OPT.FORMAT]] <- 'character'
sirius_opts[[SIR.OPT.ION]] <- 'character'
sirius_opts[[SIR.OPT.ISOTOPE]] <- 'character'
sirius_opts[[SIR.OPT.NOISE]] <- 'numeric'
sirius_opts[[SIR.OPT.CANDIDATES]] <- 'numeric'
sirius_opts[[SIR.OPT.PPM_MAX]] <- 'numeric'
sirius_opts[[SIR.OPT.PARENTMASS]] <- 'numeric'
sirius_opts[[SIR.OPT.PROFILE]] <- 'character'
sirius_opts[[SIR.OPT.MSMS]] <- 'character'
sirius_opts[[SIR.OPT.MS]] <- 'character'
sirius_opts[[SIR.OPT.FORMULA]] <- 'character'
sirius_opts[[SIR.OPT.OUTPUT]] <- 'character'

#mz, int, formula, delta_ppm, mass_theo, num (if fused))

###TODO make a function to do this shit.
#sirius_opts_values[[SIR.OPT.ELEMENTS]] <- 'CHNOPSCl'
sirius_opts_values[[SIR.OPT.ISOTOPE]] <- c('omit','filter','score')
sirius_opts_values[[SIR.OPT.PROFILE]] <- c('qtof', 'orbitrap', 'fticr')
sirius_opts_values[[SIR.OPT.FORMAT]] <- c('json', 'dot', 'sirius')


sirius_output_fields <- list()
sirius_output_fields[["mz"]] <- "mz"
sirius_output_fields[["int"]] <- "intensity"
sirius_output_fields[["delta_ppm"]] <- "massdev"
sirius_output_fields[["mass_theo"]] <-  NULL
	
	
makeOptions <- function(list_args){
	
	nArgs <- names(list_args)
	vmatch <- match(nArgs,names(sirius_opts))
	vmatchval <- match(nArgs,names(sirius_opts_values))
	strquery <- 'sirius '
	print(vmatchval)
	for(i in seq_along(list_args)){
		###Wrong arguments names
		if(is.na(vmatch[i])) stop(paste0('unknown argument passed to sirius : ',nArgs[i]))
		
		###Wrong arguments type
		if(class(list_args[[i]])!=sirius_opts[vmatch[i]])stop(paste0('wrong arguments type for',nArgs[i],' : ', class(list_args[[i]]), ' != ', sirius_opts[vmatch[i]]))
		
		###Wrong arguments values
		if(!is.na(vmatchval[i])){
			if(!(list_args[[i]] %in% sirius_opts_values[[vmatchval[i]]])){
				stop(paste0('invalid value for : ',nArgs[i]))
			}
		}
		dash_val <- ifelse(nchar(list_args[[i]]) == 1, '-', '--')
		
		####Adding to the commandline.
		strquery <- c(strquery,paste0(dash_val,nArgs[i]),list_args[[i]])
	}
	
    strquery <- paste(strquery,collapse = " ")
    return(strquery)

}

callSirius <- function(list_args,output_name){
	res <- NULL
	###Making the system call depending of the os.
	if(.Platform$OS.type == "unix") {
		res <- 	system(makeOptions(list_args), intern = TRUE)
	} else {
		res <- 	shell(makeOptions(list_args), intern = TRUE)
	}
    res
}



parseSiriusTree <- function(filename){
	res <- fromJSON(filename)
	
	peaks <- res$fragments[,unlist(sirius_output_fields)]
	colnames(peaks) <- names(sirius_output_fields)
	
	###Parsing the massdev fields
	peaks$delta_ppm<- as.numeric(sapply(peaks$delta_ppm,FUN=function(x){
		strsplit(x,split = ' ')[[1]][1]
	}))
	
	peaks
}

####The method to call the sirius program
setGeneric("callSirius", function(object, ...)
	standardGeneric("callSirius"))
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
setMethod("callSirius", "MSMSspectrum", function(object,level=c("MS2","MS1"),multiplicity=0.5, thresh=0.005, title = NULL, mzdigit = 4,intdigit = 0) {
	

###Testing on an example 
#filename <- 'temp/chelidonine_msms2.json'
#parseSiriusTree(filename)


#Testing of the function
#sirius -f C20H19NO5 -2 demo/demo-data/txt/chelidonine_msms2.txt -O json -o temp




#list_test <- list(
#	formula='C20H19NO5',
#	format='json',
#	output='temp'
#)
#list_test['2'] <- 'demo/demo-data/txt/chelidonine_msms2.txt'

###WORKING EXAMPLES
#callSirius(list_test)

#parseTree <- 
###Testing the dataset