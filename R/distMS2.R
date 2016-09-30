######################################################################
#Implementation of sp)ectral distance as described by Stein and Scott#
######################################################################


matchPpm <- function(x, y, ppm = 3, mzmin = 0) {
    if (any(is.na(y)))
        stop("NA's are not allowed in y !\n")
    ok <- !(is.na(x))
    ans <- order(x)
    keep <- seq_along(ok)[ok]
    xidx <- ans[ans %in% keep]
    xs <- x[xidx]
    yidx <- order(y)
    ys <- y[yidx]
    if (!is.double(xs))
        xs <- as.double(xs)
    if (!is.double(ys))
        ys <- as.double(ys)
    if (!is.integer(xidx))
        xidx <- as.integer(xidx)
    if (!is.integer(yidx))
        yidx <- as.integer(yidx)
    
    fm <-
        .Call(
            "closeMatchPpm",
            xs,
            ys,
            xidx,
            yidx,
            as.integer(length(x)),
            as.double(ppm),
            as.double(mzmin),
            PACKAGE = "MS2process"
        )
    fm
}


##Stein and scott values : mzexp 3 intexp 0.6
##Massbank values : mzexp 2 intexp 0.5


#' Calculate the composite distance with an adjacency term.
#' 
#' Calculate the composite distance between two MSMSspectrum object.
#' 
#' @export
#' @param s1 z MSMSpsectrum object
#' @param s1 a MSMSpsectrum object
#' @param mzexp The wieghting of the mz term.
#' @param  intexp the weighting of the intensity term.
#' @param ppm the tolerance in the matching in ppm
#' @param mzmin the minimum tolerance
#' @return The composite distance between the 2 spectra.
#' @examples
#' print("examples ot be put here")
compositeDist <- function(s1, s2, mzexp=2, intexp=0.5, ppm, mzmin = 0.005) {
    if(nrow(s1@MSMSpeaks)==0|nrow(s2@MSMSpeaks)==0) return(100)
    matchlist = matchPpm(s1@MSMSpeaks[,"mz"], s2@MSMSpeaks[,"mz"], ppm, mzmin)
    ###Weigthed intensity
    pfound = which(!sapply(matchlist, is.null, simplify = TRUE))
    matchlist=unlist(matchlist[pfound])
    ###If no peak is found.
    if (length(pfound) == 0)
        return(0)
    w1 = (s1@MSMSpeaks[,"int"] ^ intexp) * (s1@MSMSpeaks[,"mz"] ^ mzexp)
    w2 = (s2@MSMSpeaks[,"int"] ^ intexp) * (s2@MSMSpeaks[,"mz"] ^ mzexp)
    Fd = (sum(w1[pfound] * w2[matchlist])^2) / (sum(w1[pfound] ^ 2) *
                                                        sum(w2[matchlist] ^ 2))
    
    ###Adding the factor of the composed distance proposed by stein and scott.
    if (length(pfound) < 2)
        return(Fd)
    n = length(pfound)
    nU = length(s1@MSMSpeaks[,"mz"])
    vFr = ((w1[pfound][-1]) * w2[matchlist][-n]) / ((w1[pfound][-n]) *
                                                                w2[matchlist][-1])
    vSign = sapply(vFr, function(x) {
        if (x > 1) {
            return(-1)
        } else{
            return(1)
        }})
    Fr = sum(vFr ^ vSign) / n
    res = (nU * Fd + n * Fr) / (nU + n)
    res
}


#' Calculate the cosine similarity.
#' 
#' Calculate the cosine similarity btween two MSMSspectrum object.
#' 
#' @export
#' @param s1 z MSMSpsectrum object
#' @param s1 a MSMSpsectrum object
#' @param mzexp The wieghting of the mz term.
#' @param  intexp the weighting of the intensity term.
#' @param ppm the tolerance in the matching in ppm
#' @param mzmin the minimum tolerance
#' @param penality the type of the penality for missing peaks.
#' @return The cosine similarity between the 2 spectra.
#' @examples
#' print("examples ot be put here")
cosineDist <-
    function(s1,
             s2,
             mzexp=2,
             intexp=0.5,
             ppm,
             mzmin = 0.005,penality=c("logint","uniform","none")) {
        penality=match.arg(penality)
        if(nrow(s1@MSMSpeaks)==0|nrow(s2@MSMSpeaks)==0) return(0)
        matchlist = matchPpm(s1@MSMSpeaks[,"mz"], s2@MSMSpeaks[,"mz"], ppm, mzmin)
        ###Weigthed intensity
        pfound = which(!sapply(matchlist, is.null, simplify = TRUE))
        matchlist=unlist(matchlist[pfound])
        ###If no peak is found.
        if (length(pfound) == 0)
            return(0)
        w1 = s1@MSMSpeaks[,"int"] ^ intexp * s1@MSMSpeaks[,"mz"] ^ mzexp
        w2 = s2@MSMSpeaks[,"int"] ^ intexp * s2@MSMSpeaks[,"mz"] ^ mzexp
        
        ratio=1
        if(penality=="uniform"){
            nmissing=nrow(s1@MSMSpeaks)+nrow(s2@MSMSpeaks)-2*length(pfound)
            nfull=2*length(pfound)
            ratio=nfull/(nfull+nmissing)
            
        }else if(penality=="logint"){
            imissing=sum(w1[-pfound])+sum(w2[-matchlist])
            ifound=sum(w1[pfound]+w2[matchlist])
            ratio=ifound/(ifound+imissing)
        }
        
        
        valW = (sum(w1[pfound] * w2[matchlist])^2) / (sum(w1[pfound] ^
                                                                  2) * sum(w2[matchlist] ^ 2))

        valW*ratio
    }

