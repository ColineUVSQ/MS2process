#'MS-MS spectrum class to handle an MS-MS spectrum
#'
#' @export
#' @slot precursorMz The mass of the precursor.
#' @slot precursorInt The mass of the precursor.
#' @slot precursorScan The scan of the precursor.
#' @slot energy The energy on which the value have been fragmented.
#' @slot smiles The msiles correspodning to the molecules if available.
#' @slot InChi The InChi corresponding to a molecule.
#' @slot fused A boolean indicating if the molecule have been fused.
#' @slot MSpeaks A data.frame containing the MSMS-spectra in the isoltation windows
#' @slot MSMSpeaks A data.frame containing all the peaks of the msms spectrum.
#' @slot precint A vector containg the intensity of the precursor on the MS spectrum when available.
#' @slot prectime A vector containg the time of the the precursor on the MS spectrum.
#' @slot transitions a data.frame giving the transitions between all the precursor type.
#' @slot polarity the polarity of the acquisition, positive or negative.
#' @slot charge A numeric giving the charge of the molecule.
#' @slot param A list with supplementary parameters stocked.
#' @aliases MSMSspectrum-class
#' @exportClass MSMSspectrum

setClass(
    "MSMSspectrum",slot = list(
        precursorMz="numeric",
        precursorInt="numeric",
        precursorScan="numeric",
        energy="numeric",
        smiles="character",
        InChi="character",
        fused="logical",
        MSpeaks = "data.frame",
        MSMSpeaks = "data.frame",
        precint="numeric",
        prectime="numeric",
        transitions = "data.frame",
        polarity = "numeric",
        charge = "numeric",
        param="list"
    ),
    prototype = list(
        precursorMz=numeric(0),
        precursorInt=numeric(0),
        precursorScan=numeric(0),
        energy=numeric(0),
        smiles=character(0),
        Inchi=chracter(0),
        fused=FALSE,
        MSpeaks=data.frame(),
        MSMSpeaks=data.frame(),
        Precint=numeric(0),
        Prectime=numeric(0),
        transitions = data.frame(),
        polarity=numeric(0),
        charge=numeric(0),
        param=list()
    )
)


###A model for peaks is 
#mz, int, formula, delta_ppm, mass_theo, num (if fused)


setGeneric("getSpec", function(object, ...)
    standardGeneric("getSpec"))
#' Extract a spectra with the selected values filtered for a MSMS spectrum object.
#' 
#' Extract a spectra with the selected values filtered for a MSMS spectrum object.
#' @param object An MSMSspectrum object.
#' @param multiplicity a float giving the fraction of peak to be found or an integer.
#' @param The MS level spectrum to be extracted.
#' @param thresh A float giving the limit in relative intensity.
#' @return The spectra as a matrix.
#' @examples 
#' print("examples to be put here")
#' 
setMethod("getSpec", "MSMSspectrum", function(object,multiplicity= NULL,level=c("MS2","MS1"),thresh=0) {
    level=match.arg(level)
    tab=NULL
    
    ####Checking that it a multiplicity filtering spectra if possible.
    if(multiplicity)
    
    
    if(level=="MS1"){
        tab=object@MSpeaks
    }
    if(level=="MS2"){
        tab=object@MSMSpeaks
    }
    ##If it's a float.
    if(multiplicity<1){
        multiplicity=tab[,"num"]*multiplicity/max(tab[,"num"])
    }
    
    if(thresh!=round(thresh)&thresh<1){
        thresh=max(tab[,"int"])*thresh
        
    }
    pok=which(tab[,"int"]>thresh&tab[,"num"]>=multiplicity)
    if(length(pok)==0){
        toreturn=matrix(0,ncol=2,nrow=0)
        colnames(toreturn)=c("mz","int")
    }
    return(tab[pok,c("mz","int")])
})


#' Get the MS-MS scan corresponding to a mass.
#' 
#' All the MS-MS scans corresponding to a ocmpounds are placed in a list of MS2 object.
#' 
#' @export
#' @param msraw A character giving the path to an object or a mzRramp object.
#' @param mzprecursor The mass to be looked in. 
#' @param tol the tolerance un absolute value.
#' @param MSwindows The windows in which the MS scan need to be returned.
#' @param type Split the list have energy as names. Ordered list is ordered by energy. None as returned by mzR.
#' @param smoothSize the size of the smoothing windows in second. Shloud be a third of the expecte peakwidth approximately.
#' @return A list of MSMSspectrum object
#' @examples
#' print("examples to be put here")
getScanMzr<-function(msraw,mzprecursor,tol=1,MSwindows=10,type=c("split","ordered","none"),smoothSize=10){
    type <- match.arg(type)
    if(class(msraw)!="mzRramp"){
        msraw <- openMSfile(msraw)
    }
    thead <- header(msraw)
    pF <- which((thead[,"precursorMZ"]>=mzprecursor-tol)&(thead[,"precursorMZ"]<=mzprecursor+tol)&
                 thead[,"msLevel"]>=2)
    if(length(pF)==0) return(NA)
    seqMS <- thead[pF,"precursorScanNum"]
    seqTIC <- thead[which(thead[,"msLevel"]==1),"totIonCurrent"]
    seqTime <- thead[which(thead[,"msLevel"]==1),"retentionTime"]
    
    ###Smoothing the intensity to find local maxima.

    
    vMS <- mzR::peaks(msraw,scans=seqMS)
    vMS <- lapply(vMS,centroid.fwhm)
    vScans <- mzR::peaks(msraw,scans=pF)
    vScans <- lapply(vScans,centroid.fwhm)
    vRes <- vector(length=length(vScans),mode="list")
    for(i in 1:length(vScans)){
        spec <- new("MSMSspectrum")
        spec@precursorMz <- thead[pF[i],"precursorMZ"]
        spec@precursorScan <- thead[pF[i],"precursorScanNum"]
        spec@precursorInt <- thead[pF[i],"precursorIntensity"]
        spec@polarity <- thead[pF[i],"polarity"]
        spec@charge <- thead[pF[i],"precursorCharge"]
        spec@energy <- thead[pF[i],"collisionEnergy"]
        vs2 <- mzR::peaks(msraw,scans=thead[pF[i],"precursorScanNum"])
        pMS <- which((vMS[[i]][,"mz"]>=(mzprecursor-MSwindows))&(vMS[[i]][,"mz"]<=(mzprecursor+MSwindows)))
        if(length(pMS)!=0){
        spec@MSpeaks <- vMS[[i]][pMS,c("mz","int")]
        }
        if(nrow(spec@MSpeaks)==0){
            spec@MSpeaks$mz <- numeric(0)
            spec@MSpeaks$int <- numeric(0)
        }
        spec@MSMSpeaks <- vScans[[i]][,c("mz","int")]
        spec@TICtime <- seqTime
        spec@TICint <- seqTIC
        spec@param[["MSwindows"]] <- MSwindows
        spec@time <- thead[pF[i],"retentionTime"]
        vRes[[i]] <- spec
    }
    
    c("splitted","ordered","none")
    if(type=="ordered"){
        vener <- sapply(vRes,function(x){x@energy})
        vener <- order(vener)
        vRes <- vRes[vener]
    }else if(type=="split"){
        ov <- vRes
        vRes <- list()
        vener <- sapply(ov,function(x){x@energy})
        oav <- order(vener)
        tener <- table(vener)
        cumtener <- c(0,cumsum(tener))
        name_ener <- names(tener)
        for(i in 1:length(name_ener)){
            vRes[[name_ener[i]]] <- ov[oav[(cumtener[i]+1):cumtener[i+1]]]
            
        }
    }
    vRes
}



#' plot an MSMSspectrum object.
#'
#' plot an MSMSspectrum object at the MS-MS level or at the MS level.
#' 
#'
#' @export
#' @param x an MSMS spectrum object
#' @param y always null
#' @param level the level of plotting "MS1" or "MS2"
#' @param xlim an optional windows of mass to be plotted.
#' @return nothing
#' @examples
#' print("examples to be put here")

setMethod("plot", "MSMSspectrum", function(x,y=NULL,level=c("MS2","MS1"),xlim=NULL) {
    level <- match.arg(level)
    if(level=="MS2"){
        prec <- which.min(abs(x@MSMSpeaks[,"mz"]-mean(x@precursorMz)))
        if(is.null(xlim)) xlim <- c(0,max(x@MSMSpeaks[,"mz"]))
        title <- NULL
        if(x@fused){
        title <- paste(level,"fused MS-MS spectra","ener",
                    x@energy[1],
                    "precursor",
                    paste(sprintf("%0.4f",range(x@precursorMz)),collapse="-"))
            
        }else{
        title <- paste(level,"spectrum","ener",x@energy,"precursor",sprintf("%0.4f",x@precursorMz))
        }
        plot(x@MSMSpeaks[,"mz"],x@MSMSpeaks[,"int"],type="h",lwd=3,col="black",xlab="m/z",main=title,xlim=xlim,ylab="int")
        points(x@MSMSpeaks[prec,"mz"],x@MSMSpeaks[prec,"int"],lwd=3,col="red",type="h")
    }
    if(level=="MS1"){
        prec <- which.min(abs(x@MSpeaks[,"mz"]-mean(x@precursorMz)))
        if(is.null(xlim)) xlim <- range(x@MSpeaks[,"mz"])
        if(x@fused){
            title <- paste(level,"fused MS spectra","ener",
                        x@energy[1],
                        "precursor",
                        paste(sprintf("%0.4f",range(x@precursorMz)),collapse="-"))
        }else{
            title <- paste(level,"spectrum","ener",x@energy,"precursor",sprintf("%0.4f",x@precursorMz))
        }
        plot(x@MSpeaks[,"mz"],x@MSpeaks[,"int"],type="h",lwd=3,col="black",xlab="m/z",main=title,xlim=xlim,ylab="int")
        points(x@MSpeaks[prec,"mz"],x@MSpeaks[prec,"int"],lwd=3,col="red",type="h")
    }
})

#' plot an MSMSspectrum object.
#'
#' plot an MSMSspectrum object at the MS-MS level or at the MS level.
#' 
#'
#' @export
#' @param x an MSMS spectrum object
#' @param y always null
#' @param level the level of plotting "MS1" or "MS2"
#' @param xlim an optional windows of mass to be plotted.
#' @return nothing
#' @examples
#' print("examples to be put here")
plotPrecursorSeq<-function(lSpec,ylim=NULL,xlim=NULL){
    vtime <- sapply(lSpec,function(x){
        x@time
    })
    vint <- sapply(lSpec,function(x){
        x@precursorInt
    })
    vener <- sapply(lSpec,function(x){
        x@energy
    })
    if(is.null(xlim)) xlim <- range(vtime)
    if(is.null(ylim)) ylim <- c(0,max(vint))
    vu <- sort(unique(vener))
    colvec <- rainbow(length(vu))
    pM <- match(vener,vu)
    title <- paste("precursor intensity sequence")
    plot(vtime,vint,main=title,xlab="Time",ylab="Intensity",col=colvec[pM],type="h",lwd=2)
    legend("topright",as.character(vu),col <- colvec,lwd=2)
    return(invisible(data.frame(time=vtime,int=vint)))
}

findPeaksLimits<-function(seq,pos){
    a <- pos-1
    b <- pos+1
    n <- length(seq)
    #print(seq)
    while(a>1&seq[a]<seq[a+1]){
        a <- a-1
    }
    while(b<n&seq[b]<seq[b-1]){
        b <- b-1
    }
    c(max(a,1),min(b,length(seq)))
}
###Function which check that the precursor intensity peak is on the TIC peak.
havePrecursor<-function(lSpec,TICseq,TICtime,timeseq,type=c("FIA","LC"),SNR=2){
    type <- match.arg(type)
    pmax <- NULL
    seqprec <- plotPrecursorSeq(lSpec)
    seqTime <- sapply(lSpec,function(x){x@time})
    pmmaxPrecursor <- which.max(seqprec[,"int"])
    limit <- NULL
    if(type=="FIA"){
        pm <- which.max(TICseq)
        limit <- findPeaksLimits(TICseq,pm)
        ##half.max
        limit <- c(limit[1]-floor((pm-limit[1])/SNR),limit[2]+floor((limit[2]-pm)/SNR))
        # print(TICtime[limit])
        # print(seqprec[pmmaxPrecursor,"time"])
        if(seqprec[pmmaxPrecursor,"time"]>TICtime[limit[1]]&
           seqprec[pmmaxPrecursor,"time"]<TICtime[limit[2]]){
            return(limit)   
        }else{
            return(NA)
        }
    }
    if(type=="LC"){
        ##Closest pos in time
        pok <- which(seqprec[,"int"]>seqprec[pmmaxPrecursor,"int"]/SNR)
        if(length(pok)==0){
            return(NA)
        }else{
            limit <- range(pok)
            limit[1] <- which.min(abs(seqTime[limit[1]]-1-TICtime))
            limit[2] <- which.min(abs(seqTime[limit[2]]+1-TICtime))
            return(limit)
        }
    }
    ##TODO THE LC-MS case.
}

#'Combines multiples MS-MS spectra as a signle consensus spectrum.
#'
#'Return the consensus spectra of a compounds given a filename and a mass
#'with one spectra by energy.
#'
#'@export
#'@param fname Tthe mzML file to be parsed.
#'@param mz the mass of the ocmpoudns to be found.
#'@param ppm the ppm parameter used by the density function.
#'@param msWindows the size of the ms windows to search the prcursor in.
#'@param multiplicity the multiplicituy to be used for th efiltering of the peaks.
#'@return A list of spectra classified by energy.
#'@examples
#'print("examples to be put here")
combineMSMSSpectra<-function(lspec,ppmMS=3,ppmMSMS=5,SNR=2,smoothSize=10,type="FIA",...){
    ###Checking that the precursor is present
    if(length(lspec)<3){
        warning("impossible to process 3 spectra")
        return(NA)
    }
    sf=floor(smoothSize/mean(diff(lspec[[1]]@TICtime)))
    vSmooth=savgol(lspec[[1]]@TICint,sf+(sf+1)%%2)
    
    
    lim=havePrecursor(lspec,vSmooth,lspec[[1]]@TICtime,type=type,SNR)
    if(is.na(lim[1])){
        warning("no peaks corresponding to injection detected for the given compounds.")
        return(NA)
    }
    vTime=sapply(lspec,function(x){
        x@time
    })
    ToMean=which(vTime>=lspec[[1]]@TICtime[lim[1]]&vTime<=lspec[[1]]@TICtime[lim[2]])
    if(length(ToMean)==0){
        warning("no peaks corresponding to injection detected for the given compounds.")
        return(NA)
    }
    ##CHecking if the multiplicity is a fraction or an integer
    
    
    seqenergy=sapply(lspec[ToMean],function(x){x@energy})
    uener=unique(seqenergy)
    if(length(uener)!=1){
        warning("Multiples energy spectra combined.")
    }
    
    seqprecursor=sapply(lspec[ToMean],function(x){x@precursorMz})
    
    seqtime=sapply(lspec[ToMean],function(x){x@time})

    
    
    ###Getting the msSpecOfTheMeanSpec
    lMS=sapply(lspec[ToMean],function(x){x@MSpeaks},simplify=FALSE)
    lMSMS=sapply(lspec[ToMean],function(x){x@MSMSpeaks},simplify=FALSE)

    #print(seqMSInt)
    #return(seqMSInt)
    ###Obtaining a consensus MS spectra
    cMS=getUniquePeaksl(ppm,lMS,...)
    cMSMS=getUniquePeaksl(ppm,lMSMS,...)

    
    ###Creating the consensus spectra
    spec=new("MSMSspectrum")
    spec@MSpeaks=cMS
    spec@MSMSpeaks=cMSMS
    spec@energy=uener
    spec@precursorMz=seqprecursor
    spec@TICint=lspec[[1]]@TICint
    spec@TICtime=lspec[[1]]@TICtime
    spec@time=seqtime
    spec@polarity=lspec[[1]]@polarity
    spec@charge=lspec[[1]]@charge
    spec@param$ppmMS
    spec@fused=TRUE
    
    #Freturn(listSpec)
    return(spec)
}


aligned<-function(seqm,nmin,joint=TRUE){
    if(joint){
        seqm=sort(seqm)
        vd=length(which(diff(seqm)==1))
        if(vd>=nmin){
            return(TRUE)
        }else{
            return(FALSE)
        }
    }else{
        vd=length(seqm)
        if(vd>=nmin){
            return(TRUE)
        }else{
            return(FALSE)
        }
    }
}


###Grouping masses in a windows. Strongly inspired by group function.
getUniquePeaksl<-function(ppm,listcfia,inter=0.05,nPoints=512,sleep=0,intreat=c("mean","max")){
    #do.call("c",...)
    #if()
    #cat("param",multiplicity,"\n")
    vmz=NULL
    vint=NULL
    if(length(listcfia)>1&typeof(listcfia)=="list"){
        vmz=do.call("c",sapply(listcfia,function(x){x[,"mz"]},simplify=FALSE))
        vint=do.call("c",sapply(listcfia,function(x){x[,"int"]},simplify=FALSE))
        vsample=rep(1:length(listcfia),times=sapply(listcfia,nrow))
    }else{
        vmz=listcfia[[1]][,"mz"]
        vint=listcfia[[1]][,"int"]
        vsample=rep(1,nrow(listcfia[[1]]))
    }
    porder = order(vmz)
    vall=vmz[porder]
    col_val = rainbow(length(listcfia))
    mzrange=range(vall)
    mzrange = c(floor(mzrange[1] / 10) * 10,ceiling(mzrange[2] / 10) * 10)
    massInter = seq(mzrange[1],mzrange[2],inter / 2)
    masspos <- proFIA:::findEqualGreaterM(vall, massInter)
    num_group = 0;
    listgroup = vector(mode = "list",512)
    pos = 0
    previousMz = NULL
    num_group=0
    for (i in seq(length = length(massInter) - 2)) {
        if (massInter[i] %% 100 == 0)
            cat(massInter[i],":",num_group," ")
        mem = i
        start <- masspos[i]
        end <- masspos[i + 2] - 1
        if (end - start <= 0)
            next
        subpeakl = vall
        pos = pos + 1
        bw = ppm * massInter[i] * 1e-6
        den <-
            density(
                subpeakl, bw, from = massInter[i] - 3 * bw, to = massInter[i + 2] +
                    3 * bw,n = nPoints
            )
        
        ###Putting the close to 0 to 0
        maxy = max(den$y)
        den$y[which(den$y <= bw * maxy)] = 0
        plim = c(-1,0,0)
        oldMz = previousMz
        repeat {
            plim = proFIA:::findLimDensity(den$y,plim[2] + 1,plim[3])
            if (plim[1] == plim[2])
                break
            selectedPGroup = which(subpeakl >= den$x[plim[1]] &
                                       subpeakl <= den$x[plim[2]])
            seqscan=vsample[porder[selectedPGroup]]
            #print(seqscan)
            if (!aligned(seqscan,2)) {
                next
            }
            num_group = num_group + 1
            if (num_group > length(listgroup)) {
                listgroup <- c(listgroup, vector("list", length(listgroup)))
            }
            if (mean(vmz[porder[selectedPGroup]]) %in% oldMz ||
                (abs(median(subpeakl[selectedPGroup]) - massInter[i]) < bw *
                 3) ||
                (abs(median(subpeakl[selectedPGroup]) - massInter[i + 2]) <
                 bw * 3)) {
                num_group = num_group - 1
                next
            }
            listgroup[[num_group]] <-c(mean(vmz[porder[selectedPGroup]]),mean(vint[porder[selectedPGroup]]),length(unique(seqscan)))
            previousMz=c(previousMz,mean(vmz[porder[selectedPGroup]]))
        }
    }
    
    resList=vector(mode = "list",num_group)
    if(num_group<1){aa=data.frame(matrix(c(0,0,0),nrow=1,ncol=3))
    colnames(aa)=c("mz","int","num")
    return(aa)
    }
    for(i in 1:num_group){resList[[i]]=listgroup[[i]]}
    to_return=do.call(rbind,resList)
    colnames(to_return)=c("mz","int","num")
    return(as.data.frame(to_return))
}


closeMatch <- function(x,y,ppm=3,dmz=0.005,symmetric=FALSE) {
    
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
        ys<- as.double(ys)
    if (!is.integer(xidx))
        xidx <- as.integer(xidx)
    if (!is.integer(yidx))
        yidx <- as.integer(yidx)
    
    fm <- .Call("closeMatchPpm", xs, ys, xidx, yidx, as.integer(length(x)), as.double(ppm),as.double(0.005),PACKAGE = "MS2process")
    fm2 <- vector("list", length=length(fm))
    ##stop("!")
    if (symmetric){
        for (a in 1:length(fm)) {
            if (!is.null(fm[[a]][1])){
                tmp<-NULL
                for (b in 1:length(fm[[a]])){
                    if ((abs(x[a]-y[fm[[a]]][b]) == min(abs(x[a]-y[fm[[a]][b]]),
                                                        abs(x[a]  -y[fm[[a]][b]-1]),
                                                        abs(x[a]  -y[fm[[a]][b]+1]),
                                                        abs(x[a-1]-y[fm[[a]][b]]),
                                                        abs(x[a+1]-y[fm[[a]][b]]), na.rm=TRUE)
                    )) {
                        tmp<-c(tmp, fm[[a]][b])
                    }
                }
                fm2[[a]]<-tmp
            }
        }
    }else {
        fm2 <- fm}
    fm2
}



#'Extract the MS2 spectra form an FIA acquisition.
#'
#'Return the consensus spectra of a compounds given a filename and a mass
#'with one spectra by energy.
#'
#'@export
#'@param fname Tthe mzML file to be parsed.
#'@param mz the mass of the ocmpoudns to be found.
#'@param ppm the ppm parameter used by the density function.
#'@param msWindows the size of the ms windows to search the prcursor in.
#'@param multiplicity the multiplicituy to be used for th efiltering of the peaks.
#'@return A list of spectra classified by energy.
#'@examples
#'print("examples to be put here")
extractMS2spectra<-function(cfile,ppm=5,tolMS=0.5,mzlist=NULL,
                            rtlist=NULL,type=c("split","ordered","none"),
                            paramcentwave=list(ppm=5,snthresh=3,peakwidth=c(10,50),prefilter=c(3,500))
                            ){
    type=match.arg(type)
    if(require(xcms)){
        if(!all(sapply(cfile,file.exists))) stop("Missing files")
        xraw=xcmsRaw(cfile)
        paramcentwave$object=xraw
        tp=do.call("findPeaks.centWave",args=paramcentwave)
        pok=NULL
        
        if(!is.null(mzlist)){
            lFound=closeMatch(mzlist,tp[,"mz"],ppm=ppm)
            pok=sapply(lFound,is.null)
        }else{
            pok=rep(FALSE,nrow(tp))
            mzlist=c(tp[,"mz"])
        }
        print(sum(pok))
        resList=vector(length(mzlist),mode="list")
        for(k in 1:length(k)){
            if(pok[i]) next
            resList[[k]]=getScanMzr(openMSfile(cfile),mzprecursor = mzlist(pok[i]),MSwindows = tolMS,type = type)
        }
    }
    return(resList)
}

#' Export an MS2 spectra under a format.




####Examples smies
csmils <- '[H][C@@]12CC[C@H](C(C)=O)[C@@]1(C)CC[C@@]1([H])[C@@]2([H])CC[C@@]2([H])CC(=O)CC[C@]12C'
ps <- load.molecules("C:/Users/AD244905/Desktop/test_smiles.txt")

generateFormulaMSMS <- function(spectra,smiles){
	###Parsing the spectra
	
	
}