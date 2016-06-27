annotateNearest <- function(x, subject, annotate=TRUE, ...) {

    if (is.data.frame(x)) x <- makeGRangesFromDataFrame(x)
    
    if(is.data.frame(subject)) subject <- makeGRangesFromDataFrame(subject)
    
    if(class(x)!="GRanges" | class(subject)!="GRanges")
        stop("x and subject must be GRanges or a data.frame with column names chr, start and end.")
    
    dots <- list(...)
    names.dots <- names(dots)
    if ("select" %in% names.dots && dots$select == "all")
	stop("select=\"all\" not supported (yet)")
    
    # find nearests
    NN <- nearest(x, subject, ...)
    nearest_exists <- !is.na(NN)

    # full return value
    Ret <- matrix(NA, nrow=length(NN), ncol=7)
    colnames(Ret) <- c("distance", "subjectHits", "type", "amountOverlap",
	"insideDist", "size1", "size2")
    Ret <- data.frame(Ret, stringsAsFactors=FALSE)

    # queries with a nearest target
    x <- x[nearest_exists]
    # corresponding nearest targets
    y <- subject[NN[nearest_exists],]

    # partial value for when nearest exists
    ret <- matrix(NA, nrow=length(x), ncol=7)
    colnames(ret) <- c("distance", "subjectHits", "type", "amountOverlap",
	"insideDistance", "size1", "size2")
    ret[,"subjectHits"] <- NN[nearest_exists]
    ret <- data.frame(ret, stringsAsFactors=FALSE)

    x.start <- start(x)
    y.start <- start(y)
    x.end <- end(x)
    y.end <- end(y)

    ret$type[(x.start > y.start & x.end <= y.end) |
             (x.start >= y.start & x.end < y.end)] <- "inside"
    ret$type[x.start <= y.start & x.end >= y.end] <- "cover"
    ret$type[x.start > y.end] <- "disjointR"
    ret$type[x.end < y.start] <- "disjointL"
    ret$type[x.start > y.start & x.start <= y.end & x.end > y.end] <- "overlapR"    
    ret$type[x.start < y.start & x.end >= y.start & x.end < y.end] <- "overlapL"
    ret$distance <- 0L

    i <- ret$type == "disjointR"
    ret$distance[i] <- y.end[i] - x.start[i]

    i <- ret$type == "disjointL"
    ret$distance[i] <- y.start[i] - x.end[i]

    if (!annotate) {
	Ret[nearest_exists,] <- ret
	return(as.matrix(Ret[, c("distance", "subjectHits")]))
    }

    i <- ret$type == "overlapR"
    ret$amountOverlap[i] <- -1L * (y.end[i] - x.start[i] + 1L)

    i <- ret$type == "overlapL"
    ret$amountOverlap[i] <- x.end[i] - y.start[i] + 1L

    ret$type[ret$type%in%c("disjointR", "disjointL")] <- "disjoint"
    ret$type[ret$type%in%c("overlapR", "overlapL")] <- "overlap"

    ## insideDistance column:
    i <- ret$type == "inside"	#no missing ret$type at this point
    tmp0 <- cbind(x.end[i] - y.end[i], x.start[i] - y.start[i])
    tmp <- apply(abs(tmp0), 1, which.min)
    tmpinsidedist <- tmp0[, 1]
    tmpIndex <- tmp == 2
    tmpinsidedist[tmpIndex] <- tmp0[tmpIndex, 2]
    ret$insideDistance[i] <- tmpinsidedist

    ## size1 and size2 columns:
    ret$size1 <- x.end - x.start + 1L
    ret$size2 <- y.end - y.start + 1L
    
    Ret[nearest_exists,] <- ret
    Ret
}


#######################################################
#######################################################

annotateTranscripts <-function(txdb, annotationPackage=NULL, by=c("tx","gene"),codingOnly=FALSE,verbose=TRUE,requireAnnotation=FALSE){

    if( class(txdb)!="TxDb") stop("txdb must be of class TxDb")

    if(is.null(annotationPackage)){
        species <- organism(txdb)
        species <- strsplit(species," ")[[1]]
        species <- paste0(substr(species[1],1,1),
                          tolower(substr(species[2],1,1)))
        annotationPackage <-  paste("org",species,"eg.db",sep=".")
        if(verbose)
            message(sprintf("No annotationPackage supplied. Trying %s.",annotationPackage))
    }  

    if(!require(annotationPackage,character.only=TRUE)){
        if(requireAnnotation){
            stop("Can't load ",annotationPackage,".\nMake sure library is installed.\nAnd make sure the species argument follows the convention here http://www.bioconductor.org/packages/release/data/annotation/.\nFor example for human use Hs")}else{
                message("Could not load ",annotationPackage,". Will continue without annotation")
                annotationPackage <- NULL
            }
    }
    
##################################################
### Get TSS and TSE
#######################################################

    if(verbose) message("Getting TSS and TSE.")
    
    by <- match.arg(by)
        
    if(by=="tx"){
        tt <- transcriptsBy(txdb, by="gene")
        seqinfo <- seqinfo(tt)
        ## the transcript names
        LL <- tt@unlistData@elementMetadata@listData
        txx <- unlist(tt)$tx_name
    } else{
        tt <- genes(txdb)
        seqinfo <- seqinfo(tt)
        txx <- names(tt)
    }
    
    
    chr <- unlist(CharacterList(seqnames(tt)))
    strand <- unlist(CharacterList(strand(tt)))
    
    RR <- ranges(tt)
    entrez <- names(RR)
    ## starts / ends (left/right endpoints):
    TSS <- unlist(start(RR))
    TSE <- unlist(end(RR))
    
    
######################################################
### Get CSS and CSE
#######################################################
    
    if(verbose) message("Getting CSS and CSE.")

    cds <- cdsBy(txdb, by=by, use.names=(by=="tx"))
    with_coding <- which(txx %in% names(cds))
    
    ## take only the left- and right-most endpts:
    cds_ranges <- range(cds[txx[with_coding]])
    
    ## starts and ends (really left / right endpts):
    CSS <- CSE <- integer(length(txx))
    CSS[with_coding] <- unlist(start(cds_ranges))
    CSE[with_coding] <- unlist(end(cds_ranges))
    CSS[-with_coding] <- NA
    CSE[-with_coding] <- NA
    
    ##################################################
    ### Get Exons
    #######################################################

    if(verbose) message("Getting exons.")

    ee <- exonsBy(txdb, by=by, use.names=(by=="tx"))
    ## exon groups of interest in transcript order,
    ## and sorted "left to right" within each group
    Exons <- reduce(ranges(ee)[txx])
    Nexons <- elementNROWS(Exons)

    ##now annotate genes
    if(by=="tx") TT <- elementNROWS(tt) else TT <- rep(1,length(tt))
    Entrez <- Rle(entrez, TT)

    Tx <- txx

    if(!is.null(annotationPackage)){
        if(verbose) message("Annotating genes.")

        ##Annotate transcrtipt
        ##cant use select cause map is not 1-1 
        map <- get(gsub("\\.db", "SYMBOL",annotationPackage))
        which <- mappedkeys(map)
        symbols <- sapply(as.list(map[which]), paste, collapse=" ")
        genes <- symbols[entrez]	
        
        map <- get(gsub("\\.db", "REFSEQ",annotationPackage))
        which <- mappedkeys(map)
        symbols <- sapply(as.list(map[which]), paste, collapse=" ")
        refseq <- symbols[entrez]	# a handful of NA
        
        Gene <- Rle(genes, TT)
        Refseq <- Rle(refseq, TT)
    } else {
        Gene <- Rle(NA, sum(TT))
        Refseq <- Rle(NA, sum(TT))
    }

    transcripts=GRanges(ranges=IRanges(start=TSS, end=TSE),
        seqnames=chr, strand=strand,
        seqinfo=seqinfo,
        CSS, CSE, Tx, Entrez, Gene, Refseq, Nexons, Exons)
    if(codingOnly) transcripts <- transcripts[!is.na(values(transcripts)$CSS),]
    attributes(transcripts)$description <- "annotatedTranscripts"
    return(transcripts)
}


matchGenes <- function(x,subject, type=c("any","fiveprime"),
                       promoterDist=2500,
                       skipExons=FALSE,
                       verbose=TRUE){

    if(attributes(subject)$description!="annotatedTranscripts")
        stop("subject must be the output of function annotateTranscripts or have attribute(subject)$description=\"annotatedTranscripts\".")
    
    if (is.data.frame(x)) x <- makeGRangesFromDataFrame(x)
    if(class(x)!="GRanges") stop("x must be GRanges or a data.frame with column names chr, start and end.")

    type=match.arg(type)
    
    if(type=="fiveprime"){
        map <- nearest(x,resize(subject,width=1))
    } else{
        map <- nearest(x,subject)
    }
    
    ind <- which(!is.na(map))
#    dist <- rep(NA,length(map))
#    dist[ind] <- distance(x[ind,],subject[map[ind],])
#    
#    dist[ind] <- dist[ind]*ifelse(
#        (strand(subject[map[ind],])=="+" &
#         end(subject[map[ind],]) < start(x[ind,])) |
#        (strand(subject[map[ind],])=="-" &
#         start(subject[map[ind],]) > end(x[ind,])),-1,1)
    
    type=rep("",length(x))
    subtype=rep("",length(x))
    ctype=rep("",length(x))
    dist=rep(0,length(x))	# distance to 5' end of the gene
    insidedistance<-rep(NA,length(x))
    exonnumber<-rep(NA,length(x))
    nexons <- rep(NA,length(x))
    geneL=rep(0,length(x))
    codingL=rep(0,length(x))
    subdist=rep(0,length(x))
    genenames=rep("",length(x))
    geneannotation=rep("",length(x))
    entrez=rep("",length(x))
    strands=rep("",length(x))

    
    genenames[ind]<-as.character(values(subject[map[ind],])$Gene)
    geneannotation[ind]<-as.character(values(subject[map[ind],])$Refseq)
    entrez[ind]<-as.character(values(subject[map[ind],])$Entrez)
    nexons[ind] <- values(subject[map[ind],])$Nexons
    strands[ind] <- as.character(strand(subject[map[ind],]))
    
    ## This loop could be moved to C. Note: it uses nearest
    for(j in ind){
        
        i <- map[ind][j]
        
        if(verbose & j%%100==0) cat(".")
        
        TS = start(subject)[i]
        TE = end(subject)[i]
        geneL[j] = TE-TS
        if(!is.na(subject$CSS[i])){
            CS = subject$CSS[i]
            CE = subject$CSE[i]
            codingL[j]=CE-CS
        } else {
            CS <- CE <- codingL <- NA
        }
        exons <- subject[i,]$Exons[[1]]
        Exons <- cbind(start(exons), end(exons))
        
        Strand= ifelse(strand(subject[i,])=="+",1,-1)
        S = start(x[j,])
        E = end(x[j,])
        
#        type[j]=""
        
        if(S <= TS & E >= TE){
            type[j]="covers"
            subtype[j]="covers exon(s)"
            
        } else{
            if(E < TS){
                if(Strand==1){
                    type[j]="upstream" 
                    dist[j]=TS-E
                } else{
                    type[j]="downstream"
                    dist[j]=TE-E
                }
            }
            if(S > TE){
                if(Strand==-1){
                    type[j]="upstream"
                    dist[j]=S-TE
                }  else{
                    type[j]="downstream"
                    dist[j]=S-TS
                }
            }
            ## totally within gene
            if (S >= TS & E <= TE){
                type[j]="inside"
                if(Strand==-1) dist[j]=TE-E  else dist[j]=S-TS
            }
            ## overlaps exactly one side of gene ("covers" done above)
            if(type[j]==""){
                if(S < TS & E <= TE){
                    ##OVERLAP FRONT
                    if(Strand==1) type[j]="overlaps 5'" else{
                        type[j]="overlaps 3'"
                        dist[j]=TE-E
                    }
                }
                else if (S >= TS & E > TE){
                    ##OVERLAP BACK
                    if(Strand==-1) type[j]="overlaps 5'" else{
                        type[j]="overlaps 3'"
                        dist[j]=S-TS
                    }
                }
            }
        }
        
        m1=NA;m2=NA
        if( type[j]%in%c("overlaps 5'","overlaps 3'","inside") & !skipExons ){
            
            ir=IRanges(start=c(S,E),width=1)
            map2 = nearest(ir, exons)
            dist2<-distance(ir,exons[map2])
            pos <- start(ir) ##start end the same
            dist2 <- dist2*ifelse(
                (Strand==1 & 
                 end(exons[map2,]) < pos ) |
                (Strand==-1 &
                 start(exons[map2,]) > pos),-1,1)
            
            tmp<-cbind(dist2,map2)
            m1 = tmp[1,1]
            m2 = tmp[2,1]
            exon1 = tmp[1,2]
            exon2 = tmp[2,2]
            m1m2Index=which.min(abs(c(m1,m2)))
            
            if(exon1==exon2 & m1==0 & m2==0){
                subtype[j]="inside exon"
                insidedistance[j]=0
                exonnumber[j]=c(exon1,exon2)[m1m2Index]
            }
            
            if( (sign(m1)==sign(m2) & (m1!=0 & m2!=0) & exon1==exon2) |
               (sign(m1)==-1 & sign(m2)==1 & exon2-exon1==1) ){
                subtype[j]="inside intron"
                
                insidedistance[j]=c(m1,m2)[m1m2Index]
                exonnumber[j]=c(exon1,exon2)[m1m2Index]
                
            }
            if( (exon2-exon1 > 1) |
               ( (exon2-exon1 == 1) & 
                ((sign(m1)==-1 & sign(m2)==-1) |
                 (sign(m1)==1 & sign(m2)==-1) |
                 (sign(m1)==1 & sign(m2)==1) |
                 (sign(m1)==0 & sign(m2)==-1)|
                 (sign(m1)==1 & sign(m2)==0))) |
               (exon2==exon1 & sign(m1)!=sign(m2))){
                subtype[j]="covers exon(s)"
                insidedistance[j]=0
                exonnumber[j]=c(exon1,exon2)[m1m2Index]
            }
            
            if( (exon2-exon1 <= 1 & sign(m1)==-1 & sign(m2)==0) |
               (exon2==exon1 & sign(m1)==1 & sign(m2)==0)){
                if(Strand==1) subtype[j]="overlaps exon upstream" else subtype[j]="overlaps exon downstream"
                insidedistance[j]=0
                exonnumber[j]=c(exon1,exon2)[m1m2Index]
            }
            if( (exon2-exon1 <= 1 & sign(m1)==0 & sign(m2)==1) |
               (exon2==exon1 & sign(m1)==0 & sign(m2)==-1)){
                if(Strand==-1) subtype[j]="overlaps exon upstream" else subtype[j]="overlaps exon downstream"
                insidedistance[j]=0
                exonnumber[j]=c(exon1,exon2)[m1m2Index]
            }
            if( exon2-exon1 == 1 & sign(m1)==0 & sign(m2)==0){
                subtype[j]="overlaps two exons"
                insidedistance[j]=0
                exonnumber[j]=c(exon1,exon2)[m1m2Index]
            }
            
            if(Strand!=1){
                insidedistance[j]= -insidedistance[j]
                exonnumber[j] = nrow(Exons) - exonnumber[j] + 1
            }
            
            ctype[j]="inside transcription region"
            
            if(!is.na(CS)) {
                if(S<CS & E<CS){
                    if(Strand==1) ctype[j]="5' UTR" else ctype[j]="3'UTR"
                }
            }
            
            
            if(!is.na(CE)) {
                if(S>CE & E>CE){
                    if(Strand==-1) ctype[j]="5' UTR" else ctype[j]="3'UTR"
                }
            }
            
            if(!is.na(CS) & !is.na(CE)) {
                if(S<CS & E>CE){
                    ctype[j]="covers coding region"
                }
            }
            
            if(!is.na(CS)) {
                if(S<CS & E>CS){
                    if(Strand==1) ctype[j]="overlaps 5' UTR" else ctype[j]="overlaps 3'UTR"
                }
            }
            
            if(!is.na(CE)) {
                if(S<CE & E>CE){
                    if(Strand==-1) ctype[j]="overlaps 5' UTR" else ctype[j]="overlaps 3'UTR"
                }
            }
            
            
        }
    ##    if(TE-TS<10^5){##graphical check
            
    ##         plot(0,0,ylim=c(0,0.6),xlim=range(c(start(subject[i,]),end(subject[i,]),start(x[j,]),end(x[j,]))),
    ##              xlab=paste("inside distance=",dist[j],insidedistance[j],m1,m2,exonnumber[j]))
    ##         polygon(c(TS,TE,TE,TS),c(0,0,0.5,0.5),density=0,col=2)
    ##         polygon(c(CS,CE,CE,CS),c(0.1,0.1,0.4,0.4),density=0,col=3)
    ##         abline(h=0.25,lwd=2)
    ##         apply(Exons,1,function(x) polygon(c(x[1],x[2],x[2],x[1]),
    ##                                           c(0.2,0.2,0.3,0.3),col=4))
            
    ##         polygon(c(start(x[j,]),end(x[j,]),end(x[j,]),start(x[j,])),c(0.4,0.4,0.5,0.5),col=5)
    ##         lines(c(TS,TS+1000),c(0.55,0.55),lwd=3)
    ##         title(paste(j,i,Strand,type[j],subtype[j],ctype[j],dist[j],sep=":"))
    ##     }

    }
    
    type[dist<=promoterDist & type=="upstream"] <- "promoter"
    type[dist<=promoterDist & type=="downstream"] <- "close to 3'"
    
    description=type
    tmpIndex=which(description=="inside")
    description[tmpIndex] <- subtype[tmpIndex]
    tmp <- data.frame(name=I(genenames),
                      annotation=I(geneannotation),
                      description=factor(description,levels=c("upstream","promoter","overlaps 5'","inside intron","inside exon","covers exon(s)","overlaps exon upstream","overlaps exon downstream","overlaps two exons","overlaps 3'","close to 3'","downstream","covers")),
                      region=factor(type,levels=c("upstream","promoter","overlaps 5'","inside","overlaps 3'","close to 3'","downstream","covers")),
                      distance=dist,
                      subregion=factor(subtype,levels=c("inside intron","inside exon","covers exon(s)","overlaps exon upstream","overlaps exon downstream","overlaps two exons")),
                      insideDistance=insidedistance,
                      exonnumber=exonnumber,
                      nexons=nexons,
                      UTR=factor(ctype,levels=c("inside transcription region","5' UTR","overlaps 5' UTR","3'UTR","overlaps 3'UTR","covers transcription region")),
                      strand=strands,
                      geneL=geneL,
                      codingL=codingL,
                      Entrez=entrez,
                      subjectHits=map)
    return(tmp)
}

                               
