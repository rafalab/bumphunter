# matchGenes
# return values are actually nearest *transcripts*
# hg19 only, up/down/r not implemented yet
# NOTE: mc.cores>1 assumes multicore resources such as
# "-pe local N" with N >= mc.cores.  We do not verify.

	#species="human",build="hg18",r=0,up=50000000,down=50000000,
matchGenes <- function(object,promoterDist=2500,verbose=TRUE, all=FALSE,
		genes=NULL, nexons=NULL, EXONS=NULL, job=0, mc.cores=1) {

 if(is.null(genes)) {
  if(verbose) cat("Matching regions to genes.\n")
  X <- nearestgene(object, all=all)	# list of either 2 or 3
  NN <- X$N	# X$d not used here; X$N is a GRangesList
  if (all)
	stopifnot(length(NN) >= nrow(object))
  else
	stopifnot(length(NN) == nrow(object))
  metadata <- elementMetadata(NN)
  Tx <- metadata$Tx		# transcript names
  nexons <- metadata$Nexons
  EXONS <- metadata$Exons	# a CompressedIRangesList
  zero <- integer(length(NN))	# for place-holders
  starts <- object$start
  ends <- object$end
  if (all) {
	queries <- X$q
	times <- rle(queries)$lengths
	starts <- rep(starts, times)
	ends <- rep(ends, times)
  }
  genes <- data.frame(zero, zero, starts, ends, zero,
	as.vector(metadata$Gene), as.vector(metadata$Refseq), zero,
	as.vector(strand(NN)), start(NN), end(NN),
	metadata$CSS, metadata$CSE, stringsAsFactors=FALSE)
  rm(NN, X)
 }

 Ngenes = nrow(genes)
 if(Ngenes > 10000) {
	N = ceiling(Ngenes/10000)
	if(verbose) cat("Splitting the annotation into", N, "chunks of 10000 regions each\n")
	if(verbose && mc.cores>1) cat("mc.cores =", mc.cores, "\n")
	tmp = mclapply(1:N, function(i) {
		start = 10000*(i-1)+1	
		end = min(10000*i, Ngenes)
		if (mc.cores > 1) {
			job = i %% mc.cores
			if (job == 0)
				job = mc.cores
			if(verbose)
				cat(paste("job ", job, ":", sep=""), paste(start, end, sep=":"), "\n")
		}
		else
			job = 0
		II <- start:end
		matchGenes(genes=genes[II,], nexons=nexons[II], EXONS=EXONS[II], job=job)
	}, mc.cores=mc.cores)
	return(do.call(rbind, tmp))
 }

  type=rep("",nrow(genes))
  subtype=rep("",nrow(genes))
  ctype=rep("",nrow(genes))
  dist=rep(0,nrow(genes))	# distance to 5' end of the gene
  insidedistance<-rep(NA,nrow(genes))
  exonnumber<-rep(NA,nrow(genes))

  #nexons<-rep(NA,nrow(genes))
  
  geneL=rep(0,nrow(genes))
  codingL=rep(0,nrow(genes))
  subdist=rep(0,nrow(genes))

  if(verbose) {
	cat("Annotating")
	if (job > 0)
		cat("\n")
  }
  for(i in 1:nrow(genes)){
    if(verbose & i%%1000==0) {
	if (job == 0)
		cat(".")
	else
		cat(paste(job, " ", sep=""))
    }

    #if(!is.na(genes[i,10])){
    if(!is.na(genes[i,12])){
    TS = genes[i,10]
    TE = genes[i,11]
    geneL[i]=TE-TS
    CS = genes[i,12]
    CE = genes[i,13]
    codingL[i]=CE-CS

    #ES = as.numeric(strsplit(genes[i,15],",")[[1]])
    #EE = as.numeric(strsplit(genes[i,16],",")[[1]])
    #Exons= cbind(ES,EE)
    exons <- EXONS[[i]]
    Exons <- cbind(start(exons), end(exons))
    #nexons[i]=nrow(Exons)

    Strand= ifelse(genes[i,9]=="+",1,-1)
    S = genes[i,3]
    E = genes[i,4]

    type[i]=""
    #if(genes[i,3] <= TS & genes[i,4] >= TE){
    if(S <= TS & E >= TE){
      type[i]="covers"
      # dist[i] is already zero
    } else{
      #if(genes[i,3] < TS & genes[i,4] < TS){
      # region totally outside nearest gene
      # region precedes gene
      if(E < TS){
        if(Strand==1){
          type[i]="upstream" 
          dist[i]=TS-E
        } else{
          type[i]="downstream"
          dist[i]=TE-E
        }
      }
      #if(genes[i,3] > TE & genes[i,4] > TE){
      # region follows gene
      if(S > TE){
        if(Strand==-1){
          type[i]="upstream"
          dist[i]=S-TE
        }  else{
          type[i]="downstream"
          dist[i]=S-TS
        }
      }
      #if (genes[i,3]>= TS & genes[i,4] <= TE){
      # totally within gene
      if (S >= TS & E <= TE){
        type[i]="inside"
        if(Strand==-1) dist[i]=TE-E  else dist[i]=S-TS
      }
      
      # overlaps exactly one side of gene ("covers" done above)
      if(type[i]==""){
        #if(genes[i,3] < TS & genes[i,4] <=TE){
        if(S < TS & E <= TE){
          ##OVERLAP FRONT
          if(Strand==1) type[i]="overlaps 5'" else{
	  # If Rafa wants non-zero distance:
          #if(Strand==1) {type[i]="overlaps 5'"; dist[i]=min(TS-S,E-TS);} else{
            type[i]="overlaps 3'"
            dist[i]=TE-E
          }
          S=TS
        }
        #if (genes[i,3] >= TS & genes[i,4] > TE){
        else if (S >= TS & E > TE){
          ##OVERLAP BACK
          if(Strand==-1) type[i]="overlaps 5'" else{
	  # If Rafa wants non-zero distance:
          #if(Strand==-1) {type[i]="overlaps 5'"; dist[i]=min(E-TE,TE-S);} else{
            type[i]="overlaps 3'"
            dist[i]=S-TS
          }
          E=TE
        }
      }
    }

    m1=NA;m2=NA
    if(S >= TS & E <= TE){

      ##INSIDE
      #tmp1=fuzzy.match2(S,Exons)
      #tmp2=fuzzy.match2(E,Exons)
      #tmp1=fuzzy.match2(x=S,z=Exons, x.sorted=TRUE, z.sorted=TRUE)
      #tmp2=fuzzy.match2(x=E,z=Exons, x.sorted=TRUE, z.sorted=TRUE)
      tmp1=fuzzy.match2(S,Exons)
      tmp2=fuzzy.match2(E,Exons)
      
      m1=tmp1[1,1]
      m2=tmp2[1,1]
      exon1=tmp1[1,2]
      exon2=tmp2[1,2]
      m1m2Index=which.min(abs(c(m1,m2)))
      
      if(exon1==exon2 & m1==0 & m2==0){
        subtype[i]="inside exon"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }
      
      if( (sign(m1)==sign(m2) & (m1!=0 & m2!=0) & exon1==exon2) |
         (sign(m1)==-1 & sign(m2)==1 & exon2-exon1==1) ){
        subtype[i]="inside intron"

        insidedistance[i]=c(m1,m2)[m1m2Index]
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
        
      }
      if( (exon2-exon1 > 1) |
         ( (exon2-exon1 == 1) & 
          ((sign(m1)==-1 & sign(m2)==-1) |
           (sign(m1)==1 & sign(m2)==-1) |
           (sign(m1)==1 & sign(m2)==1) |
           (sign(m1)==0 & sign(m2)==-1)|
           (sign(m1)==1 & sign(m2)==0))) |
         (exon2==exon1 & sign(m1)==1 & sign(m2)==-1)){
        subtype[i]="covers exon(s)"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }
      
      if( (exon2-exon1 == 1 & sign(m1)==-1 & sign(m2)==0) |
         (exon2==exon1 & sign(m1)==1 & sign(m2)==0)){
        if(Strand==1) subtype[i]="overlaps exon upstream" else subtype[i]="overlaps exon downstream"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }
      if( (exon2-exon1 == 1 & sign(m1)==0 & sign(m2)==1) |
         (exon2==exon1 & sign(m1)==0 & sign(m2)==-1)){
        if(Strand==-1) subtype[i]="overlaps exon upstream" else subtype[i]="overlaps exon downstream"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }
      if( exon2-exon1 == 1 & sign(m1)==0 & sign(m2)==0){
        subtype[i]="overlaps two exons"
        insidedistance[i]=0
        exonnumber[i]=c(exon1,exon2)[m1m2Index]
      }

      if(Strand!=1){
        insidedistance[i]= -insidedistance[i]
        exonnumber[i] = nrow(Exons) - exonnumber[i] + 1
      }
      
      ctype[i]="inside transcription region"
      
      if(S<CS & E<CS){
        if(Strand==1) ctype[i]="5' UTR" else ctype[i]="3'UTR"
      }
      
      if(S>CE & E>CE){
        if(Strand==-1) ctype[i]="5' UTR" else ctype[i]="3'UTR"
      }
      if(S<CS & E>CE){
        #ctype[i]="covers transcription region"
        ctype[i]="covers coding region"
      }
      if(S<CS & E>CS){
        if(Strand==1) ctype[i]="overlaps 5' UTR" else ctype[i]="overlaps 3'UTR"
      }
      if(S<CE & E>CE){
        if(Strand==-1) ctype[i]="overlaps 5' UTR" else ctype[i]="overlaps 3'UTR"
      }

    }
    if(FALSE){##graphical check
      plot(0,0,ylim=c(0,1),xlim=range(genes[i,c(3:4,10:11)]),xlab=paste("inside distance=",insidedistance[i],m1,m2,exonnumber[i]))
      polygon(c(TS,TE,TE,TS),c(0,0,0.5,0.5),density=0,col=2)
      polygon(c(CS,CE,CE,CS),c(0.1,0.1,0.4,0.4),density=0,col=3)
      abline(h=0.25,lwd=2)
      apply(Exons,1,function(x) polygon(c(x[1],x[2],x[2],x[1]),
                                        c(0.2,0.2,0.3,0.3),col=4))
      polygon(as.vector(genes[i,c(3,4,4,3)]),c(0.75,0.75,0.85,0.85),col=5)
      lines(c(TS,TS+1000),c(0.65,0.65),lwd=3)
      title(paste(i,Strand,type[i],subtype[i],ctype[i],dist[i],sep=":"))
    }
  }
  }
  if(verbose) {
	if (job == 0)
		cat("Done.\n")
	else
		cat(paste("Done(", job, ").\n", sep=""))
  }

  type[dist<=promoterDist & type=="upstream"] <- "promoter"
  type[dist<=promoterDist & type=="downstream"] <- "close to 3'"

  description=type
  tmpIndex=which(description=="inside")
  description[tmpIndex] <- subtype[tmpIndex]
  tmp <- data.frame(name=I(genes[,6]),
                    annotation=I(genes[,7]),
                    description=factor(description,levels=c("upstream","promoter","overlaps 5'","inside intron","inside exon","covers exon(s)","overlaps exon upstream","overlaps exon downstream","overlaps two exons","overlaps 3'","close to 3'","downstream","covers")),
                    region=factor(type,levels=c("upstream","promoter","overlaps 5'","inside","overlaps 3'","close to 3'","downstream","covers")),
                    distance=dist,
#genes=genes,
                    subregion=factor(subtype,levels=c("inside intron","inside exon","covers exon(s)","overlaps exon upstream","overlaps exon downstream","overlaps two exons")),
                    insidedistance=insidedistance,
                    exonnumber=exonnumber,
                    nexons=nexons,
                    UTR=factor(ctype,levels=c("inside transcription region","5' UTR","overlaps 5' UTR","3'UTR","overlaps 3'UTR","covers transcription region")),
                    strand=genes[,9],
                    geneL=geneL,
                    codingL=codingL)
  if (all)
	tmp <- cbind(queries, I(Tx), tmp)
  tmp
}
