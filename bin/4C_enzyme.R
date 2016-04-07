#!/var/www/html/w4cseq/bin/R-3.1.2/bin/Rscript
args <- commandArgs(TRUE)
options(bitmapType='cairo')

path_w4CSeq <- "/var/www/html/"
path_bwa <- "/var/www/html/w4cseq/bin/bwa-0.7.12/"
path_samtools <- "/var/www/html/w4cseq/bin/samtools-1.2/"
path_bedtools <- "/var/www/html/w4cseq/bin/bedtools2-2.25.0/bin/"
path_RCircos <- "/var/www/html/w4cseq/bin/R-3.1.2/library"
path_quantsmooth <- "/var/www/html/w4cseq/bin/R-3.1.2/library"

proc <- args[1]
build <- args[3]
genome <- paste(path_w4CSeq, "/w4cseq/lib/",args[3],"/genome.fa", sep="")
work_dir <- paste(path_w4CSeq, "/w4cseq/work/", args[12], sep="")
file_in <- paste(work_dir,"/",args[2],sep="")
primer_frag <- args[4]
enzyme <- args[5]
enzyme_genome <- paste(path_w4CSeq, "/w4cseq/lib/", args[3], "/enz_sites_R_new/", args[5], "_0.bedsort", sep="")

size_inter <- args[9]
size_intra <- args[10]
window_intra <- args[11]
unzip <- args[13]
chipdata <- args[14]

bait_ch <- args[6]
bait_st <- args[7]
bait_en <- args[8]

setwd(work_dir)
#options(scipen=999)

proc
work_dir
file_in
build
genome
primer_frag
enzyme
enzyme_genome
bait_ch
bait_st
bait_en
size_inter
size_intra
window_intra
unzip
chipdata

FDR<-5

file_sel<-paste(file_in,".sel",sep="")
file_sai<-paste(file_in,".sai",sep="")
file_sam<-paste(file_in,".sam",sep="")
file_bam<-paste(file_in,".bam",sep="")
file_bed<-paste(file_in,".bed",sep="")
file_bedgraph<-paste(file_in,".bedgraph",sep="")
#file_sort_bed<-paste(file_in,".sort.bed",sep="")
file_sort_merge_bed<-paste(file_in,".sort.merge.bed",sep="")
file_sort_merge_filter2_bed<-paste(file_in,".sort.merge.filter2.bed",sep="")
file_sort_merge_filter2_realign_bed<-paste(file_in,".sort.merge.filter2.realign.bed",sep="")
file_sort_merge_filter2_realign_norm_bed<-paste(file_in,".sort.merge.filter2.realign.norm.bed",sep="")
enzyme_no_cut_bed<-paste("enzyme_no_cut.bed",sep="")
file_sort_merge_filter2_realign_all_sort_bed<-paste(file_in,".sort.merge.filter2.realign.all.sort.bed",sep="")
file_sort_merge_filter2_realign_all_sort_count_bed<-paste(file_in,".sort.merge.filter2.realign.all.sort.count.bed",sep="")

system(paste(path_w4CSeq, "/w4cseq/bin/scripts/fastq_select.pl ", file_in, " ", file_sel, " ", primer_frag, " ", enzyme, " ", unzip, sep=""))
system(paste(path_w4CSeq, "/w4cseq/bin/scripts/fastq_convert.pl ", file_sel, " > fastq_convert.fq", sep=""))
system("cp fastq_convert.fq FASTQ_FOR_MAPPING.fq")
system(paste(path_w4CSeq, "/w4cseq/bin/scripts/fastq_filter.pl fastq_convert.fq > FASTQ_FILTERED.fq", sep=""))

system(paste(path_bwa, "/bwa aln -t ", proc, " ", genome, " FASTQ_FILTERED.fq > ", file_sai, sep=""))
system(paste(path_bwa, "/bwa samse ", genome, " ", file_sai, " ", file_sel, " > ", file_sam, sep=""))
system(paste(path_samtools, "/samtools view -bq 1 ", file_sam, " -S > ", file_bam, sep=""))
system(paste(path_samtools, "/samtools sort ", file_bam, " enzyme_sort", sep=""))
system(paste(path_samtools, "/samtools index enzyme_sort.bam", sep=""))
system("cp enzyme_sort.bam MAPPED_BAM.bam")
system("cp enzyme_sort.bam.bai MAPPED_BAM.bam.bai")
system(paste(path_bedtools, "/bamToBed -i ", file_bam, " > ",file_bed, sep=""))
system(paste("cat",file_bed,"| awk '{if($6==\"+\"){print$1\"\t\"$2\"\t\"$2+6\"\t\"$5\"\t\"$6} else {print$1\"\t\"$3-6\"\t\"$3\"\t\"$5\"\t\"$6}}' >",file_bedgraph))
system(paste("sort -k1,1 -k2,2n",file_bedgraph,">", "all_reads.bed"))

system(paste(path_bedtools, "/mergeBed -i all_reads.bed -c 1 -o count -d 0 > ", file_sort_merge_bed, sep=""))
#system(paste("cat ", file_sort_merge_bed, " | awk '$4 > 1' > ", file_sort_merge_filter2_realign_bed, sep=""))
system(paste("cp ", file_sort_merge_bed, " UCSC_view.bed", sep=""))
system(paste("sed -i '1s/^/browser position ", bait_ch, ":", as.numeric(bait_st)-10000, "-", as.numeric(bait_en)+10000, "\\nbrowser hide all\\nbrowser pack refGene encodeRegions\\ntrack type=bedGraph name=\"4C signal (raw reads) (", args[12], ")\" description=\"4C read counts\" db=", build, " visibility=2 color=0,0,0 useScore=1 alwaysZero=on\\n/' UCSC_view.bed", sep=""))

system(paste("cat ", file_sort_merge_bed, " | awk '$4 > 1' > ", file_sort_merge_filter2_bed, sep=""))
system(paste(path_bedtools, "/windowBed -a ", file_sort_merge_filter2_bed, " -b ", enzyme_genome, " -u -w 0 > ",file_sort_merge_filter2_realign_bed, sep=""))
system(paste(path_bedtools, "/intersectBed -a ", enzyme_genome, " -b ", file_sort_merge_filter2_realign_bed, " -v > ", enzyme_no_cut_bed, sep=""))
system(paste("cat", file_sort_merge_filter2_realign_bed, "| awk '{print $1\"\t\"$2\"\t\"$3\"\t1\"}' >", file_sort_merge_filter2_realign_norm_bed))
system(paste("cat", file_sort_merge_filter2_realign_bed, "| awk '$1 != \"chrY\"' > DISTAL_INTERACTION_SITES.bed"))
system(paste("cat", file_sort_merge_filter2_realign_norm_bed, enzyme_no_cut_bed," | sort -k1,1 -k2,2n > ", file_sort_merge_filter2_realign_all_sort_bed))


system(paste(path_w4CSeq, "/w4cseq/bin/scripts/count_sites_binomial.pl ", file_sort_merge_filter2_realign_all_sort_bed, " ", size_inter, " ", size_intra, " ", window_intra, " ", bait_ch, " ", bait_st, " ", bait_en, " ", file_sort_merge_filter2_realign_all_sort_count_bed, sep=""))
system("cp window.bed captured_sites_in_window.bed")
system(paste("sed -i '1s/^/browser position ", bait_ch, ":1-100000000\\nbrowser hide all\\nbrowser pack refGene encodeRegions\\ntrack type=bedGraph name=\"4C signal (binarized) in window (", args[12], ")\" description=\"4C read counts summed in window\" db=", build, " visibility=2 color=255,0,0 useScore=1 alwaysZero=on\\n/' window.bed", sep=""))

system(paste(path_bedtools, "/windowBed -a ", file_sort_merge_filter2_realign_all_sort_count_bed, " -b DISTAL_INTERACTION_SITES.bed -w 0 | awk '{print $1\"\t\"$2\"\t\"$3\"\t\"$4}'> DISTAL_INTERACTION_SITES_pValue.bed", sep=""))

dat_P <- read.table("DISTAL_INTERACTION_SITES_pValue.bed")
dat_P$V5 <- p.adjust(dat_P$V4, "fdr")
dat_P$V4 <- NULL
write.table(dat_P, "DISTAL_INTERACTION_SITES_pValue_adjusted.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

system(paste("awk '$4 <=", FDR/100, "' DISTAL_INTERACTION_SITES_pValue_adjusted.bed | sort -k1,1 -k2,2n > positive_hits.bed", sep=""))
system(paste(path_bedtools, "/windowBed -a ", file_sort_merge_filter2_realign_all_sort_count_bed, " -b positive_hits.bed -w 0 | awk '{print $1\"\t\"$6\"\t\"$7}' > SIGNIFICANT_REGIONS_unmerged.bed", sep=""))
system(paste(path_bedtools, "/mergeBed -i SIGNIFICANT_REGIONS_unmerged.bed | sort -k1,1 -k2,2n > SIGNIFICANT_REGIONS.bed", sep=""))



sig_regions <- read.table("SIGNIFICANT_REGIONS.bed")
sig_regions$V4 <- bait_ch
sig_regions$V5 <- bait_st
sig_regions$V6 <- bait_en
write.table(sig_regions, append=FALSE, quote=FALSE, col.names=FALSE, file="table_for_CIRCOS.txt", row.names = FALSE,sep="\t")


# make a circos plot
library(RCircos, lib.loc=path_RCircos)

if(build=="hg19") {
  data(UCSC.HG19.Human.CytoBandIdeogram)
  cyto.info<-UCSC.HG19.Human.CytoBandIdeogram
}
if(build=="hg18") {
  cyto.info<-read.table(paste(path_w4CSeq, "/w4cseq/lib/hg18/cytobands_hg18.bed", sep=""), header=TRUE)
}
if(build=="mm10") {
  data(UCSC.Mouse.GRCm38.CytoBandIdeogram)
  cyto.info<-UCSC.Mouse.GRCm38.CytoBandIdeogram
}
if(build=="mm9") {
  cyto.info<-read.table(paste(path_w4CSeq, "/w4cseq/lib/mm9/cytobands_mm9.bed", sep=""), header=TRUE)
}

circos<-read.table("table_for_CIRCOS.txt",header=FALSE)
RCircos.Set.Core.Components(cyto.info,chr.exclude=NULL, tracks.inside=5, tracks.outside=0)
pdf(file="circos.pdf",height=8,width=8)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Link.Plot(circos,track.num=2,by.chromosome=TRUE)
dev.off()

png(file="circos.png", width = 8, height = 8, units = 'in', res = 300)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Link.Plot(circos,track.num=2,by.chromosome=TRUE)
dev.off()


#make a genome plot

#########################################################################################
# the source code below supports Mus musculus
#
# quantsmooth.r
# (c) J. Oosting 2006
.quantsmooth<-function(intensities, smooth.lambda=2, tau=0.5, ridge.kappa=0, smooth.na=TRUE) {
  m<-length(intensities)
  nas<-is.na(intensities)
  if (sum(nas)<m) {
    # fill missing values in with valid value
    intensities[nas]<-0
    E<-diag(m)
    Dif<-diff(E)
    # use zero weights for original missing values
    B<-rbind(diag(as.numeric(!nas)),smooth.lambda * Dif)
    ystar = c(intensities, rep(0, m - 1))
    if (ridge.kappa > 0) {
      B <- rbind(B, ridge.kappa * E)
      ystar <- c(ystar, rep(0, m))
    }
    myrq = try(rq.fit(B, ystar, tau=tau, method = "fn"),TRUE)
    if (class(myrq)!="try-error") {
      res<-myrq$coeff
      if (!smooth.na) res[nas]<-NA
      res
    }
    else {
      if (ridge.kappa==0) {
        warning("Problem with fit, repeated with ridge.kappa=(0.001*smooth.lambda)")
        .quantsmooth(intensities,smooth.lambda=smooth.lambda,tau=tau,ridge.kappa=smooth.lambda*0.001,smooth.na=smooth.na)
      }
      else {
        myrq  #Show error
      }
    }
  } else {
    warning("data is all NA, result is NA")
    rep(NA,m)
  }
}

quantsmooth.seg  <-  function(y, x = 1:length(y), lambda = 2, tau = 0.5,
                              kappa = 0, nb = length(x)) {
  # Quantile smoothing with smaller basis
  # Basis has nb segments
  # 
  # Paul Eilers, 2007
  # Based on function .quantsmooth() in package 'quantsmooth' (Oosting et al.)
                            
  # Remove NAs before computation
  nas = is.na(y)
  m = length(y)
  if (m == sum(nas)) {
      warning("data is all NA, result is NA")
      return(NA * y)
  }

  # Construct regression basis
  dx = (1+1e-7) * (max(x) - min(x)) / nb
  ix =  floor((x - min(x)) / dx) + 1
  if (nb == m) ix = 1:m
  B = outer(ix, 1:nb, '==')
  
  # Zero weights for NAs
  B[nas,] = 0
  y[nas] = 0

  # Augment data with penalty stuff
  E = diag(nb)
  D = diff(E)
  B = rbind(B, lambda * D)
  ystar = c(y, rep(0, nb - 1))
  if (kappa > 0) {
    B  =  rbind(B, kappa * E)
    ystar  =  c(ystar, rep(0, nb))
  }
  
  # Try quantile regression
  myrq = try(rq.fit(B, ystar, tau = tau, method = "fn"), TRUE)
  
  # If result is OK, return smooth result
  if (class(myrq) != "try-error") {
      a = myrq$coeff
      z = a[ix]
      return(z)
  }
  
  # If failure, even with kappa >0, return warning
  if (kappa > 0) return(myrq)
  
  # If failure and kappa >0, retry wih reasonable kappa
  warning("Problem with fit, repeated with kappa = 0.001 * lambda")
  z = quantsmooth.seg(y, x, lambda = lambda, tau = tau,
                       kappa = lambda * 0.001, nb = nb)
  return(z)
}


quantsmooth<-function(intensities, smooth.lambda=2, tau=0.5, ridge.kappa=0, smooth.na=TRUE, segment) {
  # if segment is set then the sequence is smoothed with overlapping segments
  # The algorhithm has steeply increasing memory needs for longer sequences
  m<-length(intensities)
  if (missing(segment)) segment<-m
  step.size<-segment %/% 2
  response<-vector(mode="numeric",length=m)
  response[1:min(m,segment)]<-.quantsmooth(intensities[1:min(m,segment)],smooth.lambda,tau,ridge.kappa,smooth.na)
  i.s<-1+step.size
  ol<-segment-step.size
  while ((i.s+step.size) < m) {
    i.e<-min(m,i.s+segment-1)
    tmp.resp<-.quantsmooth(intensities[i.s:i.e],smooth.lambda,tau,ridge.kappa,smooth.na)
    if (ol>0) {
      # set diagonal tapering on overlapping sequencing to prevent abrupt changes on start and end of overlap
      portion<-1:ol / (ol+1)
      response[i.s:(i.s+ol-1)] <- (response[i.s:(i.s+ol-1)]*(1-portion)) + (tmp.resp[1:ol] * portion)
      
      response[(i.s+ol):i.e]<-tmp.resp[(ol+1):length(tmp.resp)]
    } else {
      response[i.s:i.e]<-tmp.resp
    }
    i.s<-i.s+step.size 
  }
  response
}

quantsmooth.cv<-function(intensities, smooth.lambda=2, ridge.kappa=0) {
  m<-length(intensities)
  nas<-is.na(intensities)
  if (sum(nas)<m) {
    # fill missing values in with valid value
    intensities[nas]<-0

    E<-diag(m)
    Dif<-diff(E)
    ystar = as.vector(c(intensities, rep(0, m - 1)))
    weight.odd<-rep(c(1,0),length.out=m)*as.numeric(!nas)
    weight.even<-rep(c(0,1),length.out=m)*as.numeric(!nas)
  
    E.odd<-diag(weight.odd)
    B.odd<-rbind(E.odd,smooth.lambda * Dif)
    E.even<-diag(weight.even)
    B.even<-rbind(E.even,smooth.lambda * Dif)
    if (ridge.kappa > 0) {
      B.odd <- rbind(B.odd, ridge.kappa * E)
      B.even<- rbind(B.even, ridge.kappa * E)
      ystar <- c(ystar, rep(0, m))
    }
  
    myrq.odd = try(rq.fit(B.odd, ystar, method = "fn"),FALSE)
    if (class(myrq.odd)=="try-error") {
      warning("error in fit, result is NA")
      NA
    }
    else {  
      myrq.even = try(rq.fit(B.even, ystar, method = "fn"),FALSE)
      if (class(myrq.even)=="try-error") {
        warning("error in fit, result is NA")
        NA
      }
      else {  
        resid.odd<-intensities-myrq.odd$coefficients
        resid.even<-intensities-myrq.even$coefficients
        #sum of squares van interpolated values
        sum(resid.odd * resid.odd * weight.even,na.rm = TRUE) + sum(resid.even * resid.even * weight.odd,na.rm = TRUE)
      }
    }
  } else {
    warning("data is all NA, result is NA")
    NA
  }  
}
#
getLambdaMin<-function(intensities, lambdas, ...) {
  lambda.res<-rep(NA,length(lambdas))
  for (lambda in 1:length(lambdas)) lambda.res[lambda]<-quantsmooth.cv(intensities,lambdas[lambda],...)
  lambdas[which.min(lambda.res)]
}
#
plotSmoothed<-function(intensities, position, ylim=NULL, ylab="intensity", xlab="position", normalized.to=NULL, grid=NULL, smooth.lambda=2, interval=0.5, plotnew=TRUE, cols, cex.pts=0.6, ...) {
  # plot smoothed data
  # median line is drawn continuous
  # quantile intervals are plotted symmetrical around median ie interval 0.5 plots 0.25 and 0.75 quantiles
  # if intensities contains more than 1 column, the columns are drawn separately
  # position is single vector
  if(is.null(ylim)) ylim<-c(min(intensities,na.rm=TRUE),max(intensities,na.rm=TRUE))
  if(plotnew)plot(c(min(position),max(position)),ylim,ylab=ylab,xlab=xlab,type="n",...)
  if (!is.null(grid)) abline(v=grid,lty=2)
  if (!is.null(normalized.to)) abline(h=normalized.to)
  intensities<-as.matrix(intensities) # make sure it works if only a vector is supplied
	if(missing(cols)) cols<-1:ncol(intensities)+1
  
	idx<-order(position)
  position<-position[idx]
  intensities<-intensities[idx,,drop=FALSE]
  
  for (sample in 1:ncol(intensities)) {
	  if (cex.pts>0) points(position,intensities[,sample],col=cols[sample],pch=20,cex=cex.pts)
    if (sum(!is.na(intensities[,sample]))>10) {
      lines(position, quantsmooth(intensities[,sample],smooth.lambda,segment=150), col=cols[sample], lwd=2)
      if (length(interval)>0) {
        for (i in 1:length(interval)) {
          lines(position, quantsmooth(intensities[,sample],smooth.lambda,tau=0.5-(interval[i]/2),segment=150), col=cols[sample], lty=1+i)
          lines(position, quantsmooth(intensities[,sample],smooth.lambda,tau=0.5+(interval[i]/2),segment=150), col=cols[sample], lty=1+i)
        }
      }
    }
  }
}

getChangedIdx<-function(changed, up) {
  if (sum(changed,na.rm=TRUE)>0) {
    changed[is.na(changed)]<-FALSE
    position<-which(xor(c(FALSE,changed),c(changed,FALSE)))
    startidx<-seq(1,by=2,length.out=length(position) / 2) # odd indexes
    startpos<-position[startidx]
    endpos<-position[startidx+1]-1
    data.frame(up=up,start=startpos,end=endpos)
  } else NULL
}

getChangedRegions<-function(intensities, positions, normalized.to=1, interval, threshold, minlength=2, ...) {
  # determine regions with changes after smoothing
  # normalized.to: value to compare with
  # smooth.lambda: smoothing parameter
  # interval     : changes are defined by these smoothed boundaries crossing normalized.to
  # treshold     : changes are defined by croosing of signal outside of normalized.to + or - treshold 
  #                (only one of treshold or interval can be defined)
  # minlength    :  minimum length of a change to be listed
  #
  # value        : dataframe 3 columns up, start, end
	if (missing(positions)) positions<-1:length(intensities)
	if (!is.null(match.call()$tau)) stop("tau is set by the function")
	if (length(positions)!=length(intensities)) stop("Length of positions argument should be equal to length of intensities argument")
  if (!missing(interval)) {
    res<-rbind(getChangedIdx(quantsmooth(intensities,tau=0.5-(interval/2),...) > normalized.to,TRUE),
          getChangedIdx(quantsmooth(intensities,tau=0.5+(interval/2),...) < normalized.to,FALSE))
  } else if (!missing(threshold)) {
    smoothed<-quantsmooth(intensities,tau=0.5,...)
    res<-rbind(getChangedIdx(smoothed > (normalized.to+threshold),TRUE),
          getChangedIdx(smoothed < (normalized.to-threshold),FALSE))
  } else stop("Either treshold or interval should be defined")
  if (!is.null(res)) {
	  res[,"start"]<-positions[res[,"start"]]
	  res[,"end"]<-positions[res[,"end"]]
	} 
	res 
}

# Support functions for SnpSetIllumina
numericCHR<- function(CHR, prefix="chr") {
  # Set autosomal chromosomes to their number
  # X - 98
  # Y - 99
  # XY - 100
  # MT, M - 101
  # Variable regions - 102
  CHR<-sub(paste0("^",prefix),"",as.character(CHR))
  CHR[grep("_",CHR)]<-"102"
  CHR[CHR=="X"]<-"98"
  CHR[CHR=="Y"]<-"99"
  CHR[CHR=="XY"]<-"100"
  CHR[CHR %in% c("M","MT")]<-"101"
  as.numeric(CHR)
}
#
characterCHR<- function(CHR,prefix="") {
  CHR<-as.character(CHR)
  CHR[CHR=="98"]<-"X"
  CHR[CHR=="99"]<-"Y"
  CHR[CHR=="100"]<-"XY"
  CHR[CHR=="101"]<-"MT"
  CHR[CHR=="102"]<-"-" 
  paste(prefix,CHR,sep="")
}
#
scaleto <-function(x,fromlimits=c(0,50),tolimits=c(0.5,-0.5),adjust=TRUE) {
  if (adjust) {
    x[x>fromlimits[2]]<-fromlimits[2]
    x[x<fromlimits[1]]<-fromlimits[1]
  }  
  x<- x-fromlimits[1]
  x<- x/(fromlimits[2]-fromlimits[1])
  x<- x * (tolimits[2]-tolimits[1])
  x+tolimits[1]
}  
#
.getChrombands<-function(units) {
  if (units %in% c("cM","bases","ISCN")) {
    chrom.bands<-NULL;rm(chrom.bands) # trick to satisfy R check
    data(chrom.bands,package="quantsmooth",lib.loc=path_quantsmooth,envir=environment())
    bandpos<-switch(units,
                    cM =chrom.bands[,c("cM.top","cM.bot")],
                    bases = chrom.bands[,c("bases.top","bases.bot")],
                    ISCN =  chrom.bands[,c("ISCN.top","ISCN.bot")])
    data.frame(chr=chrom.bands$chr,segstart=bandpos[,1],segend=bandpos[,2],stain=chrom.bands$stain,band=chrom.bands$band,arm=chrom.bands$arm, stringsAsFactors=FALSE)
  } else {
    data(list=paste0("chrom.bands.",units),package="quantsmooth",lib.loc=path_quantsmooth, envir=environment())
    .convertUCSCcytoband(get(paste0("chrom.bands.",units)))
  }
}
#
.convertUCSCcytoband<- function(ucscdata) {
  data.frame(chr=characterCHR(numericCHR(ucscdata[,1])),segstart=ucscdata[,2],segend=ucscdata[,3],stain=ucscdata[,5],band=substring(ucscdata[,4],2),arm=substring(ucscdata[,4],1,1), stringsAsFactors=FALSE)
}
#
plotChromosome<-function(gendata,chrompos,chromosome,dataselection=NULL,ylim=NULL,normalized.to=NULL,grid=NULL,smooth.lambda=2,interval=0.5,...) {
  # uses gcsmoothing.R
  if (is.null(dataselection)) dataselection<-rep(TRUE,ncol(gendata))
  plotSmoothed(gendata[chrompos[,"CHR"]==chromosome,dataselection],chrompos[chrompos[,"CHR"]==chromosome,"MapInfo"],ylim=ylim,normalized.to=normalized.to,
               grid=grid,smooth.lambda=smooth.lambda,interval=interval,...)
}
#
prepareGenomePlot<-function(chrompos=NULL,cols="grey50",paintCytobands=FALSE,bleach=0,topspace=1,organism,sexChromosomes=FALSE,units="hg19",...) {
  # prepare plot with lines and axes indicating all chromosomes
  # sends extra arguments to plot function
  cytobandWidth<-0.2
  # hsa 22+ XY
  # mmu 19 + XY
  # rno 20 + XY
	par(mar=c(1,4,2,3)+0.1)

	if (!missing(organism)) {
	  organism<-match.arg(organism,c("hsa","mmu","rno"))
	  chrom.n<-switch(organism,
	                   hsa = 22,
                     mmu = 19,
                     rno = 20) 
    if (is.null(chrompos)) {
      chroms=c(1:chrom.n,if(sexChromosomes)c(98,99)else NULL)
      chrnames=characterCHR(chroms,"chr")
      mapinfo=rep(0,length(chroms))
      chrompos=cbind(CHR=chroms,MapInfo=mapinfo)
      rownames(chrompos)<-chrnames
    }
  	chrs2<-factor(numericCHR(chrompos[,"CHR"]),levels=c(1:chrom.n,if(sexChromosomes)c(98,99)else NULL))
  	if (organism %in% c("hsa","mmu"))
      lens<-lengthChromosome(levels(chrs2),units=units)
  	else
    	lens<-sapply(split(chrompos[,"MapInfo"],chrs2),function(x)max(c(0,x)))
  	names(lens)<-characterCHR(names(lens))
  	cols<-rep(cols,length.out=length(lens))
  	names(cols)<-names(lens)
  	dwidth<-NULL
  	# plot 2 columns of chromosomes, first column large->small (1-12), second column small->large (22-12)
  	for (i in 1:(chrom.n %/% 2)) dwidth[i]<-lens[i]+lens[chrom.n+1-i]
    # make sure vector length equals nr of rows in plot
  	if (chrom.n %% 2 ==1) dwidth<-c(dwidth,lens[chrom.n %/% 2 +1])
  	if (sexChromosomes) dwidth<-c(dwidth,lens["X"]+lens["Y"])
  	maxdwidth<-max(dwidth)*1.05
  	leftrow<-c(if(sexChromosomes)"X" else NULL,((chrom.n + 1) %/% 2):1)
  	rightrow<-c(if(sexChromosomes)"Y" else NULL, if (chrom.n %% 2 ==1) "" else NULL,((chrom.n + 1) %/% 2 +1):chrom.n)
  	plot(c(0,maxdwidth),c(0.5 ,0.5+length(dwidth)+topspace),type="n",ylab="Chromosome",xlab="",axes = FALSE, las = 2,...)
  	axis(2, c(1:length(dwidth)), characterCHR(leftrow), las = 2)
  	axis(4, c(1:length(dwidth)), characterCHR(rightrow), las = 2)
  	if (paintCytobands && organism %in% c("hsa","mmu")) {
    	for (i in 1:length(dwidth)) {
    	  if (lens[leftrow[i]]>0) paintCytobands(leftrow[i],c(0,i+cytobandWidth/2),units=units,width=cytobandWidth,length.out=lens[leftrow[i]],legend=FALSE,bleach=bleach)
    	  if (rightrow[i]!="" && lens[rightrow[i]]>0) paintCytobands(rightrow[i],c(maxdwidth-lens[rightrow[i]],i+cytobandWidth/2),units=units,width=cytobandWidth,length.out=lens[rightrow[i]],legend=FALSE,bleach=bleach)
  	  }
  	} else {
    	for (i in 1:length(dwidth)) {
    	  lines(c(0,lens[leftrow[i]]),c(i,i),col=cols[leftrow[i]],lwd=2)
    	  if(rightrow[i]!="") lines(c(maxdwidth-lens[rightrow[i]],maxdwidth),c(i,i),col=cols[rightrow[i]],lwd=2)
    	}
    }
    # for each locus determine postion on plot , this can be used later to fill with data
  	dchrompos<-matrix(0,nrow=nrow(chrompos),ncol=2,dimnames=list(rownames(chrompos),c("CHR","MapInfo")))
 		for (i in 1:length(rightrow)) if (rightrow[i]!="") {
 		  probes<-characterCHR(chrompos[,"CHR"])==rightrow[i]
      dchrompos[probes,2]<-chrompos[probes,"MapInfo"]+maxdwidth-lens[rightrow[i]]
      dchrompos[probes,1]<- i
    }
  	for (i in 1:length(leftrow)) {
 		  probes<-characterCHR(chrompos[,"CHR"])==leftrow[i]
			dchrompos[probes,2]<-chrompos[probes,"MapInfo"]
      dchrompos[probes,1]<- i
    }
	}
  else {
  	chrs2<-factor(numericCHR(chrompos[,"CHR"]))
  	lens<-sapply(split(chrompos[,"MapInfo"],chrs2),max)
  	m<-length(lens)
  	cols<-rep(cols,length.out=m)
    maxdwidth<-max(lens)
  	plot(c(0,maxdwidth),c(0.5,m+0.5+topspace),type="n",ylab="Chromosome",xlab="",axes = FALSE, las = 2,...)
  	axis(2, c(m:1), characterCHR(names(lens)), las = 2)
    for (i in 1:m)  lines(c(0,lens[i]),c(m+1-i,m+1-i),col=cols[as.numeric(names(lens))],lwd=2)
    dchrompos<-chrompos
    dchrompos[,1]<-m+1-as.numeric(chrs2)
  }
	dchrompos
}
#  data taken from lodplot package
#  original data available at: ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/mapview/BUILD.35.1/ideogram.gz.
qs.semicircle <- function(base.x, base.y, base.length, height=base.length, side=1, orientation=NULL,plottype="poly",...) {
  # based on lodplot package
  # David Duffy <David.Duffy@qimr.edu.au>
  # URL: http://www.qimr.edu.au/davidD
  # - col is now propagated through ..., other plotting parameters can now also be given
  # - different types poly/line 
  radius<-base.length/2
  x<-radius*seq(-1,1,length=40)
  y<-height/radius*sqrt(radius^2-x^2)
  if (is.null(orientation)) {
    co<-as.integer(cos(pi*(3-side)/2))
    so<-as.integer(sin(pi*(3-side)/2))
  }else{
    co<-cos(orientation)
    so<-sin(orientation)
  }
  tx<-co*x - so*y
  ty<-so*x + co*y
  if (is.null(orientation)) {
    if (side==1 || side==3) {
      base.x<-base.x+radius
    }else if (side==2 || side==4) {
      base.y<-base.y+radius
    }
  }
  x<-base.x+tx
  y<-base.y+ty
  switch(plottype,
    poly=polygon(x,y,...),
    line=lines(x,y,...)
  )
}

qs.grid.semicircle<-function (base.x, base.y, base.length, height = base.length, 
    side = 1, orientation = NULL, ...) 
{
    radius <- base.length/2
    x <- radius * seq(-1, 1, length = 40)
    y <- height/radius * sqrt(radius^2 - x^2)
    if (is.null(orientation)) {
        co <- as.integer(cos(pi * (3 - side)/2))
        so <- as.integer(sin(pi * (3 - side)/2))
    }
    else {
        co <- cos(orientation)
        so <- sin(orientation)
    }
    tx <- co * x - so * y
    ty <- so * x + co * y
    if (is.null(orientation)) {
        if (side == 1 || side == 3) {
            base.x <- base.x + radius
        }
        else if (side == 2 || side == 4) {
            base.y <- base.y + radius
        }
    }
    x <- base.x + tx
    y <- base.y + ty
    grid.polygon(x, y, ...)
}

#
paintCytobands<-function(chrom, pos=c(0,0), units="hg19", width=0.4, length.out, bands="major", orientation=c("h","v"), legend = TRUE, cex.leg=0.7, bleach = 0,...) {
  # Based on paint.chromosome from lodplot package
  # added:
  #  -bleach
  #  -length.out
  #  -using all of cM,bases,ISCN
  #  -using hatches for stalk, acen
  #  -legend + cex.leg
  #  -orientation
  #  extracted semicircle for general use
  bleacher<-function(x) { (x * (1-bleach)) + bleach}
  if (class(units)=="data.frame") {
    chrombands<-.convertUCSCcytoband(units)
  } else{
    chrombands<-.getChrombands(units)
  }
  chrom<-switch(as.character(chrom),
         "98"="X",
         "99"="Y",
         as.character(chrom))
  orientation<-match.arg(orientation)
  # original function only required ypos
  if (length(pos)==1) pos<-c(0,pos)
  chromdata<-subset(chrombands, chrombands$chr==chrom)
  if (nrow(chromdata)>0){
    lc<-nchar(chromdata$band)
    sel<-!(substr(chromdata$band,lc,lc) %in% letters)
    if (bands!="major") sel<-!sel
    chromdata<-chromdata[sel,]
    rm(lc,sel)

    type.b<-match(chromdata$stain,c("acen","gneg", "gpos", "gvar", "stalk","gpos25","gpos50","gpos75","gpos100","gpos33", "gpos66"))
    bandpos<-chromdata[,c("segstart","segend")]
    bandcol<-gray(bleacher(c(0.5,1,0.2,0.6,0.75,0.7,0.5,0.3,0.1,0.6,0.4)))[type.b]
    banddens<-c(30,-1,-1,-1,10,-1,-1,-1,-1,-1,-1)[type.b]
    bandbord<-gray(bleacher(c(0,0,0,0,1,0,0,0,0,0,0)))[type.b]
    if (!missing(length.out)) {
      bandpos<-(bandpos/max(bandpos))*length.out
    }
    n<-nrow(chromdata)
    centromere<-which(chromdata$arm[-n]!=chromdata$arm[-1])
    if (length(centromere==1)) {
      idx<-c(2:(centromere-1), (centromere+2):(n-1))
    } else {
      idx<-c(2:(n-1))
    }
    if (orientation=="h") {
      rect(pos[1]+bandpos[idx,1],pos[2],pos[1]+bandpos[idx,2],pos[2]-width, col=bandcol[idx], density=banddens[idx], border=bandbord[idx])
      qs.semicircle(pos[1]+bandpos[1,2], pos[2]-width, width,
                 bandpos[1,2]-bandpos[1,1], 2, col=bandcol[1], density=banddens[1], border=bandbord[1],...)
      qs.semicircle(pos[1]+bandpos[n,1], pos[2]-width, width,
                 bandpos[n,2]-bandpos[n,1], 4, col=bandcol[n], density=banddens[n], border=bandbord[n],...)
      if (length(centromere==1)) {
        qs.semicircle(pos[1]+bandpos[centromere,1], pos[2]-width, width,
                   bandpos[centromere,2]-bandpos[centromere,1],
                   4, col=bandcol[centromere], density=banddens[centromere], border=bandbord[centromere],...)
        qs.semicircle(pos[1]+bandpos[centromere+1,2], pos[2]-width, width,
                   bandpos[centromere+1,2]-bandpos[centromere+1,1],
                   2, col=bandcol[centromere+1], density=banddens[centromere+1], border=bandbord[centromere+1],...)
  
        centromere.size=0.6*0.5*width/yinch(1)
      symbols(pos[1]+bandpos[centromere,2], pos[2]-0.5*width,circles=1,inches=centromere.size, add=TRUE,fg=gray(bleacher(0)),bg="white",...)
      }
      if (legend) text(pos[1]+(bandpos[,1]+bandpos[,2])/2,pos[2]+0.5*width,paste(chromdata[,"arm"],chromdata[,"band"],sep=""),adj=c(0,0.5),srt=90,cex=cex.leg,...)
    } else {
      rect(pos[1],pos[2]-bandpos[idx,1],pos[1]-width,pos[2]-bandpos[idx,2], col=bandcol[idx], density=banddens[idx], border=bandbord[idx],...)
      qs.semicircle(pos[1]-width, pos[2]-bandpos[1,2], width,
                 bandpos[1,2]-bandpos[1,1], 3, col=bandcol[1], density=banddens[1], border=bandbord[1],...)
      qs.semicircle(pos[1]-width, pos[2]-bandpos[n,1], width,
                 bandpos[n,2]-bandpos[n,1], 1, col=bandcol[n], density=banddens[n], border=bandbord[n],...)
      if (length(centromere==1)) {
        qs.semicircle(pos[1]-width, pos[2]-bandpos[centromere,1], width,
                   bandpos[centromere,2]-bandpos[centromere,1],
                   1, col=bandcol[centromere], density=banddens[centromere], border=bandbord[centromere],...)
        qs.semicircle(pos[1]-width, pos[2]-bandpos[centromere+1,2], width,
                   bandpos[centromere+1,2]-bandpos[centromere+1,1],
                   3, col=bandcol[centromere+1], density=banddens[centromere+1], border=bandbord[centromere+1],...)
        centromere.size=0.6*0.5*width/xinch(1)
        symbols(pos[1]-0.5*width, pos[2]-bandpos[centromere,2],circles=1,inches=centromere.size, add=TRUE,fg=gray(bleacher(0)),bg="white",...)
      }
      if (legend) text(pos[1]+0.5*width,pos[2]-(bandpos[,1]+bandpos[,2])/2,paste(chromdata[,"arm"],chromdata[,"band"],sep=""),adj=c(0,0.5),srt=0,cex=cex.leg,...)
    }
  } else {
    warning(paste("Chromosome",chrom,"is not plotted because cytoband data is not available"))
  }
}

grid.chromosome<-function (chrom, side=1, units="hg19", chrom.width=0.5, length.out, 
                           bands="major", legend = c("chrom","band","none"), cex.leg=0.7, 
                           bleach = 0, ...)
{
  bleacher<-function(x) { (x * (1-bleach)) + bleach}
  if (class(units)=="data.frame") {
    chrombands<-.convertUCSCcytoband(units)
  } else{
    chrombands<-.getChrombands(units)
  }
  side<-max(1,min(side,4))
  legend<-match.arg(legend)
  chrom<-switch(as.character(chrom),
         "98"="X",
         "99"="Y",
         as.character(chrom))
    #if (new)
    #    grid.newpage()
  chromdata <- subset(chrombands, chrombands$chr == chrom)
  if (nrow(chromdata)>0){
    lc <- nchar(chromdata$band)
    sel <- !(substr(chromdata$band, lc, lc) %in% letters)
    if (bands != "major")
        sel <- !sel
    chromdata <- chromdata[sel, ]
    rm(lc, sel)
    type.b<-match(chromdata$stain,c("acen","gneg", "gpos", "gvar", "stalk","gpos25","gpos50","gpos75","gpos100","gpos33", "gpos66"))
    bandpos<-chromdata[,c("segstart","segend")]
    bandcol<-gray(bleacher(c(0.5,1,0.2,0.6,0.75,0.7,0.5,0.3,0.1,0.6,0.4)))[type.b]
    banddens<-c(30,-1,-1,-1,10,-1,-1,-1,-1,-1,-1)[type.b]
    bandbord<-gray(bleacher(c(0,0,0,0,1,0,0,0,0,0,0)))[type.b]
    if (!missing(length.out)) {
      bandpos<-(bandpos/max(bandpos))*length.out
    }
    n<-nrow(chromdata)
    centromere<-which(chromdata$arm[-n]!=chromdata$arm[-1])
    if (length(centromere==1)) {
      idx<-c(2:(centromere-1), (centromere+2):(n-1))
    } else {
      idx<-c(2:(n-1))
    }
    if (side %in% 1:2) {
      pos.bottom<-0
      pos.top<-chrom.width
      pos.chrom<-(1+chrom.width)/2
      pos.band<-chrom.width+0.1*(1-chrom.width)
    } else {                                                        
      pos.bottom<-1-chrom.width
      pos.top<-1
      pos.chrom<-(1-chrom.width)/2
      pos.band<-0.1*(1-chrom.width)
    }
    bleachblack<- gray(bleacher(0))
    if (side %in% c(1,3)) {
      pushViewport(viewport(xscale = c(bandpos[1,1] , bandpos[n,2] ), yscale = c(0, 1), clip = "on",...))
      grid.rect(x = bandpos[idx,1], y = pos.bottom, width = bandpos[idx,2] -
          bandpos[idx,1], height = chrom.width, just = c("left",
          "bottom"), default.units = "native", gp = gpar(fill = bandcol[idx],col=bandbord[idx]))
      qs.grid.semicircle(bandpos[1,2], pos.bottom, chrom.width,
          bandpos[1,2] - bandpos[1,1], 2, default.units="native", gp=gpar(fill = bandcol[1],col=bandbord[1]))
      qs.grid.semicircle(bandpos[n,1], pos.bottom, chrom.width,
          bandpos[n,2] - bandpos[n,1], 4, default.units="native", gp=gpar(fill = bandcol[n],col=bandbord[n]))
      if (length(centromere==1)) {
        qs.grid.semicircle(bandpos[centromere,1], pos.bottom, chrom.width, 
            bandpos[centromere,2] - bandpos[centromere,1],
            4, default.units="native", gp=gpar(fill = bandcol[centromere],col=bandbord[centromere]))
        qs.grid.semicircle(bandpos[centromere + 1,2], pos.bottom, chrom.width, 
            bandpos[centromere + 1,2] - bandpos[centromere + 1,1], 
            2, default.units="native", gp=gpar(fill = bandcol[centromere+1],col=bandbord[centromere+1]))
        grid.circle(bandpos[centromere,2], pos.bottom+chrom.width/2, unit(chrom.width*0.3,"npc"), default.units="native", gp = gpar(col=bleachblack, fill="white", lwd=2))
      }
      
      
      if (legend=="chrom") {
        grid.text(chrom, unit(0.5, "npc"), unit(pos.chrom,"native"), gp = gpar(cex = cex.leg))
      } else if (legend=="band") {
        grid.text(paste(chromdata[,"arm"],chromdata[,"band"],sep=""),(bandpos[,1]+bandpos[,2])/2,pos.band,default.units="native",hjust=0,vjust=0.5,rot=90,gp=gpar(cex=cex.leg))    
      }
    } else {
      pushViewport(viewport(xscale = c(0, 1), yscale = c(bandpos[n,2], bandpos[1,1] ), clip = "on",...))
      grid.rect(x = pos.bottom, y = bandpos[idx,1], width = chrom.width, height = bandpos[idx,2] -
          bandpos[idx,1], just = c("left", "bottom"), default.units = "native", gp = gpar(fill = bandcol[idx],col=bandbord[idx]))
      qs.grid.semicircle( pos.bottom, bandpos[1,2],chrom.width,
          bandpos[1,2] - bandpos[1,1],  1, default.units="native", gp=gpar(fill = bandcol[1],col=bandbord[1]))
      qs.grid.semicircle( pos.bottom, bandpos[n,1], chrom.width,
          bandpos[n,2] - bandpos[n,1], 3, default.units="native", gp=gpar(fill = bandcol[n],col=bandbord[n]))
      if (length(centromere==1)) {
        qs.grid.semicircle( pos.bottom, bandpos[centromere,1],  chrom.width,
            bandpos[centromere,2] - bandpos[centromere,1], 
            3, default.units="native", gp=gpar(fill = bandcol[centromere],col=bandbord[centromere]))
        qs.grid.semicircle( pos.bottom, bandpos[centromere + 1,2],      chrom.width,
            bandpos[centromere + 1,2] - bandpos[centromere + 1,1],  
            1, default.units="native", gp=gpar(fill = bandcol[centromere+1],col=bandbord[centromere+1]))
        grid.circle(pos.bottom+chrom.width/2, bandpos[centromere,2],  unit(chrom.width*0.3,"npc"), default.units="native", gp = gpar(col=bleachblack, fill="white", lwd=2))
      }
      if (legend=="chrom") {
        grid.text(chrom, unit(pos.chrom,"native"), unit(0.5, "npc"), gp = gpar(cex = cex.leg))
      } else if (legend=="band") {
        grid.text(paste(chromdata[,"arm"],chromdata[,"band"],sep=""),pos.band,(bandpos[,1]+bandpos[,2])/2,default.units="native",hjust=0,vjust=0.5,rot=0,gp=gpar(cex=cex.leg))    
      }
    
    }
    popViewport()
  } else {
    warning(paste("Chromosome",chrom,"is not plotted because cytoband data is not available"))
  
  }
}

lengthChromosome<-function(chrom, units="hg19") {
  if (class(units)=="data.frame") {
    chrombands<-.convertUCSCcytoband(units)
  } else{
    chrombands<-.getChrombands(units)
  }
  chrom<-characterCHR(chrom)
  chromdata<-subset(chrombands, chrombands$chr %in% chrom)
  if (nrow(chromdata)==0) {
    warning("chromosomes not found in chromosomal data")
    res=rep(NA,length(chrom))
  }else{
    chromlengths<-aggregate(chromdata$segend,list(chromdata$chr),function(x) x[length(x)])
    res<-chromlengths[match(chrom,chromlengths[,1]),2]
  }
  names(res)<-chrom
  return(res)
}


position2Cytoband<-function(chrom,position,units="hg19",bands=c("major","minor")) {
  if (class(units)=="data.frame") {
    chrombands<-.convertUCSCcytoband(units)
  } else{
    chrombands<-.getChrombands(units)
  }
  chrom<-switch(as.character(chrom),
         "98"="X",
         "99"="Y",
         as.character(chrom))
  bands<-match.arg(bands)
  chromdata<-subset(chrombands, chrombands$chr==chrom)
  if (nrow(chromdata)==0) stop("invalid chromosome:",chrom)
  lc<-nchar(chromdata$band)
  sel<-!(substr(chromdata$band,lc,lc) %in% letters)
  if (bands=="minor") 
    sel<-!sel
  chromdata<-chromdata[sel,]
  rm(lc,sel)
  res<-NULL
  for (pos1 in position) {         
    cb<-which(pos1>=chromdata$segstart & pos1<=chromdata$segend)
    if (length(cb)>0)         
      res<-c(res,paste(chromdata[cb,c(1,6,5)],collapse=""))
    else {
      warning("Position ",pos1," not valid for chromosome ",chrom,". It should be between ",chromdata$segstart[1]," and ",chromdata$segend[nrow(chromdata)])
      res<-c(res,"-")
    }
  }
  res
}


drawSimpleChrom<-function(x,y,len=3,width=1,fill,col,orientation=c("h","v"),centromere.size=0.6) {
  # put a simple drawing of a chromosome p:q = 1:2
  # events can be indictaed by fill and col fill=c("a","p","q","p1","p2","p3","q1","q2","q3")
  bandpos<-cbind(c(0,1,2,3,4,7),c(1,2,3,4,7,8))*len/8
  n<-nrow(bandpos)
  centromere<-3
  idx<-c(2:(centromere-1), (centromere+2):(n-1))
  bandcol=rep("white",6)
  if (!missing(fill)) if(length(fill)>0) for (i in 1:length(fill)) {
    if (fill[i]=="a") bandcol[1:6]<-col[i]
    else if (fill[i]=="p") bandcol[1:3]<-col[i]
    else if (fill[i]=="q") bandcol[4:6]<-col[i]
    else if (fill[i]=="p1") bandcol[3]<-col[i]
    else if (fill[i]=="p2") bandcol[2]<-col[i]
    else if (fill[i]=="p3") bandcol[1]<-col[i]
    else if (fill[i]=="q1") bandcol[4]<-col[i]
    else if (fill[i]=="q2") bandcol[5]<-col[i]
    else if (fill[i]=="q3") bandcol[6]<-col[i]
  }  
  banddens=rep(-1,6)
  if (orientation[1]=="h") {
    # draw the inside filling
    rect(x+bandpos[idx,1],y+0.5*width,x+bandpos[idx,2],y-0.5*width, col=bandcol[idx], density=banddens[idx], border=NA)
    qs.semicircle(x+bandpos[1,2], y-0.5*width, width,
               bandpos[1,2]-bandpos[1,1], 2, col=bandcol[1], density=banddens[1], border=NA)
    qs.semicircle(x+bandpos[n,1], y-0.5*width, width,
               bandpos[n,2]-bandpos[n,1], 4, col=bandcol[n], density=banddens[n], border=NA)
    qs.semicircle(x+bandpos[centromere,1], y-0.5*width, width,
               bandpos[centromere,2]-bandpos[centromere,1],
               4, col=bandcol[centromere], density=banddens[centromere], border=NA)
    qs.semicircle(x+bandpos[centromere+1,2], y-0.5*width, width,
               bandpos[centromere+1,2]-bandpos[centromere+1,1],
               2, col=bandcol[centromere+1], density=banddens[centromere+1], border=NA)
    # draw the circumference
    for (i in idx) {
      lines(x+bandpos[i,1:2],rep(y+0.5*width,2),col=1)
      lines(x+bandpos[i,1:2],rep(y-0.5*width,2),col=1)
    }
    qs.semicircle(x+bandpos[1,2], y-0.5*width, width, bandpos[1,2]-bandpos[1,1], 2, col=1, plottype="line")
    qs.semicircle(x+bandpos[n,1], y-0.5*width, width, bandpos[n,2]-bandpos[n,1], 4, col=1, plottype="line")
    qs.semicircle(x+bandpos[centromere,1], y-0.5*width, width, bandpos[centromere,2]-bandpos[centromere,1], 4, col=1, plottype="line")
    qs.semicircle(x+bandpos[centromere+1,2], y-0.5*width, width, bandpos[centromere+1,2]-bandpos[centromere+1,1], 2, col=1, plottype="line")
    # draw the centromere
    centromere.size=centromere.size*0.5*width/yinch(1)
    symbols(x+bandpos[centromere,2], y,circles=1,inches=centromere.size, add=TRUE,fg="black",bg="white")
  } else {
    # draw the inside filling
    rect(x+0.5*width,y-bandpos[idx,1],x-0.5*width,y-bandpos[idx,2], col=bandcol[idx], density=banddens[idx], border=NA)
    qs.semicircle(x-0.5*width, y-bandpos[1,2], width, bandpos[1,2]-bandpos[1,1],
               3, col=bandcol[1], density=banddens[1], border=NA)
    qs.semicircle(x-0.5*width, y-bandpos[n,1], width, bandpos[n,2]-bandpos[n,1], 
               1, col=bandcol[n], density=banddens[n], border=NA)
    qs.semicircle(x-0.5*width, y-bandpos[centromere,1], width, bandpos[centromere,2]-bandpos[centromere,1],
               1, col=bandcol[centromere], density=banddens[centromere], border=NA)
    qs.semicircle(x-0.5*width, y-bandpos[centromere+1,2], width, bandpos[centromere+1,2]-bandpos[centromere+1,1],
               3, col=bandcol[centromere+1], density=banddens[centromere+1], border=NA)
    # draw the circumference
    for (i in idx) {
      lines(rep(x+0.5*width,2),y-bandpos[i,1:2],col=1)
      lines(rep(x-0.5*width,2),y-bandpos[i,1:2],col=1)
    }
    qs.semicircle(x-0.5*width, y-bandpos[1,2], width, bandpos[1,2]-bandpos[1,1], 3, col=1, plottype="line")
    qs.semicircle(x-0.5*width, y-bandpos[n,1], width, bandpos[n,2]-bandpos[n,1], 1, col=1, plottype="line")
    qs.semicircle(x-0.5*width, y-bandpos[centromere,1], width, bandpos[centromere,2]-bandpos[centromere,1], 1, col=1, plottype="line")
    qs.semicircle(x-0.5*width, y-bandpos[centromere+1,2], width, bandpos[centromere+1,2]-bandpos[centromere+1,1], 3, col=1, plottype="line")
    # draw the centromere
    centromere.size=centromere.size*0.5*width/xinch(1)
    symbols(x, y-bandpos[centromere,2],circles=1,inches=centromere.size, add=TRUE,fg="black",bg="white")
  }
}


#########################################################################################

region <- read.table("SIGNIFICANT_REGIONS.bed")
region$V1 <- gsub("chr","",region$V1)
CHR <- region$V1
MapInfo <- region$V2
pdf(file="genome.pdf",height=8,width=8)
if (build == "hg19") {
        chrompos<-prepareGenomePlot(data.frame(CHR,MapInfo),paintCytobands = TRUE,organism="hsa",sexChromosomes = TRUE, unit = "hg19")
}
if (build == "hg18") {
        chrompos<-prepareGenomePlot(data.frame(CHR,MapInfo),paintCytobands = TRUE,organism="hsa",sexChromosomes = TRUE, unit = "hg18")
}
if (build == "mm10") {
        chrompos<-prepareGenomePlot(data.frame(CHR,MapInfo),paintCytobands = TRUE,organism="mmu",sexChromosomes = TRUE, unit = "mm10")
}
if (build == "mm9") {
        #temp <- tempfile(fileext = ".txt.gz")
        #download.file("http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/cytoBand.txt.gz",temp)
        #mm9cytobands <- read.table(temp,sep="\t")
	mm9cytobands <- read.table(paste(path_w4CSeq, "/w4cseq/lib/mm9/cytoBand.txt.gz", sep=""),sep="\t")
        chrompos<-prepareGenomePlot(data.frame(CHR,MapInfo),paintCytobands = TRUE,organism="mmu",sexChromosomes = TRUE, unit = mm9cytobands)
}
rect(chrompos[,2], chrompos[,1]+0.1, chrompos[,2]+region$V3-region$V2, chrompos[,1]+0.3, col="red", border = "red")
dev.off()

png(file="genome.png", width = 8, height = 8, units = 'in', res = 300)
if (build == "hg19") {
        chrompos<-prepareGenomePlot(data.frame(CHR,MapInfo),paintCytobands = TRUE,organism="hsa",sexChromosomes = TRUE, unit = "hg19")
}
if (build == "hg18") {
        chrompos<-prepareGenomePlot(data.frame(CHR,MapInfo),paintCytobands = TRUE,organism="hsa",sexChromosomes = TRUE, unit = "hg18")
}
if (build == "mm10") {
        chrompos<-prepareGenomePlot(data.frame(CHR,MapInfo),paintCytobands = TRUE,organism="mmu",sexChromosomes = TRUE, unit = "mm10")
}
if (build == "mm9") {
        #temp <- tempfile(fileext = ".txt.gz")
        #download.file("http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/cytoBand.txt.gz",temp)
        #mm9cytobands <- read.table(temp,sep="\t")
	mm9cytobands <- read.table(paste(path_w4CSeq, "/w4cseq/lib/mm9/cytoBand.txt.gz", sep=""),sep="\t")
        chrompos<-prepareGenomePlot(data.frame(CHR,MapInfo),paintCytobands = TRUE,organism="mmu",sexChromosomes = TRUE, unit = mm9cytobands)
}
rect(chrompos[,2], chrompos[,1]+0.1, chrompos[,2]+region$V3-region$V2, chrompos[,1]+0.3, col="red", border = "red")
dev.off()

## generate domainogram
if (build == "hg19" || build == "hg18" || build == "mm10") {
	bait_ch_len <- lengthChromosome(sub("chr", "", bait_ch), build)
}
if (build == "mm9") {
	#temp <- tempfile(fileext = ".txt.gz")
	#download.file("http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/cytoBand.txt.gz",temp)
	#mm9cytobands <- read.table(temp,sep="\t")
	mm9cytobands <- read.table(paste(path_w4CSeq, "/w4cseq/lib/mm9/cytoBand.txt.gz", sep=""), sep="\t")
	bait_ch_len <- lengthChromosome(sub("chr", "", bait_ch), mm9cytobands)
	# remove temp file
	#unlink(temp)
}


#system("/var/www/html/w4cseq/bin/bedtools2-2.25.0/bin/mergeBed -i all_reads.bed -c 1 -o count -d 0 > all_sites.bed")
#system("cat all_sites.bed | awk '$4 > 1' > all_sites_noiseRemoved.bed")
system(paste("cat ", file_sort_merge_filter2_realign_bed, " | awk '{if($1==\"", bait_ch, "\")print$1\"\t\"$2\"\t\"$3\"\t1\"}' > intra_cut.bed", sep=""))
system(paste("cat ", enzyme_genome, " | awk '$1==\"", bait_ch, "\"' > intra_all_enzyme.bed", sep=""))
system(paste(path_bedtools, "/intersectBed -a intra_all_enzyme.bed -b intra_cut.bed -v > intra_no_cut.bed", sep=""))
system("cat intra_cut.bed intra_no_cut.bed | sort -k1,1 -k2,2n > intra_sites.bed")


system(paste("cat SIGNIFICANT_REGIONS.bed | awk '$1==\"", bait_ch, "\"' > intra_domains.bed", sep=""))
system(paste("cat ", path_w4CSeq, "/w4cseq/lib/", build, "/", build, "_GENE_sorted.bed | awk '{if($1==\"", bait_ch, "\")print$2\"\t\"$3\"\t\"$6}' > genes.txt", sep=""))

#quick implementation of the running sum
running.sum <- function( x, n ){
  sum.v <- c(0,cumsum(x))
  diff(sum.v,n)
}	

#calculate the binomial p-value with a large window size to calculate the p
p.binom.variable <- function( x, window=20, large.window=3000){
  p.large <- running.sum(x > 0, large.window)/large.window
  first <- p.large[1]
  last  <- tail(p.large,1)
  p.large <- c(rep(first,large.window/2), p.large, rep(last,large.window/2-1))
  p.large <- tail(p.large, n=-window/2)
  p.large <- head(p.large, n=-window/2+1)
  
  
  #calculate the Z score
  X <- running.sum(x > 0, window)
  p.val <- pbinom(X-1, window, p.large, low=F)
  p.val
}

domainogram.bpspace.binom <- function(signal, position, window=200, plot=T, offset=NA, add=0, max.p=10, min.p=0){
  
  if(plot){
    if(is.na(offset))
      plot(c(0,max(position)), c(0,window), type='n', axes=F, xlab="", ylab="") 
    else
      plot(c(0,max(position)), c(0,offset), type='n', axes=F, xlab="", ylab="") 
  }		
  w <- length(signal)
  for(i in 10:window){
    #cat(paste("Starting to plot window size ", i, "\n", sep=""))
    p.val <- p.binom.variable(signal,window=i)
    n <- floor(i/2)
    if(i %% 2){
      x <- position[n:(w-n-1)]
    }else{
      x <- position[n:(w-n)]
    }
    p.val <- -log10(p.val)
    p.val[p.val > max.p] <- max.p
    p.val <- p.val - min.p
    p.val[p.val < 0] <- 0
    new.max.p <- max.p - min.p
    red  <- ifelse(p.val > new.max.p/2, 1, (p.val)/(new.max.p/2))
    green <- ifelse(p.val > new.max.p/2, (p.val-new.max.p/2)/(new.max.p/2), 0)
    col <- rgb(red,green,0)
    points(x,rep(i,length(x))+add, pch='.', col=col)
  }
  #cat("Finished\n")
}

#draw an x-axis given the positions of the fragment ends
drawXAxis <- function(pos, ...){
  chrom.len <- max(pos)
  mb <- chrom.len/1e6
  ats <- c(seq(0,chrom.len, by=10*1e6),chrom.len)
  label <- c(seq(0,mb, by=10),floor(mb))
  axis(1,at=ats, lab=label, cex.lab=1.5, cex.axis=1.5, ...)
}	

drawYAxis <- function(window) {
  axis(2,at=c(10,window),lab=c(10,window), cex.lab=1.5, cex.axis=1.5)
}

data <- read.table("intra_sites.bed")

pdf("domainogram.pdf", 26, 3)
par(mar=c(6,5,5,3))
domainogram.bpspace.binom(signal=data[,4] > 0, position=data[,2], plot=T)
drawXAxis(data[,2])
drawYAxis(200)
points(bait_st, 208, pch=25, cex=1.5,bg="black")
title(ylab="Window size", cex.lab=1.5)
dev.off()

png("domainogram.png", wid=2600, hei=300)
par(mar=c(6,5,5,3))
domainogram.bpspace.binom(signal=data[,4] > 0, position=data[,2], plot=T)
drawXAxis(data[,2])
drawYAxis(200)
points(bait_st, 208, pch=25, cex=2.5,bg="black")
title(ylab="Window size", cex.lab=2.5)
dev.off()


#draw a tubular shaped chromosome
drawChrom <- function(chrom.len, max.wid = 100e6, hei=10, chrom.wid=1, y.loc=5){
  r <- seq(0,2*pi, len=1000)
  chrom.wid = chrom.wid/2
  
  dim.val <- par("din")
  left  <- r[501:1000]
  right <- r[1:500]
  
  correction <- ((max.wid*chrom.wid)/(2*hei)) / (dim.val[1]/dim.val[2])
  
  x <- sin(left)*correction + correction
  y <- cos(left)*chrom.wid+y.loc
  
  x1 <- x
  y1 <- y
  
  x <- c(x, sin(right)*correction+chrom.len -correction)
  y <- c(y, cos(right)*chrom.wid+y.loc )
  
  x2 <- sin(right)*correction+chrom.len -correction
  y2 <- cos(right)*chrom.wid+y.loc
  polygon(x,y, col='white')
  
  col.seq <- seq(0.5,1,len=250)
  col.seq <- c(col.seq,rev(col.seq))
  cols <- rgb(col.seq,col.seq,col.seq)
  segments(x1,y1,x2,rev(y2), col = cols)
  polygon(x,y, lwd=2)
  
}	


#function for drawing two chromosomes
drawLocalChrom <- function(labels=c("1","2"), yloc1 = 3, yloc2 = 5, wid = 2, num.chrom=1, chrom.len){
  plot(c(0,chrom.len), c(yloc1-wid,yloc2+wid), type='n', axes=F, xlab="", ylab="", cex.lab=3)
  
  mb <- chrom.len/1e6
  #mb10 <- floor(mb/10)*10
  #draw an axis
  label <- c(seq(0,mb, by=10),floor(mb))
  #segments(label*1e6, -1e9, label*1e6, 1e9, lwd=2, lty=2, col='grey90')
  #mid.y <- (yloc1+yloc2)/2
  #segments(label*1e6, mid.y-0.15, label*1e6, mid.y+0.15, lwd=2)
  #segments(0, mid.y, mb*1e6, mid.y, lwd=2)
  
  
  #active X
  if(num.chrom==2){
    drawChrom(chrom.len=chrom.len, max.wid=chrom.len, y.loc=yloc1)
    text(-1e6,yloc1, labels[2],cex=1)
  }	
  #inactive X 
  drawChrom(chrom.len=chrom.len, max.wid=chrom.len, y.loc=yloc2)
  
  text(-1e6,yloc2, labels[1],cex=1)
  at <- 0:mb
  #axis(1, at=at*1e6, labels=NA, cex.axis=1, lwd=1,las=1)
  axis(1, at=label*1e6, labels=label, cex.axis=2.5, lwd=2,las=1)
  #axis(3, at=at*1e6, labels=NA, cex.axis=1, lwd=1,las=1)
  #axis(3, at=label*1e6, labels=label, cex.axis=1, lwd=2)
}

#draw the splines showing the interactions
drawSplines.domain <- function( dom, vp.loc, y.base, y.arc, plot=F, col='black', chrom.size=166e6, relative=F){
  if(nrow(dom) == 0)
    return
  for(i in 1:nrow(dom)){
    #start <- min(dom[i,1], dom[i,2])
    #end <- max(dom[i,1], dom[i,2])
    start <- dom[i,1]
    end <- dom[i,2]
    xspline(c(vp.loc, (end + vp.loc)/2, end, start, (end + vp.loc)/2, vp.loc), c(y.base,y.arc,y.base,y.base,y.arc,y.base), open=F, shape=c(0,1,0,0,1,0), col=col, border=col)
  }
}

#dom:       matrix or data.frame with two columns, start and end position of the interactions
#vp.loc:    the position of the viewpoint
#color:     color of the interaction splines
#gene:      data.frame containing the columns for the start, end and strand of gene
#labels:    what should be put on the left side of the plot
#chrom.len: the length of the chromosome

makeSpiderGramSingle <- function( dom, vp.loc, color='black', gene="", labels=c(""), chrom.len ){
  
  #draw two chromosomes
  drawLocalChrom( labels = labels, num.chrom=1, chrom.len = chrom.len )
  
  #and the splines
  drawSplines.domain(dom, vp.loc=vp.loc, plot=F, relative=F, y.arc=8, y.base=5.5, col=color)
  
  #draw the genes
  if(! is.null(nrow(gene))){
    gene <- gene[gene[,2]-gene[,1] < 2e5,]
    rect(gene[,1],5, gene[,2], ifelse(gene[,3]=='+',5.5,4.5), col='black')
  }	
  
}
#cat SIGNIFICANT_REGIONS.bed | awk '$1=="chr5"' > intra_domains.bed
data <- read.table("intra_domains.bed")
dom <- data[,c(2,3)]
#cat mm10_GENE_sorted.bed | awk '{if($1=="chr5")print$2"\t"$3"\t"$6}' > genes.txt
genes<- read.table("genes.txt")

pdf("spider.pdf", 26, 3)
makeSpiderGramSingle(dom=dom, vp.loc=as.numeric(bait_st), color='purple', gene=genes, chrom.len=bait_ch_len)
dev.off()

png("spider.png", wid=2600, hei=300)
makeSpiderGramSingle(dom=dom, vp.loc=as.numeric(bait_st), color='purple', gene=genes, chrom.len=bait_ch_len)
dev.off()


#generate a summary report
#total_reads <- read.table("all_interact.bed",header=FALSE)
#distal_reads <- read.table("distal_interact.bed",header=FALSE)
sink("summary_report.txt")
cat("Here is a summary report of your W4CSEQ result\n")
sink()
system("wc fastq_convert.fq | awk '{print \"Total number of bait-containing reads = \",$1/4}' >> summary_report.txt")
system("wc FASTQ_FILTERED.fq | awk '{print \"Total number of bait-containing reads with good base quality = \",$1/4}' >> summary_report.txt")
#system("wc all_interact.bed | awk '{print \"Total number of interaction reads after removing randomly aligned reads = \",$1}' >> summary_report.txt")
#system("wc distal_interact.bed | awk '{print \"The number of distal interaction reads after removing randomly aligned reads = \",$1}' >> summary_report.txt")
system("wc DISTAL_INTERACTION_SITES.bed | awk '{print \"The number of distal interacting sites = \",$1}' >> summary_report.txt")
system("wc SIGNIFICANT_REGIONS.bed | awk '{print \"The number of significant interacting regions = \",$1}' >> summary_report.txt")


#generate a gene distance distribution plot
system(paste(path_bedtools, "/closestBed -a DISTAL_INTERACTION_SITES.bed -b ", path_w4CSeq, "/w4cseq/lib/", build, "/", build, "_GENE_sorted.bed -D a > 4C_GENE.txt",sep=""))
system(paste(path_bedtools, "/closestBed -a ", enzyme_genome, " -b ", path_w4CSeq, "/w4cseq/lib/",build,"/",build,"_GENE_sorted.bed -D a | awk '$1 != \"chrY\"' > All_GENE.txt",sep=""))

system(paste(path_bedtools, "/closestBed -a DISTAL_INTERACTION_SITES.bed -b ", path_w4CSeq, "/w4cseq/lib/",build,"/",build,"_TSS_sorted.bed -D a > 4C_TSS.txt",sep=""))
system(paste(path_bedtools, "/closestBed -a ", enzyme_genome," -b ", path_w4CSeq, "/w4cseq/lib/",build,"/",build,"_TSS_sorted.bed -D a | awk '$1 != \"chrY\"' > All_TSS.txt",sep=""))

system(paste(path_bedtools, "/closestBed -a DISTAL_INTERACTION_SITES.bed -b ", path_w4CSeq, "/w4cseq/lib/",build,"/",build,"_TTS_sorted.bed -D a > 4C_TTS.txt",sep=""))
system(paste(path_bedtools, "/closestBed -a ", enzyme_genome," -b ", path_w4CSeq, "/w4cseq/lib/",build,"/",build,"_TTS_sorted.bed -D a | awk '$1 != \"chrY\"' > All_TTS.txt",sep=""))

system(paste(path_bedtools, "/closestBed -a DISTAL_INTERACTION_SITES.bed -b ", path_w4CSeq, "/w4cseq/lib/",build,"/",build,"_CPG_sorted.bed -D a > 4C_CPG.txt",sep=""))
system(paste(path_bedtools, "/closestBed -a ", enzyme_genome," -b ", path_w4CSeq, "/w4cseq/lib/",build,"/",build,"_CPG_sorted.bed -D a | awk '$1 != \"chrY\"' > All_CPG.txt",sep=""))

GENE_4C<-read.table("4C_GENE.txt",header=FALSE)
TSS_4C<-read.table("4C_TSS.txt",header=FALSE)
TTS_4C<-read.table("4C_TTS.txt",header=FALSE)
CPG_4C<-read.table("4C_CPG.txt",header=FALSE)

GENE_all<-read.table("All_GENE.txt",header=FALSE)
TSS_all<-read.table("All_TSS.txt",header=FALSE)
TTS_all<-read.table("All_TTS.txt",header=FALSE)
CPG_all<-read.table("All_CPG.txt",header=FALSE)

pdf(file="distance.pdf",height=8,width=8)

par(mfrow = c(2,2))
marks<-c(-500000, -400000, -300000, -200000, -100000, 0, 100000, 200000, 300000, 400000, 500000)

max_gene=max(max(density(GENE_4C$V11,n=16384,bw=10000)$y),max(density(GENE_all$V11,n=16384,bw=10000)$y))*1.2
plot(density(GENE_4C$V11,n=16384,bw=10000), xaxt='n', xlim=c(-500000,500000), ylim=c(0,max_gene), xlab="Distance to the closest gene", ylab="Density", col="red", lwd=4,main="Distance Distribution to Reference Genes")
axis(1,at=marks,labels=format(marks,scientific=FALSE))
lines(density(GENE_all$V11,n=16384,bw=10000),col='blue',lwd=4)
legend("topright",c("4C enzyme sites", "Genome enzyme sites"), cex=0.7, fill=c("red","blue"))

max_tss=max(max(density(TSS_4C$V11,n=16384,bw=10000)$y),max(density(TSS_all$V11,n=16384,bw=10000)$y))*1.2
plot(density(TSS_4C$V11,n=16384,bw=10000), xaxt='n', xaxs="i", xlim=c(-500000,500000), ylim=c(0,max_tss), xlab="Distance to the closet TSS", ylab="Density", col="red", lwd=4,main="Distance Distribution to TSS")
axis(1,at=marks,labels=format(marks,scientific=FALSE))
lines(density(TSS_all$V11,n=16384,bw=10000),col='blue',lwd=4)
legend("topright",c("4C enzyme sites", "Genome enzyme sites"), cex=0.7, fill=c("red","blue"))

max_tts=max(max(density(TTS_4C$V11,n=16384,bw=10000)$y),max(density(TTS_all$V11,n=16384,bw=10000)$y))*1.2
plot(density(TTS_4C$V11,n=16384,bw=10000), xaxt='n', xaxs='i', xlim=c(-500000,500000),ylim=c(0,max_tts),xlab="Distance to the closest TTS",ylab="Density", col="red", lwd=4,main="Distance Distribution to TTS")
axis(1,at=marks,labels=format(marks,scientific=FALSE))
lines(density(TTS_all$V11,n=16384,bw=10000),col='blue',lwd=4)
legend("topright",c("4C enzyme sites", "Genome enzyme sites"), cex=0.7, fill=c("red","blue"))

max_cpg=max(max(density(CPG_4C$V9,n=16384,bw=10000)$y),max(density(CPG_all$V9,n=16384,bw=10000)$y))*1.2
plot(density(CPG_4C$V9,n=16384,bw=10000), xaxt='n', xaxs='i', xlim=c(-500000,500000),ylim=c(0,max_cpg),xlab="Distance to the closest CpG islands",ylab="Density", col="red", lwd=4,main="Distance Distribution to CpG islands")
axis(1,at=marks,labels=format(marks,scientific=FALSE))
lines(density(CPG_all$V9,n=16384,bw=10000),col='blue',lwd=4)
legend("topright",c("4C enzyme sites", "Genome enzyme sites"), cex=0.7, fill=c("red","blue"))
dev.off()

png(file="distance.png", width = 8, height = 8, units = 'in', res = 300)
par(mfrow = c(2,2))
marks<-c(-500000, -400000, -300000, -200000, -100000, 0, 100000, 200000, 300000, 400000, 500000)

max_gene=max(max(density(GENE_4C$V11,n=16384,bw=10000)$y),max(density(GENE_all$V11,n=16384,bw=10000)$y))*1.2
plot(density(GENE_4C$V11,n=16384,bw=10000), xaxt='n', xlim=c(-500000,500000), ylim=c(0,max_gene), xlab="Distance to the closest gene", ylab="Density", col="red", lwd=4,main="Distance Distribution to Reference Genes")
axis(1,at=marks,labels=format(marks,scientific=FALSE))
lines(density(GENE_all$V11,n=16384,bw=10000),col='blue',lwd=4)
legend("topright",c("4C enzyme sites", "Genome enzyme sites"), cex=0.7, fill=c("red","blue"))

max_tss=max(max(density(TSS_4C$V11,n=16384,bw=10000)$y),max(density(TSS_all$V11,n=16384,bw=10000)$y))*1.2
plot(density(TSS_4C$V11,n=16384,bw=10000), xaxt='n', xaxs="i", xlim=c(-500000,500000), ylim=c(0,max_tss), xlab="Distance to the closet TSS", ylab="Density", col="red", lwd=4,main="Distance Distribution to TSS")
axis(1,at=marks,labels=format(marks,scientific=FALSE))
lines(density(TSS_all$V11,n=16384,bw=10000),col='blue',lwd=4)
legend("topright",c("4C enzyme sites", "Genome enzyme sites"), cex=0.7, fill=c("red","blue"))

max_tts=max(max(density(TTS_4C$V11,n=16384,bw=10000)$y),max(density(TTS_all$V11,n=16384,bw=10000)$y))*1.2
plot(density(TTS_4C$V11,n=16384,bw=10000), xaxt='n', xaxs='i', xlim=c(-500000,500000),ylim=c(0,max_tts),xlab="Distance to the closest TTS",ylab="Density", col="red", lwd=4,main="Distance Distribution to TTS")
axis(1,at=marks,labels=format(marks,scientific=FALSE))
lines(density(TTS_all$V11,n=16384,bw=10000),col='blue',lwd=4)
legend("topright",c("4C enzyme sites", "Genome enzyme sites"), cex=0.7, fill=c("red","blue"))

max_cpg=max(max(density(CPG_4C$V9,n=16384,bw=10000)$y),max(density(CPG_all$V9,n=16384,bw=10000)$y))*1.2
plot(density(CPG_4C$V9,n=16384,bw=10000), xaxt='n', xaxs='i', xlim=c(-500000,500000),ylim=c(0,max_cpg),xlab="Distance to the closest CpG islands",ylab="Density", col="red", lwd=4,main="Distance Distribution to CpG islands")
axis(1,at=marks,labels=format(marks,scientific=FALSE))
lines(density(CPG_all$V9,n=16384,bw=10000),col='blue',lwd=4)
legend("topright",c("4C enzyme sites", "Genome enzyme sites"), cex=0.7, fill=c("red","blue"))
dev.off()


#Generate DNA replication analysis
RD_frame<-data.frame(var=1:500000)
RD.tot = 10;
if(build =="hg19" || build =="hg18") {
	cell_type<-c("ESC_BG01","ESC_BG02r1","ESC_BG02r2","ESC_H7","ESC_H9","iPSC_4r1","iPSC_4r2","iPSC_5r1","iPSC_5r2","NPC_BG01r1")
}
if(build =="mm10" || build =="mm9") {
	cell_type<-c("ESC_46C","ESC_D3","ESC_TT2","iPSC","iPSC_1D4","iPSC_2D4","EPL_D3","EBM3_D3","EpiSC5","EpiSC7") 
}
for(RD in (1:RD.tot)){
	system(paste("cp ", path_w4CSeq, "/w4cseq/lib/", build, "/RD/RD", RD ,"_GSM*.bed RD.bed",sep=""))
 	RD_data<-read.table("RD.bed",header=FALSE)	
 	system(paste(path_bedtools, "/windowBed -a RD.bed -b SIGNIFICANT_REGIONS.bed -u -w 0 > RD_4C.bed", sep=""))
 	RD_4C_data<-read.table("RD_4C.bed",header=FALSE)
 	if(RD==1){
		RD_frame<-data.frame(c(RD_4C_data$V4,rep(NA,nrow(RD_frame)-length(RD_4C_data$V4))),c(RD_data$V4,rep(NA,nrow(RD_frame)-length(RD_data$V4))))
	}
 	if(RD> 1){
		RD_frame<-data.frame(RD_frame,c(RD_4C_data$V4,rep(NA,nrow(RD_frame)-length(RD_4C_data$V4))),c(RD_data$V4,rep(NA,nrow(RD_frame)-length(RD_data$V4))))
	}
}
RD_frame$var<-NULL
pdf(file="DNA_replication.pdf",height=8,width=8)
par(pin=c(6.4,4))
par(mar=c(6,5,5,3))
boxplot(RD_frame,at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29),col=c("red","grey"),ylim=c(-6,6),xaxt='n',frame.plot=FALSE,ylab="Log2 transformed loess normalized early/late replication timing ratio",main="DNA replication timing of interacting regions", cex.lab=1, cex.axis=1.5,cex.main=1.5)
axis(1,seq(1,RD.tot)*3-1.5,cell_type,las=2, cex.axis=1)
legend("topright",c("Significant 4C interacting regions", "The whole genome"), cex=0.8, fill=c("red","grey"))
dev.off()

png(file="DNA_replication.png", width = 8, height = 8, units = 'in', res = 300)
par(pin=c(6.4,4))
par(mar=c(6,5,5,3))
boxplot(RD_frame,at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29),col=c("red","grey"),ylim=c(-6,6),xaxt='n',frame.plot=FALSE,ylab="DNA replication timing",main="DNA replication timing of interacting regions", cex.lab=1, cex.axis=1.5,cex.main=1.5)
axis(1,seq(1,RD.tot)*3-1.5,cell_type,las=2, cex.axis=1)
legend("topright",c("Significant 4C interacting regions", "The whole genome"), cex=0.8, fill=c("red","grey"))
dev.off()


#generate ChIP-Seq enrichment analysis
enzyme_genome_100k <- "enzyme_genome_100k.bedsort0"
system(paste("sort -R ", enzyme_genome, " | awk '$1 != \"chrY\"' | head -n 100000 | sort -k1,1 -k2,2n > enzyme_genome_100k.bedsort0", sep=""))
if(chipdata == "yes" && file.info("chip_name.txt")$size > 0) {
	chip_frame <- data.frame(var2=1:100000)
	control = "no"
	chip.tot = 0
	chip_file = ""	

	if (file.exists("chipc.bed")) {
		control = "yes"
	}
	
	if (control == "yes") {
		chipfile <- c("chipc","chip1","chip2","chip3","chip4","chip5","chip6","chip7","chip8","chip9","chip10")
		for (i in 1:11) {
			chip_file <- paste(chipfile[i], ".bed", sep="")
			if (file.exists(chip_file)) {
 				system(paste(path_bedtools, "/windowBed -a DISTAL_INTERACTION_SITES.bed -b ", chipfile[i], ".bed -c > ", chipfile[i], ".txt", sep=""))
				system(paste(path_bedtools, "/windowBed -a ", enzyme_genome_100k, " -b ", chipfile[i], ".bed -c > ", chipfile[i],"_rand.txt", sep=""))
				system(paste("wc -l ", chipfile[i], ".bed | awk '{print $1}' > number", sep=""))
				nrow <- read.table("number")[1,1]
					
				chip_name<-paste(chipfile[i],".txt",sep="")
 				chip_rand_name<-paste(chipfile[i],"_rand.txt",sep="")
 				chip_table<-read.table(chip_name,header=FALSE,sep="\t")
 				chip_rand_table<-read.table(chip_rand_name,header=FALSE,sep="\t")
 		
				if(is.na(chip_table[1,5]) == "FALSE" && i==1) {
					chipc_value <- chip_table[,5]*10000000/nrow
	 				chipc_rand_value <- chip_rand_table[,5]*10000000/nrow
	 				#control="yes"
	 				next
				}

 				if(is.na(chip_table[1,5]) == "FALSE") {
					chip_value <- chip_table[,5]*10000000/nrow
                                	chip_rand_value <- chip_rand_table[,5]*10000000/nrow
					chip_frame<-data.frame(chip_frame,c(chip_value-chipc_value,rep(NA,nrow(chip_frame)-length(chip_table[,5]))),c(chip_rand_value-chipc_rand_value,rep(NA,nrow(chip_frame)-length(chip_rand_table[,5]))))
					chip.tot<-chip.tot+1
				}
			}
		}
	}
	if (control == "no") {
		chipfile <- c("chip1","chip2","chip3","chip4","chip5","chip6","chip7","chip8","chip9","chip10")
		for (i in 1:10) {
			chip_file <- paste(chipfile[i], ".bed", sep="")
                        if (file.exists(chip_file)) {
                        	system(paste(path_bedtools, "/windowBed -a DISTAL_INTERACTION_SITES.bed -b ", chipfile[i], ".bed -c > ", chipfile[i], ".txt", sep=""))
                        	system(paste(path_bedtools, "/windowBed -a ", enzyme_genome_100k, " -b ", chipfile[i], ".bed -c > ", chipfile[i],"_rand.txt", sep=""))
                		system(paste("wc -l ", chipfile[i], ".bed | awk '{print $1}' > number", sep=""))
                        	nrow <- read.table("number")[1,1]

                        	chip_name<-paste(chipfile[i],".txt",sep="")
                        	chip_rand_name<-paste(chipfile[i],"_rand.txt",sep="")
                        	chip_table<-read.table(chip_name,header=FALSE,sep="\t")
                        	chip_rand_table<-read.table(chip_rand_name,header=FALSE,sep="\t")
                

                        	if(is.na(chip_table[1,5]) == "FALSE") {
					chip_value <- chip_table[,5]*10000000/nrow
                                	chip_rand_value <- chip_rand_table[,5]*10000000/nrow
                                	chip_frame<-data.frame(chip_frame,c(chip_value,rep(NA,nrow(chip_frame)-length(chip_table[,5]))),c(chip_rand_value,rep(NA,nrow(chip_frame)-length(chip_rand_table[,5]))))
					chip.tot<-chip.tot+1
                       		}
                	}
	
		}
	}
	
	chip_frame$var2<-NULL
	chip_type=read.table("chip_name.txt",header=FALSE)
	pdf(file="ChIP-Seq.pdf",height=8,width=8)
	par(pin=c(6.4,6))
	par(mar=c(6,5,5,3))
	boxplot(chip_frame,col=c("red","grey"),at=setdiff(1:(3*chip.tot), seq(0,3*chip.tot,3)),ylim=c(-45,45),xaxt='n',frame.plot=FALSE,ylab="ChIP-Seq Tag Density",main="ChIP-Seq signal intensity around interacting sites")
	axis(1,seq(1,3*chip.tot,3)+0.5,chip_type$V1,las=2)
	legend("topright",c("Significant 4C sites", "Random sites"), cex=0.8, fill=c("red","grey"), bg="white")
	dev.off()

	png(file="ChIP-Seq.png", width = 8, height = 8, units = 'in', res = 300)
	par(pin=c(6.4,6))
	par(mar=c(6,5,5,3))
        boxplot(chip_frame,col=c("red","grey"),at=setdiff(1:(3*chip.tot), seq(0,3*chip.tot,3)),ylim=c(-45,45),xaxt='n',frame.plot=FALSE,ylab="ChIP-Seq Tag Density",main="ChIP-Seq signal intensity around interacting sites")
        axis(1,seq(1,3*chip.tot,3)+0.5,chip_type$V1,las=2)
        legend("topright",c("Significant 4C sites", "Random sites"), cex=0.8, fill=c("red","grey"), bg="white")
        dev.off()

}


