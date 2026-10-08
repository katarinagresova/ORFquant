# ORFquant, a software to perform splice-aware 
# quantification of ORF translation using Ribo-seq data
#
# Authors: 
# Lorenzo Calviello (calviello.l.bio@gmail.com)
# Uwe Ohler (Uwe.Ohler@mdc-berlin.de)
#
# This software is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software. If not, see
# <http://www.gnu.org/licenses/>.

# Step 2: the P-sites and the junction reads of a sample.


#' Offset spliced reads on plus strand
#'
#' This function calculates P-sites positions for spliced reads on the plus strand
#' @keywords Ribo-seQC, ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param x a \code{GAlignments} object with a cigar string
#' @param cutoff number representing the offset value
#' @return a \code{GRanges} object with offset reads
#' @seealso \code{\link{prepare_for_ORFquant}}
#' @export

get_ps_fromspliceplus<-function(x,cutoff){
  rang<-cigarillo::cigars_as_ranges_along_ref(cigar(x), lmmpos=start(x),ops="M")
  cs<-lapply(rang,function(x){cumsum(x@width)})
  rangok<-lapply(which(IntegerList(cs)>cutoff),"[[",1)
  rangok<-unlist(rangok)
  gr<-as(x,"GRanges")
  
  ones<-rangok==1
  mores<-rangok>1
  psmores<-GRanges()
  psones<-GRanges()
  
  if(sum(ones)>0){
    psones<-shift(resize(gr[ones],width=1,fix="start"),shift=cutoff)
  }
  
  if(sum(mores)>0){
    rangmore<-rang[mores]
    rangok<-rangok[mores]
    cms<-cumsum(width(rangmore))
    shft<-c()
    for(i in 1:length(rangok)){shft<-c(shft,cutoff-cms[[i]][rangok[i]-1])}
    stt<-start(rangmore)
    stok<-c()
    for(i in 1:length(shft)){stok<-c(stok,stt[[i]][rangok[i]]+shft[i])}
    psmores<-GRanges(IRanges(start=stok,width = 1),seqnames = seqnames(x[mores]),strand=strand(x[mores]),seqlengths=seqlengths(x[mores]))
    
  }
  ps<-sort(c(psones,psmores))
  return(ps)
}

#' Offset spliced reads on minus strand
#'
#' This function calculates P-sites positions for spliced reads on the minus strand
#' @keywords Ribo-seQC, ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param x a \code{GAlignments} object with a cigar string
#' @param cutoff number representing the offset value
#' @return a \code{GRanges} object with offset reads
#' @seealso \code{\link{prepare_for_ORFquant}}
#' @export

get_ps_fromsplicemin<-function(x,cutoff){
  rang<-cigarillo::cigars_as_ranges_along_ref(cigar(x), lmmpos=start(x),ops="M")
  rang<-endoapply(rang,rev)
  cs<-lapply(rang,function(x){cumsum(x@width)})
  rangok<-lapply(which(IntegerList(cs)>cutoff),"[[",1)
  rangok<-unlist(rangok)
  gr<-as(x,"GRanges")
  
  ones<-rangok==1
  mores<-rangok>1
  psmores<-GRanges()
  psones<-GRanges()
  
  
  if(sum(ones)>0){
    psones<-shift(resize(gr[ones],width=1,fix="start"),shift=-cutoff)
  }
  
  
  if(sum(mores)>0){
    rangmore<-rang[mores]
    rangok<-rangok[mores]
    cms<-cumsum(width(rangmore))
    shft<-c()
    for(i in 1:length(rangok)){shft<-c(shft,cutoff-cms[[i]][rangok[i]-1])}
    #start?
    stt<-end(rangmore)
    stok<-c()
    for(i in 1:length(shft)){stok<-c(stok,stt[[i]][rangok[i]]-shft[i])}
    psmores<-GRanges(IRanges(start=stok,width = 1),seqnames = seqnames(x[mores]),strand=strand(x[mores]),seqlengths=seqlengths(x[mores]))
    
  }
  ps<-sort(c(psones,psmores))
  return(ps)
}

#' Prepare the "for_ORFquant" file
#'
#' 
#' @details This function uses a list of pre-determined read lengths, cutoffs and compartments to calculate P_sites positions.\cr
#' Alternatively, bigwig files containing P_sites position for each strand can be specified. Optional bigwig files for uniquely mapping P_sites position (with and without mismatches)
#' can be specified to obtain more statistics on the ORFquant-identified ORFs
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param annotation_file Full path to the annotation file (*Rannot)
#' @param bam_file Full path to the bam file
#' @param chunk_size the number of alignments to read at each iteration, defaults to 5000000, increase when more RAM is available
#' @param path_to_rl_cutoff_file path to the rl_cutoff_file file specifying in 3 columns the read lengths, cutoffs and compartments ("nucl" for standard chromosomes)
#' @param path_to_P_sites_plus_bw path to a bigwig file containing P_sites positions on the plus strand
#' @param path_to_P_sites_minus_bw path to a bigwig file containing P_sites positions on the minus strand
#' @param path_to_P_sites_uniq_plus_bw (Optional) path to a bigwig file containing uniquely mapping P_sites positions on the plus strand
#' @param path_to_P_sites_uniq_minus_bw (Optional) path to a bigwig file containing uniquely mapping P_sites positions on the minus strand
#' @param path_to_P_sites_uniq_mm_plus_bw (Optional) path to a bigwig file containing uniquely mapping (with mismatches) P_sites positions on the plus strand
#' @param path_to_P_sites_uniq_mm_minus_bw (Optional) path to a bigwig file containing uniquely mapping (with mismatches) P_sites positions on the minus strand
#' @param dest_name prefix to use for the output files. Defaults to same as \code{bam_file} (appends "for_ORFquant" to its filename)
#' @param n_cores number of cores to use, each for one chunk of \code{chunk_size} alignments at a time. Defaults to 1
#' @seealso \code{\link{run_ORFquant}}
#' @export

prepare_for_ORFquant<-function(annotation_file,bam_file,path_to_rl_cutoff_file=NA,chunk_size=5000000,path_to_P_sites_plus_bw=NA,
                               path_to_P_sites_minus_bw=NA,path_to_P_sites_uniq_plus_bw=NA,path_to_P_sites_uniq_minus_bw=NA,
                               path_to_P_sites_uniq_mm_plus_bw=NA,path_to_P_sites_uniq_mm_minus_bw=NA,
                               dest_name=NA,n_cores=1){
  
  if(is.na(dest_name)){dest_name=bam_file}
  
  check_n_cores(n_cores)
  
  if(is.na(path_to_rl_cutoff_file) & is.na(path_to_P_sites_plus_bw) & is.na(path_to_P_sites_minus_bw)){
    stop(paste("Please input either the paths to the P_sites bw files, or the path a suitable rl_cutoff table! ", date(),sep=""))
  }
  
  if(!is.na(path_to_rl_cutoff_file) & (!is.na(path_to_P_sites_plus_bw) | !is.na(path_to_P_sites_minus_bw))){
    stop(paste("Please input either the paths to the P_sites bw files, or the path a suitable rl_cutoff table! ", date(),sep=""))
  }
  
  if(xor(is.na(path_to_P_sites_plus_bw),is.na(path_to_P_sites_minus_bw))){
    stop(paste("Please input both path_to_P_sites_plus_bw and path_to_P_sites_minus_bw, or neither! ", date(),sep=""))
  }
  
  if(xor(is.na(path_to_P_sites_uniq_plus_bw),is.na(path_to_P_sites_uniq_minus_bw))){
    stop(paste("Please input both path_to_P_sites_uniq_plus_bw and path_to_P_sites_uniq_minus_bw, or neither! ", date(),sep=""))
  }
  
  if(xor(is.na(path_to_P_sites_uniq_mm_plus_bw),is.na(path_to_P_sites_uniq_mm_minus_bw))){
    stop(paste("Please input both path_to_P_sites_uniq_mm_plus_bw and path_to_P_sites_uniq_mm_minus_bw, or neither! ", date(),sep=""))
  }
  
  check_files_exist(c(annotation_file,bam_file,path_to_rl_cutoff_file,path_to_P_sites_plus_bw,path_to_P_sites_minus_bw,
                      path_to_P_sites_uniq_plus_bw,path_to_P_sites_uniq_minus_bw,
                      path_to_P_sites_uniq_mm_plus_bw,path_to_P_sites_uniq_mm_minus_bw))
  check_output_dir(dest_name)
  
  load_annotation(annotation_file)
  
  if(!is.na(path_to_rl_cutoff_file)){
    rl_cutoff<-read.table(path_to_rl_cutoff_file,header = T,sep = "\t",stringsAsFactors = F)
    
    if(!all(c("read_length","cutoff","compartment")%in%colnames(rl_cutoff))){stop(
      paste("Error: please format the rl_cutoff file correctly, using a tab-separated table with 'read_length', 'cutoff' and 'compartment' as column names! ",date(),sep="")
    )}
    
    if(nrow(rl_cutoff)==0){stop("The rl_cutoff file has no rows: it needs the P-site cutoff of each read length")}
    if(anyNA(suppressWarnings(as.numeric(rl_cutoff$read_length))) || anyNA(suppressWarnings(as.numeric(rl_cutoff$cutoff)))){
      stop("The columns 'read_length' and 'cutoff' of the rl_cutoff file must be numbers")
    }
    #a compartment is "nucl" (all chromosomes that aren't circular) or the name of a (circular) chromosome
    comps_ok<-unique(rl_cutoff$compartment)%in%c("nucl",seqlevels(GTF_annotation$seqinfo))
    if(!any(comps_ok)){
      stop("No compartment of the rl_cutoff file is 'nucl' or the name of a chromosome of the annotation, so no read would give a P-site; its compartments are: ",
           paste(unique(rl_cutoff$compartment),collapse=", "))
    }
    if(!all(comps_ok)){
      warning("These compartments of the rl_cutoff file are not 'nucl' or the name of a chromosome of the annotation, so their rows are not used: ",
              paste(unique(rl_cutoff$compartment)[!comps_ok],collapse=", "))
    }
    
    rl_cutoffs_comp<-split(rl_cutoff,rl_cutoff$compartment)
    compnms<-names(rl_cutoffs_comp)
    
    for(compar in compnms){
      message("Using ",paste(rl_cutoffs_comp[[compar]]$read_length,collapse=","), " nt long footprints with ",paste(rl_cutoffs_comp[[compar]]$cutoff,collapse=",")," as cutoffs, '", compar,"' compartment ... ")
    }
  }
  
  opts <- BamFile(file=bam_file, yieldSize=chunk_size) 
  circs_seq<-seqnames(GTF_annotation$seqinfo)[which(isCircular(GTF_annotation$seqinfo))]
  seqs <- seqinfo(opts)
  if(!any(seqlevels(seqs)%in%seqlevels(GTF_annotation$seqinfo))){
    stop("No chromosome of the BAM file is in the annotation. The BAM file has e.g. ",paste(head(seqlevels(seqs),3),collapse=", "),
         " and the annotation has e.g. ",paste(head(seqlevels(GTF_annotation$seqinfo),3),collapse=", "))
  }
  circs <- seqs@seqnames[which(seqs@seqnames%in%circs_seq)]
  
  param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),what=c("mapq"),tag = "MD")
  seqllll<-seqlevels(GTF_annotation$seqinfo)
  seqleee<-seqlengths(GTF_annotation$seqinfo)
  
  input_P_sites<-GRanges()
  seqlevels(input_P_sites)<-seqllll
  seqlengths(input_P_sites)<-seqleee
  input_P_sites_uniq<-input_P_sites
  input_P_sites_uniq_mm<-input_P_sites
  
  if(!is.na(path_to_P_sites_plus_bw)){
    input_P_sites_mn<-GRanges()
    input_P_sites_pl<-import(path_to_P_sites_plus_bw)
    strand(input_P_sites_pl)<-"+"
    if(!is.na(path_to_P_sites_minus_bw)){
      input_P_sites_mn<-import(path_to_P_sites_minus_bw)
      strand(input_P_sites_mn)<-"-"
    }
    suppressWarnings(input_P_sites<-c(input_P_sites_pl,input_P_sites_mn))
    seqlevels(input_P_sites,pruning.mode="coarse")<-seqllll
    seqlengths(input_P_sites)<-seqleee
    input_P_sites<-sort(input_P_sites)
    
  }
  
  if(!is.na(path_to_P_sites_uniq_plus_bw)){
    input_P_sites_uniq_mn<-GRanges()
    input_P_sites_uniq_pl<-import(path_to_P_sites_uniq_plus_bw)
    strand(input_P_sites_uniq_pl)<-"+"
    if(!is.na(path_to_P_sites_uniq_minus_bw)){
      input_P_sites_uniq_mn<-import(path_to_P_sites_uniq_minus_bw)
      strand(input_P_sites_uniq_mn)<-"-"
    }
    suppressWarnings(input_P_sites_uniq<-c(input_P_sites_uniq_pl,input_P_sites_uniq_mn))
    seqlevels(input_P_sites_uniq,pruning.mode="coarse")<-seqllll
    seqlengths(input_P_sites_uniq)<-seqleee
    input_P_sites_uniq<-sort(input_P_sites_uniq)
    
  }
  
  if(!is.na(path_to_P_sites_uniq_mm_plus_bw)){
    input_P_sites_uniq_mm_mn<-GRanges()
    input_P_sites_uniq_mm_pl<-import(path_to_P_sites_uniq_mm_plus_bw)
    strand(input_P_sites_uniq_mm_pl)<-"+"
    if(!is.na(path_to_P_sites_uniq_mm_minus_bw)){
      input_P_sites_uniq_mm_mn<-import(path_to_P_sites_uniq_mm_minus_bw)
      strand(input_P_sites_uniq_mm_mn)<-"-"
    }
    suppressWarnings(input_P_sites_uniq_mm<-c(input_P_sites_uniq_mm_pl,input_P_sites_uniq_mm_mn))
    seqlevels(input_P_sites_uniq_mm,pruning.mode="coarse")<-seqllll
    seqlengths(input_P_sites_uniq_mm)<-seqleee
    input_P_sites_uniq_mm<-sort(input_P_sites_uniq_mm)
    
  }
  
  #what to do with each chunk (read as alignment file)
  
  yiel<-function(x){
    readGAlignments(x,param = param)
  }
  
  #operations on the chunk (here count reads and whatnot)
  
  mapp<-function(x){
    
    mcols(x)$MD[which(is.na(mcols(x)$MD))]<-"NO"
    
    x<-x[seqnames(x)%in%seqnames(GTF_annotation$seqinfo)]
    seqlevels(x)<-seqlevels(GTF_annotation$seqinfo)
    x_I<-x[grep("I",cigar(x))]
    
    if(length(x_I)>0){
      x<-x[grep("I",cigar(x),invert=T)]
      
    }
    x_D<-x[grep("D",cigar(x))]
    if(length(x_D)>0){
      x<-x[grep("D",cigar(x),invert=T)]
      
    }
    
    # softclipping
    
    clipp <- width(cigarillo::cigars_as_ranges_along_query(x@cigar, ops="S"))
    clipp[elementNROWS(clipp)==0] <- 0
    len_adj <- qwidth(x)-sum(clipp)
    mcols(x)$len_adj <- len_adj
    
    # Remove S from Cigar (read positions/length are already adjusted)
    # it helps calculating P-sites positions for spliced reads
    
    cigg<-cigar(x)
    cigg_s<-grep(cigg,pattern = "S")
    if(length(cigg_s)>0){
      cigs<-cigg[cigg_s]
      cigs<-gsub(cigs,pattern = "^[0-9]+S",replacement = "")
      cigs<-gsub(cigs,pattern = "[0-9]+S$",replacement = "")
      cigg[cigg_s]<-cigs
      x@cigar<-cigg
    }
    mcols(x)$cigar_str<-x@cigar
    x_uniq<-x[x@elementMetadata$mapq>50]
    
    pos<-x[strand(x)=="+"]
    neg<-x[strand(x)=="-"]
    
    
    
    # junctions
    
    juns<-summarizeJunctions(x)
    juns_pos<-juns
    juns_neg<-juns
    mcols(juns_pos)<-NULL
    mcols(juns_neg)<-NULL
    juns_pos$reads<-juns$plus_score
    juns_neg$reads<-juns$minus_score
    strand(juns_pos)<-"+"
    strand(juns_neg)<-"-"
    juns<-sort(c(juns_pos,juns_neg))
    juns<-juns[juns$reads>0]
    
    uniq_juns<-summarizeJunctions(x_uniq)
    uniq_juns_pos<-uniq_juns
    uniq_juns_neg<-uniq_juns
    mcols(uniq_juns_pos)<-NULL
    mcols(uniq_juns_neg)<-NULL
    uniq_juns_pos$reads<-uniq_juns$plus_score
    uniq_juns_neg$reads<-uniq_juns$minus_score
    strand(uniq_juns_pos)<-"+"
    strand(uniq_juns_neg)<-"-"
    uniq_juns<-sort(c(uniq_juns_pos,uniq_juns_neg))
    uniq_juns<-uniq_juns[uniq_juns$reads>0]
    if(length(juns)>0){
      juns$unique_reads<-0
      mat<-match(uniq_juns,juns)
      juns$unique_reads[mat]<-uniq_juns$reads
    }
    rang_jun<-GTF_annotation$junctions
    rang_jun$reads<-0
    rang_jun$unique_reads<-0
    if(length(juns)>0){
      mat<-match(juns,rang_jun)
      juns<-juns[!is.na(mat)]
      mat<-mat[!is.na(mat)]
      rang_jun$reads[mat]<-juns$reads
      rang_jun$unique_reads[mat]<-juns$unique_reads
    }
    
    
    # P-sites calculation
    all_ps_comps<-GRanges()
    seqlevels(all_ps_comps)<-seqllll
    seqlengths(all_ps_comps)<-seqleee
    
    uniq_ps_comps<-GRanges()
    seqlevels(uniq_ps_comps)<-seqllll
    seqlengths(uniq_ps_comps)<-seqleee
    
    uniq_mm_ps_comps<-GRanges()
    seqlevels(uniq_mm_ps_comps)<-seqllll
    seqlengths(uniq_mm_ps_comps)<-seqleee
    
    if(!is.na(path_to_rl_cutoff_file)){
      list_pss<-list()
      for(comp in names(rl_cutoffs_comp)){
        all_rl_ps<-GRangesList()
        uniq_rl_ps<-GRangesList()
        uniq_rl_mm_ps<-GRangesList()
        
        seqlevels(all_rl_ps)<-seqllll
        seqlevels(uniq_rl_ps)<-seqllll
        seqlevels(uniq_rl_mm_ps)<-seqllll
        
        seqlengths(all_rl_ps)<-seqleee
        seqlengths(uniq_rl_ps)<-seqleee
        seqlengths(uniq_rl_mm_ps)<-seqleee
        
        
        chroms<-comp
        
        if(comp=="nucl"){chroms=seqlevels(x)[!seqlevels(x)%in%circs]}
        resul<-rl_cutoffs_comp[[comp]]
        
        for(i in seq_along(resul$read_length)){
          
          all_ps<-GRangesList()
          uniq_ps<-GRangesList()
          uniq_mm_ps<-GRangesList()
          seqlevels(all_ps)<-seqllll
          seqlevels(uniq_ps)<-seqllll
          seqlevels(uniq_mm_ps)<-seqllll
          
          seqlengths(all_ps)<-seqleee
          seqlengths(uniq_ps)<-seqleee
          seqlengths(uniq_mm_ps)<-seqleee
          
          rl<-as.numeric(resul$read_length[i])
          ct<-as.numeric(resul$cutoff[i])
          ok_reads<-pos[mcols(pos)$len_adj%in%rl]
          ok_reads<-ok_reads[as.vector(seqnames(ok_reads))%in%chroms]
          
          ps_plus<-GRanges()
          seqlevels(ps_plus)<-seqllll
          seqlengths(ps_plus)<-seqleee
          ps_plus_uniq<-ps_plus
          ps_plus_uniq_mm<-ps_plus
          
          if(length(ok_reads)>0){
            unspl<-ok_reads[grep(pattern="N",x=cigar(ok_reads),invert=T)]
            
            ps_unspl<-shift(resize(GRanges(unspl),width=1,fix="start"),shift=ct)
            
            spl<-ok_reads[grep(pattern="N",x=cigar(ok_reads))]
            firstb<-as.numeric(sapply(strsplit(cigar(spl),"M"),"[[",1))
            lastb<-as.numeric(sapply(strsplit(cigar(spl),"M"),function(x){gsub(x[length(x)],pattern="^[^_]*N",replacement="")}))
            firstok<-spl[firstb>ct]
            firstok<-shift(resize(GRanges(firstok),width=1,fix="start"),shift=ct)
            
            lastok<-spl[lastb>=rl-ct]
            lastok<-shift(resize(GRanges(lastok),width=1,fix="end"),shift=-(rl-ct-1))
            
            
            multi<-spl[firstb<=ct & lastb<rl-ct]
            
            
            ps_spl<-GRanges()
            seqlevels(ps_spl)<-seqllll
            seqlengths(ps_spl)<-seqleee
            
            if(length(multi)>0){
              ps_spl<-get_ps_fromspliceplus(multi,cutoff=ct)
              
            }
            mcols(ps_spl)<-mcols(multi)
            
            seqlevels(firstok)<-seqllll
            seqlevels(lastok)<-seqllll
            seqlevels(ps_unspl)<-seqllll
            seqlevels(ps_spl)<-seqllll
            
            seqlengths(firstok)<-seqleee
            seqlengths(lastok)<-seqleee
            seqlengths(ps_unspl)<-seqleee
            seqlengths(ps_spl)<-seqleee
            
            
            ps_plus<-c(ps_unspl,firstok,lastok,ps_spl)
            ps_plus_uniq<-ps_plus[mcols(ps_plus)$mapq>50]
            ps_plus_uniq_mm<-ps_plus[mcols(ps_plus)$mapq>50 & nchar(mcols(ps_plus)$MD)>3]
            
            mcols(ps_plus)<-NULL
            mcols(ps_plus_uniq)<-NULL
            mcols(ps_plus_uniq_mm)<-NULL
            
          }
          ok_reads<-neg[mcols(neg)$len_adj%in%rl]
          ok_reads<-ok_reads[as.vector(seqnames(ok_reads))%in%chroms]
          
          ps_neg<-GRanges()
          seqlevels(ps_neg)<-seqllll
          seqlengths(ps_neg)<-seqleee
          ps_neg_uniq<-ps_neg
          ps_neg_uniq_mm<-ps_neg
          
          if(length(ok_reads)>0){
            unspl<-ok_reads[grep(pattern="N",x=cigar(ok_reads),invert=T)]
            
            ps_unspl<-shift(resize(GRanges(unspl),width=1,fix="start"),shift=-ct)
            
            spl<-ok_reads[grep(pattern="N",x=cigar(ok_reads))]
            
            firstb<-as.numeric(sapply(strsplit(cigar(spl),"M"),"[[",1))
            lastb<-as.numeric(sapply(strsplit(cigar(spl),"M"),function(x){gsub(x[length(x)],pattern="^[^_]*N",replacement="")}))
            lastok<-spl[lastb>ct]
            lastok<-shift(resize(GRanges(lastok),width=1,fix="start"),shift=-ct)
            
            firstok<-spl[firstb>=rl-ct]
            firstok<-shift(resize(GRanges(firstok),width=1,fix="end"),shift=(rl-ct-1))
            
            multi<-spl[firstb<rl-ct & lastb<=ct]
            
            
            ps_spl<-GRanges()
            seqlevels(ps_spl)<-seqllll
            seqlengths(ps_spl)<-seqleee
            
            
            if(length(multi)>0){
              ps_spl<-get_ps_fromsplicemin(multi,cutoff=ct)
            }
            mcols(ps_spl)<-mcols(multi)
            
            seqlevels(firstok)<-seqllll
            seqlevels(lastok)<-seqllll
            seqlevels(ps_unspl)<-seqllll
            seqlevels(ps_spl)<-seqllll
            
            seqlengths(firstok)<-seqleee
            seqlengths(lastok)<-seqleee
            seqlengths(ps_unspl)<-seqleee
            seqlengths(ps_spl)<-seqleee
            
            ps_neg<-c(ps_unspl,firstok,lastok,ps_spl)
            ps_neg_uniq<-ps_neg[mcols(ps_neg)$mapq>50]
            ps_neg_uniq_mm<-ps_neg[mcols(ps_neg)$mapq>50 & nchar(mcols(ps_neg)$MD)>3]
            
            mcols(ps_neg)<-NULL
            mcols(ps_neg_uniq)<-NULL
            mcols(ps_neg_uniq_mm)<-NULL
          }
          
          all_ps<-sort(c(ps_plus,ps_neg))
          uniq_ps<-sort(c(ps_plus_uniq,ps_neg_uniq))
          uniq_mm_ps<-sort(c(ps_plus_uniq_mm,ps_neg_uniq_mm))
          if(length(all_ps)>0){
            ps_res<-unique(all_ps)
            ps_res$score<-countOverlaps(ps_res,all_ps,type="equal")
            all_ps<-ps_res
          }
          if(length(uniq_ps)>0){
            ps_res<-unique(uniq_ps)
            ps_res$score<-countOverlaps(ps_res,uniq_ps,type="equal")
            uniq_ps<-ps_res
            
          }
          if(length(uniq_mm_ps)>0){
            ps_res<-unique(uniq_mm_ps)
            ps_res$score<-countOverlaps(ps_res,uniq_mm_ps,type="equal")
            uniq_mm_ps<-ps_res
            
          }
          all_rl_ps[[as.character(rl)]]<-all_ps
          uniq_rl_ps[[as.character(rl)]]<-uniq_ps
          uniq_rl_mm_ps[[as.character(rl)]]<-uniq_mm_ps
          
          
        }
        #here comps
        list_rlct<-list(all_rl_ps,uniq_rl_ps,uniq_rl_mm_ps)
        names(list_rlct)<-c("P_sites_all","P_sites_uniq","P_sites_uniq_mm")
        list_pss[[comp]]<-list_rlct
      }
      #for rl, merge psites
      
      
      all_ps_comps<-GRangesList()
      seqlevels(all_ps_comps)<-seqllll
      seqlengths(all_ps_comps)<-seqleee
      rls_comps<-unique(unlist(lapply(list_pss,FUN=function(x) names(x[["P_sites_all"]]) )))
      for(rl in rls_comps){
        reads_rl_comp<-GRanges()
        seqlevels(reads_rl_comp)<-seqllll
        seqlengths(reads_rl_comp)<-seqleee
        for(comp in names(list_pss)){
          if(sum(rl%in%names(list_pss[[comp]][["P_sites_all"]]))>0){
            oth<-list_pss[[comp]][["P_sites_all"]][[rl]]
            
            if(!is.null(oth)){
              seqlevels(oth)<-seqllll
              seqlengths(oth)<-seqleee
              reads_rl_comp<-c(reads_rl_comp,oth)
            }
          }          
          all_ps_comps[[rl]]<-reads_rl_comp
        }
        
      }
      
      uniq_ps_comps<-GRangesList()
      seqlevels(uniq_ps_comps)<-seqllll
      seqlengths(uniq_ps_comps)<-seqleee
      rls_comps<-unique(unlist(lapply(list_pss,FUN=function(x) names(x[["P_sites_uniq"]]) )))
      for(rl in rls_comps){
        reads_rl_comp<-GRanges()
        seqlevels(reads_rl_comp)<-seqllll
        seqlengths(reads_rl_comp)<-seqleee
        for(comp in names(list_pss)){
          if(sum(rl%in%names(list_pss[[comp]][["P_sites_uniq"]]))>0){
            oth<-list_pss[[comp]][["P_sites_uniq"]][[rl]]
            
            if(!is.null(oth)){
              seqlevels(oth)<-seqllll
              seqlengths(oth)<-seqleee
              reads_rl_comp<-c(reads_rl_comp,oth)
            }
          }          
          uniq_ps_comps[[rl]]<-reads_rl_comp
        }
        
      }
      
      uniq_mm_ps_comps<-GRangesList()
      seqlevels(uniq_mm_ps_comps)<-seqllll
      seqlengths(uniq_mm_ps_comps)<-seqleee
      rls_comps<-unique(unlist(lapply(list_pss,FUN=function(x) names(x[["P_sites_uniq_mm"]]) )))
      for(rl in rls_comps){
        reads_rl_comp<-GRanges()
        seqlevels(reads_rl_comp)<-seqllll
        seqlengths(reads_rl_comp)<-seqleee
        for(comp in names(list_pss)){
          if(sum(rl%in%names(list_pss[[comp]][["P_sites_uniq_mm"]]))>0){
            oth<-list_pss[[comp]][["P_sites_uniq_mm"]][[rl]]
            if(!is.null(oth)){
              seqlevels(oth)<-seqllll
              seqlengths(oth)<-seqleee
              reads_rl_comp<-c(reads_rl_comp,oth)
            }
          }          
          uniq_mm_ps_comps[[rl]]<-reads_rl_comp
        }
        
      }
      
    }
    
    list_res<-list(all_ps_comps,uniq_ps_comps,uniq_mm_ps_comps,rang_jun)
    names(list_res)<-c("P_sites_all","P_sites_uniq","P_sites_uniq_mm","junctions")
    
    return(list_res)
  }
  
  message("Calculating P-sites positions and junctions ... ", date())
  
  #read the BAM in chunks of chunk_size alignments, same as GenomicFiles::reduceByYield. With
  #n_cores>1, forked processes run mapp() on the chunks, at most n_cores at a time, while the next
  #ones are read. The junction reads are added up chunk by chunk, and the P-sites of all chunks once,
  #below: adding each chunk to the P-sites of the earlier ones took longer with each chunk. The sums
  #don't depend on the chunks.

  open(opts)
  on.exit(if(Rsamtools::isOpen(opts)){close(opts)})
  for_ORFquant<-list()
  chunks_ps<-list()
  jobs<-list()
  eof<-FALSE
  while(!eof || length(jobs)>0){
    if(!eof && length(jobs)<n_cores){
      chunk<-yiel(opts)
      eof<-length(chunk)==0
      if(!eof){jobs[[length(jobs)+1]]<-mc_job(mapp(chunk),n_cores)}
      rm(chunk)
    }
    #the results of the oldest jobs that are done, so that their processes end; wait for the oldest
    #when n_cores are running or all chunks are read
    while(length(jobs)>0){
      res<-mc_value(jobs[[1]],wait=eof || length(jobs)>=n_cores)
      if(is.null(res)){break}
      jobs<-jobs[-1]
      if(length(for_ORFquant)==0){
        for_ORFquant<-res
      }else{
        for_ORFquant$junctions$reads<-for_ORFquant$junctions$reads+res$junctions$reads
        for_ORFquant$junctions$unique_reads<-for_ORFquant$junctions$unique_reads+res$junctions$unique_reads
      }
      chunks_ps[[length(chunks_ps)+1]]<-res[c("P_sites_all","P_sites_uniq","P_sites_uniq_mm")]
    }
  }
  close(opts)
  if(length(chunks_ps)>1){
    for(n in names(chunks_ps[[1]])){
      ps<-lapply(chunks_ps,function(x){unlist(x[[n]])})
      ps<-ps[lengths(ps)>0]
      if(length(ps)>0){for_ORFquant[[n]]<-do.call(c,unname(ps))}
    }
  }
  rm(chunks_ps)
  
  if(length(for_ORFquant$P_sites_all)>0){
    merged_all_ps<-unlist(for_ORFquant$P_sites_all)
    
    if(length(merged_all_ps)>0){
      covv_pl<-coverage(merged_all_ps[strand(merged_all_ps)=="+"],weight = merged_all_ps[strand(merged_all_ps)=="+"]$score)
      covv_pl<-GRanges(covv_pl)
      covv_pl<-covv_pl[covv_pl$score>0]
      covv_min<-coverage(merged_all_ps[strand(merged_all_ps)=="-"],weight = merged_all_ps[strand(merged_all_ps)=="-"]$score)
      covv_min<-GRanges(covv_min)
      covv_min<-covv_min[covv_min$score>0]
      strand(covv_pl)<-"+"
      strand(covv_min)<-"-"
      
      merged_all_ps<-sort(c(covv_pl,covv_min))
      
    }
    for_ORFquant$P_sites_all<-merged_all_ps
  }
  
  
  if(length(for_ORFquant$P_sites_uniq)>0){
    merged_uniq_ps<-unlist(for_ORFquant$P_sites_uniq)
    if(length(merged_uniq_ps)>0){
      
      covv_pl<-coverage(merged_uniq_ps[strand(merged_uniq_ps)=="+"],weight = merged_uniq_ps[strand(merged_uniq_ps)=="+"]$score)
      covv_pl<-GRanges(covv_pl)
      covv_pl<-covv_pl[covv_pl$score>0]
      
      covv_min<-coverage(merged_uniq_ps[strand(merged_uniq_ps)=="-"],weight = merged_uniq_ps[strand(merged_uniq_ps)=="-"]$score)
      covv_min<-GRanges(covv_min)
      covv_min<-covv_min[covv_min$score>0]
      strand(covv_pl)<-"+"
      strand(covv_min)<-"-"
      merged_uniq_ps<-sort(c(covv_pl,covv_min))
      
    }
    for_ORFquant$P_sites_uniq<-merged_uniq_ps
  }
  
  if(length(for_ORFquant$P_sites_uniq_mm)>0){
    merged_uniq_mm_ps<-unlist(for_ORFquant$P_sites_uniq_mm)
    if(length(merged_uniq_mm_ps)>0){
      
      covv_pl<-coverage(merged_uniq_mm_ps[strand(merged_uniq_mm_ps)=="+"],weight = merged_uniq_mm_ps[strand(merged_uniq_mm_ps)=="+"]$score)
      covv_pl<-GRanges(covv_pl)
      covv_pl<-covv_pl[covv_pl$score>0]
      
      covv_min<-coverage(merged_uniq_mm_ps[strand(merged_uniq_mm_ps)=="-"],weight = merged_uniq_mm_ps[strand(merged_uniq_mm_ps)=="-"]$score)
      covv_min<-GRanges(covv_min)
      covv_min<-covv_min[covv_min$score>0]
      strand(covv_pl)<-"+"
      strand(covv_min)<-"-"
      merged_uniq_mm_ps<-sort(c(covv_pl,covv_min))
    }
    for_ORFquant$P_sites_uniq_mm<-merged_uniq_mm_ps
  }
  
  
  if(!is.na(path_to_P_sites_plus_bw) | !is.na(path_to_P_sites_minus_bw) ){
    for_ORFquant$P_sites_all<-input_P_sites
  }
  
  if(!is.na(path_to_P_sites_uniq_plus_bw) | !is.na(path_to_P_sites_uniq_minus_bw) ){
    for_ORFquant$P_sites_uniq<-input_P_sites_uniq
  }
  
  if(!is.na(path_to_P_sites_uniq_mm_plus_bw) | !is.na(path_to_P_sites_uniq_mm_minus_bw) ){
    for_ORFquant$P_sites_uniq_mm<-input_P_sites_uniq_mm
  }
  
  message("Calculating P-sites positions and junctions --- Done! ", date())
  
  if(length(for_ORFquant$P_sites_all)==0){
    stop("No P-sites were found. Check that the read lengths of the rl_cutoff file are those of the reads of the BAM file (",bam_file,")")
  }
  
  save(for_ORFquant,file = paste(dest_name,"for_ORFquant",sep = "_"))
  invisible(paste(dest_name,"for_ORFquant",sep = "_"))
}
