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

# Detection of translated ORFs in a transcript: the ORF search, the choice of the start
# codon, and the multitaper and p-value tests.


#' Find ATG-starting ORFs in a sequence
#'
#' This function loads the annotation created by the \code{prepare_annotation_files function}
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param tx_name transcript_id
#' @param sequence DNAString object containing the sequence of the transcript
#' @param get_all_starts Output all possible start codons? Defaults to \code{TRUE}
#' @param Stop_Stop Also report the parts of a frame that no ORF with a start codon covers, as ORFs without a start codon (\code{type} \code{frame_ok})? Each goes from the start of the transcript or a stop codon to the next stop codon, the first start codon or the end of the transcript. Only frames with an ORF with a start codon have them. Defaults to \code{FALSE}
#' @param scores Deprecated
#' @param genetic_code_table GENETIC_CODE table to use
#' @return \code{GRanges} object containing coordinates for the detected ORFs
#' @seealso \code{\link{detect_translated_orfs}}
#' @export


get_orfs<-function(tx_name,sequence,get_all_starts=T,Stop_Stop=F,scores=c(1,.5),genetic_code_table){
  list_frames<-list()
  length<-nchar(sequence)
  for(u in 0:2){
    pept<-NA
    pept<-unlist(strsplit(as.character(suppressWarnings(translate_solve(subseq(sequence,start=u+1),genetic.code = genetic_code_table))),split=""))
    
    
    starts<-pept=="M"
    
    stops<-pept=="*"
    
    start_pos<-((1:length(pept))[starts])*3
    if(length(start_pos)>0){
      start_pos<-start_pos+u-2
    } else {start_pos<-NA}
    
    stop_pos<-((1:length(pept))[stops])*3-3
    if(length(stop_pos)>0){
      stop_pos<-stop_pos+u
    } else {stop_pos<-NA}
    
    st2vect<-c()
    for(h in 1:length(start_pos)){
      st1<-start_pos[h]
      diff<-stop_pos-st1
      diff<-diff[diff>0]
      if(length(diff)>0){st2<-st1+min(diff)}
      if(length(diff)==0){st2<-NA}
      st2vect[h]<-st2
      
    }
    st_st<-data.frame(cbind(start_pos,st2vect),stringsAsFactors=F)
    st_st<-st_st[!is.na(st_st[,1]),]
    st_st<-st_st[!is.na(st_st[,2]),]
    if(dim(st_st)[1]==0){list_frames[[paste("frame",u,sep = "_")]]<-GRanges(); next}
    if(get_all_starts==T){
      gra_orf<-GRanges(seqnames = paste(tx_name,"frame",u,sep = "_"),strand = "+",ranges = IRanges(start = st_st[,1],end = st_st[,2]))
    }
    if(get_all_starts==F){
      gra_orf<-reduce(GRanges(seqnames = paste(tx_name,"frame",u,sep = "_"),strand = "*",ranges = IRanges(start = st_st[,1],end = st_st[,2])))
    }
    if(length(gra_orf)>0){
      # columns are set on mcols() and put back once: a GRanges $<- also runs updateObject() (~20 ms)
      cols<-mcols(gra_orf)
      cols$type="ORF"
      cols$score<-scores[1]
      mcols(gra_orf)<-cols
    }
    
    gra_par<-GRanges()
    if(Stop_Stop==T){
      stops<-pept=="*"
      stop_pos<-((1:length(pept))[stops])*3
      if(length(stop_pos)>0){
        stop_pos<-stop_pos+u
      } else {stop_pos<-NA}
      
      
      vector_starts<-c()
      vector_stops<-c()
      
      if(sum(stops)==0){vector_starts<-u+1; vector_stops<-length}
      
      if(sum(stops)==1){
        vector_starts[1]<-u+1
        vector_starts[2]<-stop_pos+1
        vector_stops[1]<-stop_pos-3
        vector_stops[2]<-length
        
      }
      
      if(sum(stops)>0){
        for(stoppo in 0:(length(stop_pos))){
          start<-stop_pos[stoppo]+1
          if(stoppo==0){
            start<-u+1
          }
          end<-stop_pos[stoppo+1]-3
          if(stoppo==length(stop_pos)){
            end<-length
          }
          
          if((end-start)<5){
            next
          }
          vector_starts[stoppo+1]<-start
          vector_stops[stoppo+1]<-end
          
        }
        
        
      }
      
      st_st<-data.frame(cbind(vector_starts,vector_stops))
      st_st<-st_st[!is.na(st_st[,"vector_starts"]),]
      st_st<-st_st[!is.na(st_st[,"vector_stops"]),]
      #st_st<-st_st[(st_st[,2]-st_st[,1]>11),]
      gra_par<-GRanges(seqnames = paste(tx_name,"frame",u,sep = "_"),strand = "+",ranges = IRanges(start = st_st[,1],end = st_st[,2]))
      if(length(gra_par)>0){
        gra_par<-setdiff(gra_par,gra_orf)
        
      }
      if(length(gra_par)>0){
        gra_par$type="frame_ok"
        gra_par$score<-scores[2]                        
      }
    }
    gra<-c(gra_orf,gra_par)
    
    list_frames[[paste("frame",u,sep = "_")]]<-gra
  }
  GRangesList(list_frames)
}

#' Extract output from multitaper analysis of a signal
#'
#' This function uses the multitaper tool to extract F-values and multitaper spectral coefficients
#' @details Values reported correspond to the closest frequency to 1/3 (same parameters as in RiboTaper). \cr
#' Padding to a minimum length of 1024 is performed to increase spectral resolution.
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param x numeric signal to analyze
#' @param n_tapers n of tapers to use
#' @param time_bw time_bw parameter
#' @param slepians_values set of calculated slepian functions to use in the multitaper analysis
#' @return two numeric values representing the F-value for the multitaper test and its corresponding spectral coefficient at the closest frequency to 1/3
#' @seealso \code{\link{detect_translated_orfs}}, \code{\link[multitaper]{spec.mtm}}, \code{\link[multitaper]{dpss}}
#' @export

take_Fvals_spect<-function(x,n_tapers,time_bw,slepians_values){
  if(length(x)<25){
    remain<-50-length(x)
    x<-c(rep(0,as.integer(remain/2)),x,rep(0,remain%%2+as.integer(remain/2)))
  }
  if(length(x)<1024/2){padding<-1024}
  if(length(x)>=1024/2){padding<-"default"}
  resSpec1 <- spec.mtm(as.ts(x), k=n_tapers, nw=time_bw, nFFT = padding, centreWithSlepians = TRUE, Ftest = TRUE, maxAdaptiveIterations = 100,returnZeroFreq=F,plot=F,dpssIN=slepians_values)
  
  Fmax_3nt<-resSpec1$mtm$Ftest[which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))]
  spect_3nt<-resSpec1$spec[which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))]
  return(c(Fmax_3nt,spect_3nt))
  
}

#' Select start codon 
#'
#' This function selects the start codon for ORFs in the same transcript
#' @details ORFs are divided based on stop codon and Ribo-seq signal between start codons is used to select one.\cr
#' When more than \code{cutoff_ave} fraction of codons is in-frame between two candidate start codons, the most upstream
#' is selected.
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param ORFs Set of detected ORFs
#' @param P_sites_rle Rle signal of P_sites along the transcript
#' @param cutoff cutoff of total in-frame signal between start codons (sensitive to outliers). Defaults to NA
#' @param cutoff_ave cutoff for frequency of in-frame codons between two start codons (less sensitive to outliers). Defaults to .5
#' @return Set of detected ORFs, including info about the possible longest ORF for that frame.
#' @seealso \code{\link{detect_translated_orfs}}, \code{\link{get_orfs}}
#' @export

select_start<-function(ORFs,P_sites_rle,cutoff=NA,cutoff_ave=.5){
  ORFs<-ORFs[order(end(ORFs),start(ORFs),decreasing = F)]
  names(ORFs)<-NULL
  longest_ORF<-split(ORFs,end(ORFs))
  maxo<-which.max(width(longest_ORF))
  longest_ORF<-unlist(longest_ORF[splitAsList(unname(maxo), names(maxo))])
  # the P_sites of each ORF are a slice of one plain vector: P_sites_rle[ORFs], an RleList
  # with one element per ORF, costs S4 calls per ORF
  psit_all<-as.vector(P_sites_rle)
  st<-start(ORFs)
  en<-end(ORFs)
  ok<-sapply(seq_along(ORFs),function(i){sum(psit_all[st[i]:en[i]]>0)>2})
  if(length(ok)==0){return(GRanges())}
  ORFs<-ORFs[ok]
  st<-st[ok]
  en<-en[ok]
  rrr<-lapply(seq_along(ORFs),function(i){
    fr<-suppressWarnings(matrix(psit_all[st[i]:en[i]],nrow = 3))
    fr<-fr[,colSums(fr)>0,drop=F]
    if(dim(fr)[2]==0){return(GRanges())}
    fra<-apply(fr,2,function(y){y/sum(y)})
    infr_freq<-rowMeans(fra)[1]
    infr<-((rowSums(fr))[1])/sum(fr)
    c(round(infr_freq*100,digits = 4),round(infr*100,digits = 4))
  })
  # one DataFrame for all ORFs: building one per ORF and rbind-ing them is slow
  rrr<-matrix(as.numeric(unlist(rrr)),ncol = 2,byrow = T)
  mcols(ORFs)<-DataFrame(ave_pct_fr=rrr[,1],pct_fr=rrr[,2],ave_pct_fr_st=rep(NA,nrow(rrr)),pct_fr_st=rep(NA,nrow(rrr)))
  
  if(!is.na(cutoff)){
    ORFs<-ORFs[ORFs$pct_fr>=cutoff]
  }
  if(!is.na(cutoff_ave)){
    ORFs<-ORFs[ORFs$ave_pct_fr>=cutoff_ave]
  }
  if(length(ORFs)==0){return(GRanges())}
  # the ORFs with the same stop are indices into ORFs (sorted by end and start), and the
  # ORF kept for each stop is taken once at the end: endoapply() on split(ORFs) costs S4
  # calls per stop
  st<-start(ORFs)
  en<-end(ORFs)
  ave_pct_fr_st<-ORFs$ave_pct_fr_st
  pct_fr_st<-ORFs$pct_fr_st
  keep<-integer(0)
  for(x in split(seq_along(ORFs),en)){
    if(length(x)==1){keep<-c(keep,x);next}
    okorfa<-c()
    for(cnt in 1:length(x)){
      psit<-psit_all[st[x[cnt]]:en[x[cnt]]]
      
      if(cnt<length(x)){
        psit<-psit[1:(st[x[cnt+1]]-st[x[cnt]])]
      }
      if(cnt==length(x)){okorfa<-cnt}
      
      frames<-suppressWarnings(matrix(as.vector(psit),nrow = 3))
      frames<-frames[,colSums(frames)>0,drop=F]
      if(dim(frames)[2]==0){next}
      frames<-apply(frames,2,function(y){y/sum(y)})
      infr_freq<-rowMeans(frames)[1]
      infr<-sum(psit[seq(1,length(psit),by=3)])/sum(psit)
      if(!is.na(cutoff)){
        if(infr>=cutoff & !is.na(infr)){
          ave_pct_fr_st[x[cnt]]<-round(infr_freq,digits = 4)
          pct_fr_st[x[cnt]]<-round(infr,digits = 4)
          okorfa<-cnt
          break
        }
      }
      
      if(!is.na(cutoff_ave)){
        if(infr_freq>=cutoff_ave & !is.na(infr_freq)){
          ave_pct_fr_st[x[cnt]]<-round(infr_freq*100,digits = 4)
          pct_fr_st[x[cnt]]<-round(infr*100,digits = 4)
          okorfa<-cnt
          break
        }
        
      }
      
      
    }
    keep<-c(keep,x[okorfa])
  }
  # columns are set on mcols() and put back once: a GRanges $<- also runs updateObject() (~20 ms)
  cols<-mcols(ORFs)
  cols$ave_pct_fr_st<-ave_pct_fr_st
  cols$pct_fr_st<-pct_fr_st
  mcols(ORFs)<-cols
  ORFs<-ORFs[keep]
  mcols(ORFs)$longest_ORF<-longest_ORF[match(end(ORFs),names(longest_ORF))]
  return(ORFs)
  
}

#' Collect ORF Ribo-seq statistics
#'
#' This function calculates statistics for the analysis of P_sites profiles for each ORF
#' @details Number of P_sites (uniquely mapping or all), frame percentage and multitaper test 
#' statistics are collected for each ORF. The parameter space for the multitaper analysis was explored in the RiboTaper paper.
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param ORFs Set of detected ORFs
#' @param P_sites_rle Rle signal of P_sites along the transcript
#' @param P_sites_uniq_rle Rle signal of uniquely mapping P_sites along the transcript
#' @param P_sites_uniq_mm_rle Rle signal of uniquely mapping P_sites with mismatches along the transcript
#' @param cutoff cutoff of the fraction of the ORF's P_sites that are in frame (\code{pct_fr}): only ORFs above it, with P_sites at more than 2 positions, are tested. Defaults to .5
#' @param tapers Number of tapers to use in the multitaper analysis. Defaults to 24
#' @param bw time_bw parameter to use in the multitaper analysis. Defaults to 12
#' @return Set of detected ORFs, with the columns \code{pval}, \code{pval_uniq}, \code{P_sites_raw}, \code{P_sites_raw_uniq}, \code{P_sites_raw_uniq_mm}, \code{pct_fr} (a fraction, replacing the percentage from \code{select_start}) and \code{ORF_id_tr}; \code{pval} is \code{NA} for ORFs not tested.
#' @seealso \code{\link{detect_translated_orfs}}, \code{\link{get_orfs}}, \code{\link{take_Fvals_spect}}
#' @export

calc_orf_pval<-function(ORFs,P_sites_rle,P_sites_uniq_rle,P_sites_uniq_mm_rle,cutoff=.5,tapers=24,bw=12){
  # columns are filled as plain vectors and set once at the end: assigning
  # one element of a GRanges column (ORFs$pval[i]<-) costs ~6 ms
  pval<-rep(NA,length(ORFs))
  pval_uniq<-rep(NA,length(ORFs))
  P_sites_raw<-rep(NA,length(ORFs))
  P_sites_raw_uniq<-rep(NA,length(ORFs))
  P_sites_raw_uniq_mm<-rep(NA,length(ORFs))
  pct_fr<-rep(NA,length(ORFs))
  ORF_id_tr<-ORFs$ORF_id_tr
  # the P_sites of each ORF are a slice of plain vectors, and its coordinates are taken
  # once: ORFs[i] and its accessors cost S4 calls per ORF
  psit_all<-as.vector(P_sites_rle)
  psit_uniq_all<-as.vector(P_sites_uniq_rle)
  psit_uniq_mm_all<-as.vector(P_sites_uniq_mm_rle)
  st<-start(ORFs)
  en<-end(ORFs)
  sn<-as.character(seqnames(ORFs))
  
  
  for(i in 1:length(ORFs)){
    psit<-psit_all[st[i]:en[i]]
    psit_uniq<-psit_uniq_all[st[i]:en[i]]
    psit_uniq_mm<-psit_uniq_mm_all[st[i]:en[i]]
    P_sites_raw[i]<-sum(psit)
    P_sites_raw_uniq[i]<-sum(psit_uniq)
    
    P_sites_raw_uniq_mm[i]<-sum(psit_uniq_mm)
    ORF_id_tr[i]<-paste(sn[i],st[i],en[i],sep = "_")
    if(sum(psit)>0){
      infr<-round(sum(psit[seq(1,length(psit),by=3)])/sum(psit),digits = 4)
      pct_fr[i]<-infr
    }
    if(sum(psit>0)>2){
      if(infr>cutoff){
        if(length(psit)<25){slepians<-dpss(n=length(psit)+(50-length(psit)),k=tapers,nw=bw)}
        if(length(psit)>=25){slepians<-dpss(n=length(psit),k=tapers,nw=bw)}
        vals<-take_Fvals_spect(x = psit,n_tapers = tapers,time_bw = bw,slepians_values = slepians)
        pval[i]<-pf(q=vals[1],df1=2,df2=(2*tapers)-2,lower.tail=F)
        vals<-take_Fvals_spect(x = psit_uniq,n_tapers = tapers,time_bw = bw,slepians_values = slepians)
        pval_uniq[i]<-pf(q=vals[1],df1=2,df2=(2*tapers)-2,lower.tail=F)
        
        
      }
    }
  }
  # columns are set on mcols() and put back once: a GRanges $<- also runs updateObject() (~20 ms)
  cols<-mcols(ORFs)
  cols$pval<-pval
  cols$pval_uniq<-pval_uniq
  cols$P_sites_raw<-P_sites_raw
  cols$P_sites_raw_uniq<-P_sites_raw_uniq
  cols$P_sites_raw_uniq_mm<-P_sites_raw_uniq_mm
  cols$pct_fr<-pct_fr
  cols$ORF_id_tr<-ORF_id_tr
  mcols(ORFs)<-cols
  return(ORFs)
}


#' Detect actively translated ORFs 
#'
#' This function detects translated ORFs
#' @details A set of transcripts, together with genome sequence and Ribo-signal are analyzed to extract translated ORFs
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param selected_txs set of selected transcripts, output from \code{select_txs}
#' @param genome_sequence BSgenome object
#' @param annotation Rannot object containing annotation of CDS and transcript structures (see \code{prepare_annotation_files})
#' @param P_sites GRanges object with P_sites positions
#' @param P_sites_uniq GRanges object with uniquely mapping P_sites positions
#' @param P_sites_uniq_mm GRanges object with uniquely mapping (with mismatches) P_sites positions
#' @param genomic_region GRanges object with genomic coordinates of the genomic region analyzed
#' @param genetic_code GENETIC_CODE table to use
#' @param all_starts \code{get_all_starts} parameter for the \code{get_orfs} function
#' @param nostarts \code{Stop_Stop} parameter for the \code{get_orfs} function
#' @param start_sel_cutoff \code{cutoff} parameter for the \code{select_start} function
#' @param start_sel_cutoff_ave \code{cutoff_ave} parameter for the \code{select_start} function
#' @param cutoff_fr_ave \code{cutoff} parameter for the  \code{calc_orf_pval} functions
#' @param uniq_signal Use only signal from uniquely mapping reads? Defaults to \code{FALSE}.
#' @return A list with transcript coordinates, exonic coordinates and statistics for each ORF exonic bin and junction(from \code{select_txs}).\cr\cr
#' The value for each column is as follows:\cr\cr
#' \code{ave_pct_fr}: average percentage of in-frame reads for each codon in the ORF
#' \code{pct_fr}: fraction (0 to 1) of in-frame reads in the ORF
#' \code{ave_pct_fr_st}: average percentage of in-frame reads per each codon between the selected start codon and the next candidate one
#' \code{pct_fr_st}: percentage of in-frame reads between the selected start codon and the next candidate one
#' \code{longest_ORF}: GRanges coordinates for the longest ORF with the same stop codon
#' \code{pval}: P-value for the multitaper F-test at 1/3 using the ORF P_sites profile
#' \code{pval_uniq}: P-value for the multitaper F-test at 1/3 using the ORF P_sites profile (only uniquely mapping reads)
#' \code{P_sites_raw}: Raw number of P_sites mapping to the ORF
#' \code{P_sites_raw_uniq}: Uniquely mapping P_sites mapping to the ORF
#' \code{P_sites_raw_uniq_mm}: Uniquely mapping P_sites with mismatches mapping to the ORF
#' \code{ORF_id_tr}: ORF id containing <tx_id>_<start>_<end>
#' \code{Protein}: AAString sequence of the translated protein
#' \code{region}: Genomic coordinates of the analyzed region
#' \code{gene_id}: gene_id for the corresponding analyzed transcript
#' \code{gene_biotype}: gene biotype for the corresponding analyzed transcript
#' \code{gene_name}: gene name for the corresponding analyzed transcript
#' \code{transcript_id}: transcript_id for the corresponding analyzed ORF
#' \code{transcript_biotype}: transcript biotype for the corresponding analyzed ORF
#' @seealso \code{\link{select_txs}}, \code{\link{get_orfs}}, \code{\link{take_Fvals_spect}}, \code{\link{select_start}}, \code{\link{prepare_annotation_files}}, \code{\link{ORFquant_output}}
#' @export

detect_translated_orfs<-function(selected_txs,genome_sequence,annotation,P_sites,P_sites_uniq,P_sites_uniq_mm,genomic_region,genetic_code,
                                 all_starts=T,nostarts=F,start_sel_cutoff=NA,start_sel_cutoff_ave=.5,
                                 cutoff_fr_ave=.5,uniq_signal=F){
  
  orfs_gr<-list()
  orfs_gen_gr<-GRangesList()
  intr_txs<-annotation$introns_txs
  tr_gen<-annotation$trann
  txs_sels<-unique(unlist(selected_txs$txs_selected))
  annot_sels<-annotation$exons_txs[txs_sels]
  # the selected transcripts' introns, looked up below: a lookup by name in the genome-wide list costs ~20 ms
  intr_sels<-intr_txs[names(intr_txs)%in%txs_sels]
  mapp<-mapToTranscripts(P_sites,annot_sels)
  # mcols(x)$col<- rather than x$col<-: a GRanges $<- also runs updateObject() (~20 ms)
  mcols(mapp)$reads<-P_sites$score[mapp$xHits]
  lens_sels<-sum(width(annot_sels))
  seqm<-seqlengths(mapp)
  seqlengths(mapp)<-lens_sels[match(names(seqm),names(lens_sels))]
  cov_txs<-coverage(mapp,weight = mapp$reads)
  
  mapp<-mapToTranscripts(P_sites_uniq,annot_sels)
  if(length(P_sites_uniq)>0){
    mcols(mapp)$reads<-P_sites_uniq$score[mapp$xHits]
    seqm<-seqlengths(mapp)
    seqlengths(mapp)<-lens_sels[match(names(seqm),names(lens_sels))]
    cov_uniq_txs<-coverage(mapp,weight = mapp$reads)
  }
  
  if(length(P_sites_uniq)==0){cov_uniq_txs<-coverage(mapp)}
  
  mapp<-mapToTranscripts(P_sites_uniq_mm,annot_sels)
  
  if(length(P_sites_uniq_mm)==0){cov_uniq_mm_txs<-coverage(mapp)}
  
  if(length(P_sites_uniq_mm)>0){
    mcols(mapp)$reads<-P_sites_uniq_mm$score[mapp$xHits]
    seqm<-seqlengths(mapp)
    seqlengths(mapp)<-lens_sels[match(names(seqm),names(lens_sels))]
    cov_uniq_mm_txs<-coverage(mapp,weight = mapp$reads)
  }
  txs_seqs<-extractTranscriptSeqs(genome_sequence,annot_sels)
  
  for(tx in txs_sels){
    
    ex_txs<-annot_sels[tx]
    intr_tx<-intr_sels[[tx]]
    #map cds in tx space
    
    
    covtx<-cov_txs[[tx]]
    covtx_uniq<-cov_uniq_txs[[tx]]
    covtx_uniq_mm<-cov_uniq_mm_txs[[tx]]
    
    seq_tx<-txs_seqs[[tx]]
    if(length(covtx)==0){next}
    orfs<-unlist(get_orfs(tx_name = tx,sequence = seq_tx,get_all_starts=all_starts,Stop_Stop = nostarts,genetic_code_table=genetic_code))
    #no need to know frame for now
    #names(orfs)<-NULL
    if(length(orfs)==0){next}
    
    orfs<-GRanges(seqnames = tx,ranges = orfs@ranges,strand = strand(orfs))
    orfs<-select_start(ORFs = orfs,P_sites_rle = covtx,cutoff = start_sel_cutoff,cutoff_ave = start_sel_cutoff_ave)
    if(length(orfs)==0){next}
    orfs<-calc_orf_pval(ORFs = orfs,P_sites_rle = covtx,P_sites_uniq_rle = covtx_uniq,P_sites_uniq_mm_rle = covtx_uniq_mm,cutoff = cutoff_fr_ave)
    if(length(orfs)==0){next}
    #uniq_flag
    if(uniq_signal){
      orfs<-subset(orfs,pval_uniq<.05)
    }
    if(!uniq_signal){
      orfs<-subset(orfs,pval<.05)
    }
    if(length(orfs)==0){next}
    #orfs$gene_id<-mapIds(keys = tx,x = annot,column = "GENEID",keytype = "TXNAME")
    # proteins are collected in an AAStringSet and set once: a GRanges $<- per ORF costs ~20 ms
    Protein<-AAStringSet(rep("NA",length(orfs)))
    for(h in 1:length(orfs)){
      Protein[h]<-AAStringSet(as.character(translate_solve(seq_tx[orfs@ranges[h]],genetic.code = genetic_code)))
    }
    orfs$Protein<-Protein
    #orfs$Protein<-AAStringSet(orfs$Protein)
    orfs_gen<-from_tx_togen(ORFs = orfs,exons = ex_txs,introns = intr_tx)
    #here I should annotate
    
    tr_gen_tx<-tr_gen[as.vector(match(seqnames(orfs)[1],tr_gen[,"transcript_id"])),]
    
    cols<-mcols(orfs)
    cols$region<-genomic_region
    cols$gene_id<-unique(as.character(tr_gen_tx[,"gene_id"]))
    cols$gene_biotype<-unique(as.character(tr_gen_tx[,"gene_biotype"]))
    cols$gene_name<-unique(as.character(tr_gen_tx[,"gene_name"]))
    cols$transcript_id<-unique(as.character(tr_gen_tx[,"transcript_id"]))
    cols$transcript_biotype<-unique(as.character(tr_gen_tx[,"transcript_biotype"]))
    
    #must add the other compatible txs, to avoid calculating same stuff
    cols$compatible_with<-NA
    mcols(orfs)<-cols
    for(w in 1:length(orfs)){
      orf<-orfs[w]
      nam<-orf$ORF_id_tr
      orfs_gr[[nam]]<-orf
      orfs_gen_gr[[nam]]<-orfs_gen[[nam]]
      
    }
    orfs_gr<-GRangesList(orfs_gr)
    
  }
  if(length(orfs_gr)==0){
    return(list())
  }
  tx_orfs<-unique(sapply(strsplit(names(orfs_gr),split="_"),function(x){
    len<-length(x)
    x[-c(len-1,len)]
  }))
  a<-lapply(selected_txs$txs,function(x){x[x%in%tx_orfs]})
  cols<-mcols(selected_txs)
  cols$txs_orfs<-CharacterList(a)
  
  check<-sapply(cols$txs_orfs,FUN = length)
  use<-rep("shared",length(check))
  use[check==1]<-"unique"
  use[check==0]<-"absent"
  use[check>1]<-"shared"
  cols$use_ORFs<-use
  mcols(selected_txs)<-cols
  
  orfs_unq_gr<-list()
  
  for(j in names(orfs_gr)){
    
    orf_gen<-orfs_gen_gr[[j]]
    orf_tx<-orfs_gr[[j]]
    featexs<-selected_txs[selected_txs$type=="E"]
    featjuns<-selected_txs[selected_txs$type=="J"]
    
    over_tx<-featexs[featexs%over%orf_gen]
    if(length(featjuns)>0){
      over_tx<-sort(c(over_tx,featjuns[which(featjuns%in%gaps(orf_gen))]))
    }
    
    a<-sapply(over_tx$txs,function(x){length(x[x%in%as.character(seqnames(orf_tx))])})
    over_tx<-over_tx[a>0]
    
    orfs_unq_gr[[j]]<-over_tx
  }
  orfs_unq_gr<-GRangesList(orfs_unq_gr)
  
  if(length(orfs_gr)>1){
    
    ident_mat<-matrix(FALSE,nrow=length(orfs_gr),ncol=length(orfs_gr))
    for(i in names(orfs_gr)){
      gen<-orfs_gen_gr[[i]]
      
      ident<-c()
      for(j in 1:length(orfs_gen_gr)){
        ident[j]<-identical(gen,orfs_gen_gr[[j]])
      }
      ident_mat[which(names(orfs_gr)==i),]<-ident
    }
    #diag(ident_mat)<-NA
    ident_mat[lower.tri(ident_mat,diag=T)]<-NA
    ide<-which(ident_mat,arr.ind=T)
    if(length(ide)>0){
      rems<-c()
      ide<-split(ide,ide[,1])
      for(i in 1:length(ide)){
        iden<-ide[[i]]
        iden<-iden[!iden==as.numeric(names(ide)[i])]
        if(sum(iden%in%rems)>0){next}
        ind<-names(orfs_gr)[as.numeric(names(ide)[i])]
        ch<-names(orfs_gr)[iden]
        rems<-unique(c(rems,iden))
        #new_cat<-paste(unlist(orfs_gr[ch])$category_tx,collapse=";")
        new_id<-paste(unlist(orfs_gr[ch])$ORF_id_tr,collapse=";")
        orfs_gr[[ind]]$compatible_with<-new_id
        #if(length(unique(unlist(orfs_gr[ch])$category_tx))==1){next}
        #orfs_gr[[ind]]$compatible_categories<-new_cat
        #orfs_gr[[ind]]$category_tx<-"multiple"
      }
      orfs_gr<-orfs_gr[-rems]
      orfs_gen_gr<-orfs_gen_gr[-rems]
      orfs_unq_gr<-orfs_unq_gr[-rems]
    }
  }
  list_res<-list(orfs_gr,orfs_gen_gr,orfs_unq_gr)
  names(list_res)<-c("ORFs_tx_position","ORFs_genomic_position","ORFs_features")
  return(list_res)
}
