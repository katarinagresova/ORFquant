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

# The annotation of the detected ORFs: splice features, transcript and genome space.


#' Annotate splice features of detected ORFs
#'
#' This function detects usage of different exons and exonic boundaries of one ORF with respect to a reference ORF.
#' @details each exon is aligned to the closest one to match acceptor and donor sites, or to annotate missing exons.
#' \code{5ss} and \code{3ss} indicate exon 5' and 3', respectively. \code{CDS_spanning} indicates retained intron;
#' \code{missing_CDS} indicates no overlapping exon (missed or included); \code{monoCDS} indicates a single-exon ORF; 
#' \code{firstCDS} and \code{lastCDS} indicate first CDS exon or last CDS exon.
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param orf_gen Exon structure of a detected ORF
#' @param ref_cds Exon structure of a reference ORF
#' @return Exon structure of detected ORF including possible missing exons from reference, together with a \code{spl_type} column
#' including the annotation for each exon (e.g. alternative acceptors or donor).
#' @seealso \code{\link{detect_translated_orfs}}, \code{\link{annotate_ORFs}}
#' @export

annotate_splicing<-function(orf_gen,ref_cds){
  
  # the overlaps are found once, on orf_gen before it is sorted below: two %over% and a
  # findOverlaps() on the sorted exons cost several ms each
  ov<-findOverlaps(orf_gen,ref_cds)
  overref<-seq_along(orf_gen)%in%queryHits(ov)
  refover<-seq_along(ref_cds)%in%subjectHits(ov)
  spl_ran<-GRanges()
  
  if(sum(!refover)>0){
    spl_ran<-c(spl_ran,ref_cds[!refover])
    refgrl<-ref_cds[!refover]
    # one exon per element, made at once: [[<- copies the GRangesList for each exon
    grliss<-relist(refgrl,IRanges::PartitioningByEnd(seq_along(refgrl)))
    # columns are set on mcols() and put back once: a GRanges $<- also runs updateObject() (~20 ms)
    cols<-mcols(spl_ran)
    cols$ref<-grliss
    cols$spl_type<-"missing_CDS"
    cols$cds_id<-NULL
    cols$cds_name<-NULL
    cols$exon_rank<-NULL
    mcols(spl_ran)<-cols
    
  }
  
  o<-order(orf_gen)
  orf_gen<-sort(orf_gen)
  if(length(orf_gen)>0){
    # overlaps found and exons combined once, not per exon: each GRanges op costs several ms.
    # The hits of ov, with the queries numbered as in the sorted orf_gen (sort() orders as order())
    hq<-match(queryHits(ov),o)
    hs<-subjectHits(ov)
    refs<-list()
    spl_types<-list()
    for(f in 1:length(orf_gen)){
      ran<-orf_gen[f]
      # ref and spl_type are kept, and set on all exons at once below: per exon, GRangesList() and
      # each GRanges $<- cost several ms, and binding the exons' GRangesList columns even more
      ref<-NULL
      spl_type<-NULL
      last_ex<-length(orf_gen)
      if(overref[f]==T){
        #sorted: in ref_cds order, as ref_cds[ref_cds%over%ran]; findOverlaps() may give an exon's hits by position
        ref_over<-ref_cds[sort(hs[hq==f])]
        #annotate for 5' and 3'; porcoddio
        
        if(length(ref_over)>1){
          ref<-ref_over
          spl_type<-"CDS_spanning"
          if(as.vector(strand(orf_gen[1]))=="+"){
            if(start(ran)==min(start(ref_over))){
              if(end(ran)==max(end(ref_over))){
                spl_type<-"CDS_spanning;same_5ss;same_3ss"
                if(f==last_ex){spl_type<-"CDS_spanning;same_5ss;same_lastCDS"}
                if(f==1){spl_type<-"CDS_spanning;same_firstCDS;same_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;same_5monoCDS;same_3monoCDS"}
              }
              if(end(ran)>max(end(ref_over))){
                spl_type<-"CDS_spanning;same_5ss;down_3ss"
                if(f==last_ex){spl_type<-"CDS_spanning;same_5ss;down_lastCDS"}
                if(f==1){spl_type<-"CDS_spanning;same_firstCDS;down_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;same_5monoCDS;down_3monoCDS"}
              }
              if(end(ran)<max(end(ref_over))){
                spl_type<-"CDS_spanning;same_5ss;up_3ss"
                if(f==last_ex){spl_type<-"CDS_spanning;same_5ss;up_lastCDS"}
                if(f==1){spl_type<-"CDS_spanning;same_firstCDS;up_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;same_5monoCDS;up_3monoCDS"}
              }
              
            }
            if(end(ran)==max(end(ref_over))){
              if(start(ran)>min(start(ref_over))){
                spl_type<-"CDS_spanning;down_5ss;same_3ss"
                if(f==last_ex){spl_type<-"CDS_spanning;down_5ss;same_lastCDS"}
                if(f==1){spl_type<-"CDS_spanning;down_firstCDS;same_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;down_5monoCDS;same_3monoCDS"}
              }
              if(start(ran)<min(start(ref_over))){
                spl_type<-"CDS_spanning;up_5ss;same_3ss"
                if(f==last_ex){spl_type<-"CDS_spanning;up_5ss;same_lastCDS"}
                if(f==1){spl_type<-"CDS_spanning;up_firstCDS;same_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;up_5monoCDS;same_3monoCDS"}
              }
              
            }
            
            if(end(ran)>max(end(ref_over))){
              if(start(ran)>min(start(ref_over))){
                spl_type<-"CDS_spanning;down_5ss;down_3ss"
                if(f==last_ex){spl_type<-"CDS_spanning;down_5ss;down_lastCDS"}
                if(f==1){spl_type<-"CDS_spanning;down_firstCDS;down_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;down_5monoCDS;down_3monoCDS"}
              }
              if(start(ran)<min(start(ref_over))){
                spl_type<-"CDS_spanning;up_5ss;down_3ss"
                if(f==last_ex){spl_type<-"CDS_spanning;up_5ss;down_lastCDS"}
                if(f==1){spl_type<-"CDS_spanning;up_firstCDS;down_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;up_5monoCDS;down_3monoCDS"}
              }
              
            }
            
            if(end(ran)<max(end(ref_over))){
              if(start(ran)>min(start(ref_over))){
                spl_type<-"CDS_spanning;down_5ss;up_3ss"
                if(f==last_ex){spl_type<-"CDS_spanning;down_5ss;up_lastCDS"}
                if(f==1){spl_type<-"CDS_spanning;down_firstCDS;up_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;down_5monoCDS;up_3monoCDS"}
              }
              if(start(ran)<min(start(ref_over))){
                spl_type<-"CDS_spanning;up_5ss;up_3ss"
                if(f==last_ex){spl_type<-"CDS_spanning;up_5ss;up_lastCDS"}
                if(f==1){spl_type<-"CDS_spanning;up_firstCDS;up_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;up_5monoCDS;up_3monoCDS"}
              }
              
            }
            
            
            
            
          }
          
          #if - and spanning
          
          if(as.vector(strand(orf_gen[1]))=="-"){
            if(start(ran)==min(start(ref_over))){
              if(end(ran)==max(end(ref_over))){
                spl_type<-"CDS_spanning;same_5ss;same_3ss"
                if(f==1){spl_type<-"CDS_spanning;same_5ss;same_lastCDS"}
                if(f==last_ex){spl_type<-"CDS_spanning;same_firstCDS;same_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;same_5monoCDS;same_3monoCDS"}
              }
              if(end(ran)>max(end(ref_over))){
                spl_type<-"CDS_spanning;up_5ss;same_3ss"
                if(f==1){spl_type<-"CDS_spanning;up_5ss;same_lastCDS"}
                if(f==last_ex){spl_type<-"CDS_spanning;up_firstCDS;same_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;up_5monoCDS;same_3monoCDS"}
              }
              if(end(ran)<max(end(ref_over))){
                spl_type<-"CDS_spanning;down_5ss;same_3ss"
                if(f==1){spl_type<-"CDS_spanning;down_5ss;same_lastCDS"}
                if(f==last_ex){spl_type<-"CDS_spanning;down_firstCDS;same_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;down_5monoCDS;same_3monoCDS"}
              }
              
            }
            if(end(ran)==max(end(ref_over))){
              if(start(ran)>min(start(ref_over))){
                spl_type<-"CDS_spanning;same_5ss;up_3ss"
                if(f==1){spl_type<-"CDS_spanning;same_5ss;up_lastCDS"}
                if(f==last_ex){spl_type<-"CDS_spanning;same_firstCDS;up_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;same_5monoCDS;up_3monoCDS"}
              }
              if(start(ran)<min(start(ref_over))){
                spl_type<-"CDS_spanning;same_5ss;down_3ss"
                if(f==1){spl_type<-"CDS_spanning;same_5ss;down_lastCDS"}
                if(f==last_ex){spl_type<-"CDS_spanning;same_firstCDS;down_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;same_5monoCDS;down_3monoCDS"}
              }
              
            }
            
            if(end(ran)>max(end(ref_over))){
              if(start(ran)>min(start(ref_over))){
                spl_type<-"CDS_spanning;up_5ss;up_3ss"
                if(f==1){spl_type<-"CDS_spanning;up_5ss;up_lastCDS"}
                if(f==last_ex){spl_type<-"CDS_spanning;up_firstCDS;up_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;up_5monoCDS;up_3monoCDS"}
              }
              if(start(ran)<min(start(ref_over))){
                spl_type<-"CDS_spanning;up_5ss;down_3ss"
                if(f==1){spl_type<-"CDS_spanning;up_5ss;down_lastCDS"}
                if(f==last_ex){spl_type<-"CDS_spanning;up_firstCDS;down_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;up_5monoCDS;down_3monoCDS"}
              }
              
            }
            
            if(end(ran)<max(end(ref_over))){
              if(start(ran)>min(start(ref_over))){
                spl_type<-"CDS_spanning;down_5ss;up_3ss"
                if(f==1){spl_type<-"CDS_spanning;down_5ss;up_lastCDS"}
                if(f==last_ex){spl_type<-"CDS_spanning;down_firstCDS;up_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;down_5monoCDS;up_3monoCDS"}
              }
              if(start(ran)<min(start(ref_over))){
                spl_type<-"CDS_spanning;down_5ss;down_3ss"
                if(f==1){spl_type<-"CDS_spanning;down_5ss;down_lastCDS"}
                if(f==last_ex){spl_type<-"CDS_spanning;down_firstCDS;down_3ss"}
                if(1==last_ex){spl_type<-"CDS_spanning;down_5monoCDS;down_3monoCDS"}
              }
              
            }
          }
          
          
        }
        if(length(ref_over)==1){
          ref<-ref_over
          spl_type<-NA
          if(as.vector(strand(orf_gen[1]))=="+"){
            if(start(ran)==(start(ref_over))){
              if(end(ran)==(end(ref_over))){
                spl_type<-"same_5ss;same_3ss"
                if(f==last_ex){spl_type<-"same_5ss;same_lastCDS"}
                if(f==1){spl_type<-"same_firstCDS;same_3ss"}
                if(1==last_ex){spl_type<-"same_5monoCDS;same_3monoCDS"}
              }
              if(end(ran)>(end(ref_over))){
                spl_type<-"same_5ss;down_3ss"
                if(f==last_ex){spl_type<-"same_5ss;down_lastCDS"}
                if(f==1){spl_type<-"same_firstCDS;down_3ss"}
                if(1==last_ex){spl_type<-"same_5monoCDS;down_3monoCDS"}
              }
              if(end(ran)<(end(ref_over))){
                spl_type<-"same_5ss;up_3ss"
                if(f==last_ex){spl_type<-"same_5ss;up_lastCDS"}
                if(f==1){spl_type<-"same_firstCDS;up_3ss"}
                if(1==last_ex){spl_type<-"same_5monoCDS;up_3monoCDS"}
              }
              
            }
            if(end(ran)==(end(ref_over))){
              if(start(ran)>(start(ref_over))){
                spl_type<-"down_5ss;same_3ss"
                if(f==last_ex){spl_type<-"down_5ss;same_lastCDS"}
                if(f==1){spl_type<-"down_firstCDS;same_3ss"}
                if(1==last_ex){spl_type<-"down_5monoCDS;same_3monoCDS"}
              }
              if(start(ran)<(start(ref_over))){
                spl_type<-"up_5ss;same_3ss"
                if(f==last_ex){spl_type<-"up_5ss;same_lastCDS"}
                if(f==1){spl_type<-"up_firstCDS;same_3ss"}
                if(1==last_ex){spl_type<-"up_5monoCDS;same_3monoCDS"}
              }
              
            }
            
            if(end(ran)>(end(ref_over))){
              if(start(ran)>(start(ref_over))){
                spl_type<-"down_5ss;down_3ss"
                if(f==last_ex){spl_type<-"down_5ss;down_lastCDS"}
                if(f==1){spl_type<-"down_firstCDS;down_3ss"}
                if(1==last_ex){spl_type<-"down_5monoCDS;down_3monoCDS"}
              }
              if(start(ran)<(start(ref_over))){
                spl_type<-"up_5ss;down_3ss"
                if(f==last_ex){spl_type<-"up_5ss;down_lastCDS"}
                if(f==1){spl_type<-"up_firstCDS;down_3ss"}
                if(1==last_ex){spl_type<-"up_5monoCDS;down_3monoCDS"}
              }
              
            }
            
            if(end(ran)<(end(ref_over))){
              if(start(ran)>(start(ref_over))){
                spl_type<-"down_5ss;up_3ss"
                if(f==last_ex){spl_type<-"down_5ss;up_lastCDS"}
                if(f==1){spl_type<-"down_firstCDS;up_3ss"}
                if(1==last_ex){spl_type<-"down_5monoCDS;up_3monoCDS"}
              }
              if(start(ran)<(start(ref_over))){
                spl_type<-"up_5ss;up_3ss"
                if(f==last_ex){spl_type<-"up_5ss;up_lastCDS"}
                if(f==1){spl_type<-"up_firstCDS;up_3ss"}
                if(1==last_ex){spl_type<-"up_5monoCDS;up_3monoCDS"}
              }
              
            }
            
            
            
            
          }
          
          if(as.vector(strand(orf_gen[1]))=="-"){
            if(start(ran)==(start(ref_over))){
              if(end(ran)==(end(ref_over))){
                spl_type<-"same_5ss;same_3ss"
                if(f==1){spl_type<-"same_5ss;same_lastCDS"}
                if(f==last_ex){spl_type<-"same_firstCDS;same_3ss"}
                if(1==last_ex){spl_type<-"same_5monoCDS;same_3monoCDS"}
              }
              if(end(ran)>(end(ref_over))){
                spl_type<-"up_5ss;same_3ss"
                if(f==1){spl_type<-"up_5ss;same_lastCDS"}
                if(f==last_ex){spl_type<-"up_firstCDS;same_3ss"}
                if(1==last_ex){spl_type<-"up_5monoCDS;same_3monoCDS"}
              }
              if(end(ran)<(end(ref_over))){
                spl_type<-"down_5ss;same_3ss"
                if(f==1){spl_type<-"down_5ss;same_lastCDS"}
                if(f==last_ex){spl_type<-"down_firstCDS;same_3ss"}
                if(1==last_ex){spl_type<-"down_5monoCDS;same_3monoCDS"}
              }
              
            }
            if(end(ran)==(end(ref_over))){
              if(start(ran)>(start(ref_over))){
                spl_type<-"same_5ss;up_3ss"
                if(f==1){spl_type<-"same_5ss;up_lastCDS"}
                if(f==last_ex){spl_type<-"same_firstCDS;up_3ss"}
                if(1==last_ex){spl_type<-"same_5monoCDS;up_3monoCDS"}
              }
              if(start(ran)<(start(ref_over))){
                spl_type<-"same_5ss;down_3ss"
                if(f==1){spl_type<-"same_5ss;down_lastCDS"}
                if(f==last_ex){spl_type<-"same_firstCDS;down_3ss"}
                if(1==last_ex){spl_type<-"same_5monoCDS;down_3monoCDS"}
              }
              
            }
            
            if(end(ran)>(end(ref_over))){
              if(start(ran)>(start(ref_over))){
                spl_type<-"up_5ss;up_3ss"
                if(f==1){spl_type<-"up_5ss;up_lastCDS"}
                if(f==last_ex){spl_type<-"up_firstCDS;up_3ss"}
                if(1==last_ex){spl_type<-"up_5monoCDS;up_3monoCDS"}
              }
              if(start(ran)<(start(ref_over))){
                spl_type<-"up_5ss;down_3ss"
                if(f==1){spl_type<-"up_5ss;down_lastCDS"}
                if(f==last_ex){spl_type<-"up_firstCDS;down_3ss"}
                if(1==last_ex){spl_type<-"up_5monoCDS;down_3monoCDS"}
              }
              
            }
            
            if(end(ran)<(end(ref_over))){
              if(start(ran)>(start(ref_over))){
                spl_type<-"down_5ss;up_3ss"
                if(f==1){spl_type<-"down_5ss;up_lastCDS"}
                if(f==last_ex){spl_type<-"down_firstCDS;up_3ss"}
                if(1==last_ex){spl_type<-"down_5monoCDS;up_3monoCDS"}
              }
              if(start(ran)<(start(ref_over))){
                spl_type<-"down_5ss;down_3ss"
                if(f==1){spl_type<-"down_5ss;down_lastCDS"}
                if(f==last_ex){spl_type<-"down_firstCDS;down_3ss"}
                if(1==last_ex){spl_type<-"down_5monoCDS;down_3monoCDS"}
              }
              
            }
          }
          
          
          
          
        }
        
      }
      if(overref[f]==F){
        if(length(ref_cds)>0){ref<-ref_cds[nearest(x=ran,subject=ref_cds)]}
        if(length(ref_cds)==0){ref<-GRanges()}
        spl_type<-"new_CDS"
        if(f==last_ex){
          if(as.vector(strand(orf_gen[1]))=="+"){
            spl_type<-"new_lastCDS"
          }
          if(as.vector(strand(orf_gen[1]))=="-"){
            spl_type<-"new_firstCDS"
          }
        }
        if(f==1){
          
          if(as.vector(strand(orf_gen[1]))=="-"){
            spl_type<-"new_lastCDS"
          }
          if(as.vector(strand(orf_gen[1]))=="+"){
            spl_type<-"new_firstCDS"
          }                                       
        }
        if(1==last_ex){
          spl_type<-"new_monoCDS"
          
        }
      }
      refs[f]<-list(ref)
      spl_types[f]<-list(spl_type)
      
      
    }
    if(!any(vapply(refs,is.null,NA)) & !any(vapply(spl_types,is.null,NA))){
      ran<-orf_gen
      cols<-mcols(ran)
      cols$ref<-GRangesList_fast(refs)
      cols$spl_type<-unlist(spl_types)
      mcols(ran)<-cols
      spl_ran<-sort(c(spl_ran,ran))
    }else{
      #an exon without ref and spl_type gets no such columns: exons combined one by one, as before
      rans<-lapply(1:length(orf_gen),function(f){
        ran<-orf_gen[f]
        if(!is.null(refs[[f]])){ran$ref<-GRangesList(refs[[f]])}
        ran$spl_type<-spl_types[[f]]
        ran
      })
      spl_ran<-sort(do.call(c,c(list(spl_ran),rans)))
    }
    
  }
  
  newspl<-spl_ran$spl_type
  newspl[grep(spl_ran$spl_type,pattern = "new")]<-"new_miss"
  newspl[grep(spl_ran$spl_type,pattern = "missing")]<-"new_miss"
  rlesp<-Rle(newspl)
  
  #change new and missing
  cols<-mcols(spl_ran)
  if(runValue(rlesp)[1]=="new_miss"){
    if(as.character(strand(spl_ran)[1])=="+"){cols$spl_type[1:runLength(rlesp)[1]]<-paste(cols$spl_type[1:runLength(rlesp)[1]],"_5prime",sep = "")}
    if(as.character(strand(spl_ran)[1])=="-"){cols$spl_type[1:runLength(rlesp)[1]]<-paste(cols$spl_type[1:runLength(rlesp)[1]],"_3prime",sep = "")}
  }
  lenna<-length(runValue(rlesp))
  if(runValue(rlesp)[lenna]=="new_miss"){
    if(as.character(strand(spl_ran)[1])=="+"){cols$spl_type[(length(spl_ran)-(runLength(rlesp)[lenna]-1)):length(spl_ran)]<-paste(cols$spl_type[(length(spl_ran)-(runLength(rlesp)[lenna]-1)):length(spl_ran)],"_3prime",sep = "")}
    if(as.character(strand(spl_ran)[1])=="-"){paste(cols$spl_type[(length(spl_ran)-(runLength(rlesp)[lenna]-1)):length(spl_ran)],"_5prime",sep = "")}
  }
  cols$spl_type<-gsub(cols$spl_type,pattern = "_5prime_3prime",replacement = "_notoverl")
  mcols(spl_ran)<-cols
  spl_ran
}

# ORF_category_Tx of an ORF (start sta and end sto in the transcript, without the stop codon) with
# respect to the annotated CDS (start ann_sta, and ann_sto = CDS end - 3). It replaces a cascade of
# if where the last true one won (uORF over overl_uORF, dORF over overl_dORF), and gives the same
# label for all inputs, also for a CDS shorter than 4 nt (ann_sta>ann_sto)
tx_category<-function(sta,sto,ann_sta,ann_sto){
  if(sto==ann_sto){
    if(sta==ann_sta){"ORF_annotated"} else if(sta<ann_sta){"N_extension"} else {"N_truncation"}
  } else if(sta==ann_sta){
    if(sto<ann_sto){"C_truncation"} else {"C_extension"}
  } else if(sto<ann_sto){
    if(sta>ann_sta){"nested_ORF"} else if(sto<ann_sta){"uORF"} else {"overl_uORF"}
  } else {
    if(sta>ann_sto){"dORF"} else if(sta<ann_sta){"NC_extension"} else {"overl_dORF"}
  }
}

#' Annotate detected ORFs in transcript and genome space
#'
#' This function annotates quantified ORFs with respect to other detected ORFs and annotated ones, in both genome and transcript space.
#' @details As multiple transcripts can contain the same ORF,
#' all the transcript and transcript biotypes are indicated, with a preference for protein_coding transcripts in the "compatible" 
#' columns (to be conservative when assessing translation of non-protein coding transcripts). Such compatibility is also output considering the most upstream
#' start codon for that ORF. \cr
#' Splice features of each orf is annotated with respect to the longest coding transcripts and to the highest translated ORF in that gene.\cr
#' Variants in N or C terminus of the translated proteins are also indicated (Beta).\cr
#' ORF annotation with respect to the annotated transcript is also indicated, as follows:\cr\cr
#' \code{novel}: no ORF annotated in the transcript.\cr
#' \code{ORF_annotated}: same exact ORF as annotated.\cr
#' \code{N_extension}: N terminal extension.\cr
#' \code{N_truncation}: N terminal truncation.\cr
#' \code{uORF}: upstream ORF.\cr
#' \code{overl_uORF}: upstream overlappin uORF.\cr
#' \code{NC_extension}: N and C termini extension.\cr
#' \code{dORF}: downstream ORF.\cr
#' \code{overl_dORF}: downstream overlapping ORF.\cr
#' \code{nested_ORF}: nested ORF.\cr
#' \code{C_truncation}: C terminal truncation.\cr
#' \code{C_extension}: C terminal extension.\cr\cr
#' As transcipt-specific annotation can be misleading due to a plethora of different transcripts, it is important to distinguish ORFs
#' also on the basis of their overlap with know CDS regions.
#' ORFs that overlap no CDS exon of the analyzed genomic region are annotated in genome space as follows:\cr\cr
#' \code{novel}: No CDS region is annotated in the entire region.\cr
#' \code{novel_Upstream}: ORF is upstream of annotated CDS regions (does not overlap).\cr
#' \code{novel_Downstream}: ORF is downstream of annotated CDS regions (does not overlap).\cr
#' \code{novel_Internal}: genomic location of the ORF is present between the start of the first,
#'  and the end of the last CDS region (does not overlap).\cr\cr
#' The start and stop codons of the other ORFs are compared with those of the CDS of \code{ref_id} (the longest annotated CDS of the gene; see \code{\link{ORFquant_output}} for genes without one):\cr\cr
#' \code{exact_start_stop}: Same start and end locations.\cr
#' \code{Alt5_start}: Different start region, upstream.\cr
#' \code{Alt3_start}: Different start region, downstream.\cr
#' \code{Alt5_stop}: Different end region, upstream.\cr
#' \code{Alt3_stop}: Different end region, downstream.\cr
#' \code{Alt5_start_Alt5_stop}, \code{Alt5_start_Alt3_stop}, \code{Alt3_start_Alt5_stop}, \code{Alt3_start_Alt3_stop}: Different start and end regions.\cr\cr
#' Another layer of annotation is performed by checking the position of the ORF stop codon 
#' with respect to the last exon-exon junction.
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param results_ORFs Full list of detected ORFs, from \code{select_quantify_ORFs}
#' @param Annotation Rannot object containing annotation of CDS and transcript structures (see \code{prepare_annotation_files}
#' @param genome_sequence BSgenome object
#' @param region genomic region being analyzed
#' @param genetic_code GENETIC_CODE table to use
#' @return \code{results_ORFs} with a new element \code{ORFs_splice_feats}: a list of two \code{GRangesList}, \code{annotation_wrt_longest} and \code{annotation_wrt_maxORF},
#' with, for each ORF, its exon structure compared with the CDS of \code{ref_id} or with the ORF \code{ref_id_maxORF} (see \code{annotate_splicing}):
#' the exons, including possible missing exons from reference, with a \code{spl_type} column
#' including the annotation for each exon (e.g. alternative acceptors or donor).\cr\cr
#' Additional columns are added to the ORFs_tx object:\cr
#' \code{compatible_with}: Set of ORF ids (\code{<transcript_id>_<start>_<end>}), one on each transcript containing the entire ORF structure, including \code{ORF_id_tr}.\cr
#' \code{compatible_biotype}: Compatible transcript biotype; if a protein coding transcript can contain 
#' the ORF, this is set to protein_coding.\cr
#' \code{compatible_tx}: One selected compatible transcript (preference if protein_coding).\cr
#' \code{compatible_ORF_id_tr}: ORF_id_tr id if selecting the compatible transcript.\cr
#' \code{compatible_with_longest}: Same as \code{compatible_with} but using the most upstream start codon.\cr
#' \code{compatible_ORF_id_tr_longest}: Same as \code{compatible_ORF_id_tr} but using the most upstream start codon .\cr
#' \code{ref_id}: transcript_id of the transcript used to annotate splicing (longest) .\cr
#' \code{ref_id_maxORF}: ORF_id_tr of the ORF used to annotated splicing (most translated of the gene).\cr
#' \code{NC_protein_isoform}: Annotation of possible N or C termini variant compared with the CDS of \code{ref_id} (when the ORF overlaps an annotated CDS; \code{NA} otherwise).\cr
#' \code{ORF_category_Tx}: ORF annotation with respect to ORF position in the transcript .\cr
#' \code{ORF_category_Tx_compatible}: ORF annotation with respect to ORF position in the transcript, using the \code{compatible_ORF_id_tr} .\cr
#' \code{ORF_category_Gen}: ORF annotation with respect to its genomic position .\cr
#' \code{NMD_candidate}: TRUE or FALSE, depending on the presence of an additional exon-exon junction downstream the stop codon.\cr
#' \code{NMD_candidate_compatible_txs}: same as NMD_candidate, but for all transcripts compatible with the ORF structure.\cr
#' \code{Distance_to_lastExEx}: Distance (in nt) between the last exon-exon junction and the stop codon.\cr
#' \code{Distance_to_lastExEx_compatible_txs}: same as  Distance_to_lastExEx, but for all transcripts compatible with the ORF structure.
#' @seealso \code{\link{select_quantify_ORFs}}, \code{\link{annotate_splicing}}, \code{\link{ORFquant_output}}
#' @export

annotate_ORFs<-function(results_ORFs,Annotation,genome_sequence,region,genetic_code){
  
  annotated_cds<-Annotation$cds_genes[Annotation$cds_genes%over%region]
  annotated_cds_tx<-Annotation$cds_txs[Annotation$cds_txs%over%region]
  annotated_cds_tx_genes<-Annotation$trann$gene_id[match(names(annotated_cds_tx),Annotation$trann$transcript_id)]
  annotated_exons_tx<-Annotation$exons_txs[Annotation$exons_txs%over%region]
  
  ORFs_tx<-results_ORFs$ORFs_tx_position
  ORFs_gen<-results_ORFs$ORFs_genomic_position
  
  orfssss_tx<-unlist(GRangesList(ORFs_tx))
  orfssss_tx<-orfssss_tx[!is.na(orfssss_tx$ORF_pct_P_sites)]
  maxORF_orf<-c()
  max_cds<-c()
  max_cdsok<-GRanges()
  
  for(i in unique(orfssss_tx$gene_id)){
    okorfsss<-orfssss_tx[orfssss_tx$gene_id==i]
    maxORF_orf[i]<-okorfsss$ORF_id_tr[which.max(okorfsss$ORF_pct_P_sites)]
    anncdss<-annotated_cds_tx[annotated_cds_tx_genes==i]
    if(length(anncdss)>0){
      max_cds[i]<-names(which.max(sum(width(anncdss))))
    }
    
    if(length(anncdss)==0){
      orf_gen<-reduce(unlist(ORFs_gen[okorfsss$ORF_id_tr]))
      
      moret<-sapply(annotated_cds_tx,FUN=function(x){
        sum(width(setdiff(orf_gen,x)))
      })
      moret<-t(data.frame(moret,stringsAsFactors=F))
      max_cds[i]<-colnames(moret)[which.min(colSums(moret))]
    }
    
    
  }
  
  #compatibilities
  
  for(i in names(ORFs_gen)){
    x<-ORFs_gen[[i]]
    comp<-c()
    # one mapToTranscripts for all transcripts, then split: a call per transcript costs ~50 ms
    mapp<-mapToTranscripts(x,transcripts = annotated_exons_tx)
    mapp<-reduce(split(mapp,seqnames(mapp)))
    mapp<-unlist(mapp[sum(width(mapp))==sum(width(x)) & elementNROWS(mapp)==1])
    if(length(mapp)>0){
      comp<-paste(seqnames(mapp),start(mapp),end(mapp),sep="_")
    }
    names(comp)<-NULL
    # columns are set on mcols() and put back once per ORF: a GRanges $<- also runs updateObject() (~20 ms)
    cols<-mcols(ORFs_tx[[i]])
    cols$compatible_with<-NULL
    cols$compatible_with<-CharacterList(unlist(comp))
    
    
    cols$compatible_biotype<-cols$transcript_biotype
    cols$compatible_tx<-cols$transcript_id
    cols$compatible_ORF_id_tr<-cols$ORF_id_tr
    compats<-elementNROWS(cols$compatible_with)>1
    comp_txs<-cols$compatible_with[compats]
    if(length(comp_txs)>0){
      compid<-sapply(comp_txs,function(x){
        txs<-sapply(strsplit(x,split = "_"),function(x){paste(x[-((length(x)-1):length(x))],collapse="_")})
        btps<-Annotation$trann$transcript_biotype[match(txs,Annotation$trann$transcript_id)]
        btps[is.na(btps)]<-"Not_found"
        pcd<-btps=="protein_coding"
        if(sum(pcd,na.rm = T)>0){c(sort(x[pcd])[1],sort(txs[pcd])[1],"protein_coding")}else{c(x[1],txs[1],btps[1])}
      })
      cols$compatible_ORF_id_tr[compats]<-t(compid)[,1]
      cols$compatible_tx[compats]<-t(compid)[,2]
      cols$compatible_biotype[compats]<-t(compid)[,3]
      
    }
    cols$compatible_with_longest<-cols$compatible_with
    cols$compatible_biotype_longest<-cols$compatible_biotype
    cols$compatible_tx_longest<-cols$compatible_tx
    cols$compatible_ORF_id_tr_longest<-cols$compatible_ORF_id_tr
    lng<-GRanges(cols$longest_ORF)
    mcols(lng)$ORF_id_tr<-paste(seqnames(lng),start(lng),end(lng),sep="_")
    if(start(ORFs_tx[[i]])!=start(lng)){
      x<-from_tx_togen(ORFs = lng,exons = Annotation$exons_txs[cols$transcript_id],introns = Annotation$introns_txs[[cols$transcript_id]])[[1]]
      mapp<-mapToTranscripts(x,transcripts = annotated_exons_tx)
      redmapp<-reduce(split(mapp,seqnames(mapp)))
      
      redmapp<-redmapp[which(sum(width(redmapp))==sum(width(x)))]
      redmapp<-redmapp[elementNROWS(redmapp)==1]
      
      comp<-sapply(redmapp,function(x){paste(seqnames(x),start(x),end(x),sep="_")})
      names(comp)<-NULL
      cols$compatible_with_longest<-CharacterList(unlist(comp))
      cols$compatible_biotype_longest<-cols$transcript_biotype
      cols$compatible_tx_longest<-cols$transcript_id
      cols$compatible_ORF_id_tr_longest<-paste(as.character(seqnames(lng)[1]),start(lng),end(lng),sep = "_")
      
      compats<-elementNROWS(cols$compatible_with)>1
      comp_txs<-cols$compatible_with_longest[compats]
      if(length(comp_txs)>0){
        compid<-sapply(comp_txs,function(x){
          txs<-sapply(strsplit(x,split = "_"),function(x){paste(x[-((length(x)-1):length(x))],collapse="_")})
          btps<-Annotation$trann$transcript_biotype[match(txs,Annotation$trann$transcript_id)]
          btps[is.na(btps)]<-"Not_found"
          pcd<-btps=="protein_coding"
          if(sum(pcd,na.rm = T)>0){c(sort(x[pcd])[1],sort(txs[pcd])[1],"protein_coding")}else{c(x[1],txs[1],btps[1])}
        })
        cols$compatible_ORF_id_tr_longest[compats]<-t(compid)[,1]
        cols$compatible_tx_longest[compats]<-t(compid)[,2]
        cols$compatible_biotype_longest[compats]<-t(compid)[,3]
        
      }        
      
    }
    mcols(ORFs_tx[[i]])<-cols
    
  }
  
  for(i in names(ORFs_gen)){
    
    cols<-mcols(ORFs_tx[[i]])
    compss<-cols$compatible_with[[1]]
    compss_txs<-sapply(strsplit(compss,split = "_"),function(x){paste(x[-((length(x)-1):length(x))],collapse="_")})
    ok_id<-compss[compss_txs==cols$compatible_tx_longest]
    comp_ln<-cols$compatible_biotype_longest
    comp_prev<-cols$compatible_biotype
    if(comp_ln%in%"protein_coding" | (!comp_prev%in%"protein_coding" & !comp_ln%in%"protein_coding") ){
      cols$compatible_ORF_id_tr<-ok_id
      cols$compatible_tx<-cols$compatible_tx_longest
      cols$compatible_biotype<-cols$compatible_biotype_longest
    }
    cols$compatible_tx_longest<-NULL
    cols$compatible_biotype_longest<-NULL
    mcols(ORFs_tx[[i]])<-cols
    
  }
  
  #annotate NMD candidates
  exs<-annotated_exons_tx[unique(unlist(GRangesList(ORFs_tx))$transcript_id)]
  
  strands_exs<-sapply(strand(exs),function(x){x@values[1]})
  exs_pos<-exs[strands_exs=="+"]
  exs_neg<-exs[strands_exs=="-"]
  
  
  last_ex_pos<-which.max(start(exs_pos))
  last_ex_pos<-exs_pos[splitAsList(unname(last_ex_pos), names(last_ex_pos))]
  last_ex_pos<-last_ex_pos[match(names(exs_pos),names(last_ex_pos))]
  txs_pos<-pmapToTranscripts(last_ex_pos,transcripts = exs_pos)
  
  
  last_ex_neg<-which.min(start(exs_neg))
  last_ex_neg<-exs_neg[splitAsList(unname(last_ex_neg), names(last_ex_neg))]
  last_ex_neg<-last_ex_neg[match(names(exs_neg),names(last_ex_neg))]
  txs_neg<-pmapToTranscripts(last_ex_neg,transcripts = exs_neg)
  txs_all<-unlist(c(txs_pos,txs_neg))
  last_exexs<-start(txs_all)[match(unlist(GRangesList(ORFs_tx))$transcript_id,names(txs_all))]
  Distance_EJCs<-last_exexs-end(unlist(GRangesList(ORFs_tx)))
  
  #NMD_compat
  unlORFs<-unlist(GRangesList(ORFs_tx))
  unlORFs_comp<-sapply(strsplit(unlist(unlORFs$compatible_with),"_"),function(x){paste(x[-((length(x)-1):length(x))],collapse="_")})
  exs<-annotated_exons_tx[unlORFs_comp]
  
  strands_exs<-sapply(strand(exs),function(x){x@values[1]})
  exs_pos<-exs[strands_exs=="+"]
  exs_neg<-exs[strands_exs=="-"]
  
  last_ex_pos<-which.max(start(exs_pos))
  last_ex_pos<-exs_pos[splitAsList(unname(last_ex_pos), names(last_ex_pos))]
  last_ex_pos<-last_ex_pos[match(names(exs_pos),names(last_ex_pos))]
  txs_pos<-pmapToTranscripts(last_ex_pos,transcripts = exs_pos)
  
  
  last_ex_neg<-which.min(start(exs_neg))
  last_ex_neg<-exs_neg[splitAsList(unname(last_ex_neg), names(last_ex_neg))]
  last_ex_neg<-last_ex_neg[match(names(exs_neg),names(last_ex_neg))]
  txs_neg<-pmapToTranscripts(last_ex_neg,transcripts = exs_neg)
  txs_all<-unlist(c(txs_pos,txs_neg)[unlORFs_comp])
  
  last_exexs<-start(txs_all)
  
  endcompat<-as.numeric(unlist(lapply(strsplit(unlist(unlORFs$compatible_with),"_"),function(x){x[length(x)]})))
  
  Distance_EJCs_compat<-last_exexs-endcompat
  nrws<-elementNROWS(unlORFs$compatible_with)
  fcttr<-rep(1:length(nrws),nrws)
  Distance_EJCs_compat<-CharacterList(unname(split(Distance_EJCs_compat,fcttr)))
  
  ORFs_tx<-lapply(ORFs_tx,function(x){
    cols<-mcols(x)
    cols[,c("ref_id","ref_id_maxORF","NC_protein_isoform","ORF_category_Tx","ORF_category_Tx_compatible","ORF_category_Gen",
            "NMD_candidate","Distance_to_lastExEx")]<-NA
    cols[,c("NMD_candidate_compatible_txs","Distance_to_lastExEx_compatible_txs")]<-CharacterList("")
    mcols(x)<-cols
    x
  })
  
  for(i in 1:length(ORFs_tx)){
    cols<-mcols(ORFs_tx[[i]])
    nmd<-FALSE
    if(Distance_EJCs[i]>0){
      nmd<-TRUE
    }
    cols$NMD_candidate<-nmd
    cols$Distance_to_lastExEx<-Distance_EJCs[i]
    
    nmd<-Distance_EJCs_compat[i]>0
    cols$NMD_candidate_compatible_txs<-nmd
    cols$Distance_to_lastExEx_compatible_txs<-Distance_EJCs_compat[i]
    mcols(ORFs_tx[[i]])<-cols
  }
  
  
  
  #here bulk of work
  ORFs_splice_feats<-list()
  ORFs_splice_feats_tomaxORF<-list()
  
  for(i in 1:length(ORFs_tx)){
    orf_tx<-ORFs_tx[[i]]
    cols<-mcols(ORFs_tx[[i]])
    ORFs_splice_feats[[orf_tx$ORF_id_tr]]<-GRanges()
    ORFs_splice_feats_tomaxORF[[orf_tx$ORF_id_tr]]<-GRanges()
    
    #annotate tx position
    
    annotated_ORF<-Annotation$cds_txs_coords[as.character(seqnames(Annotation$cds_txs_coords))==orf_tx$transcript_id]
    annotated_ORF_compatible<-Annotation$cds_txs_coords[as.character(seqnames(Annotation$cds_txs_coords))==orf_tx$compatible_tx]
    
    if(length(annotated_ORF)==0){
      cols$ORF_category_Tx<-"novel"
    }
    if(length(annotated_ORF_compatible)==0){
      cols$ORF_category_Tx_compatible<-"novel"
    }
    
    if(length(annotated_ORF)>0){
      
      ann_sta<-start(annotated_ORF)
      ann_sto<-end(annotated_ORF)-3
      sta<-start(orf_tx)
      sto<-end(orf_tx)
      cols$ORF_category_Tx<-tx_category(sta,sto,ann_sta,ann_sto)
    }
    
    
    if(length(annotated_ORF_compatible)>0){
      
      ann_sta<-start(annotated_ORF_compatible)
      ann_sto<-end(annotated_ORF_compatible)-3
      #change
      sta<-as.numeric(sapply(strsplit(orf_tx$compatible_ORF_id_tr,"_"),function(x){x[length(x)-1]}))
      sto<-as.numeric(sapply(strsplit(orf_tx$compatible_ORF_id_tr,"_"),function(x){x[length(x)]}))
      cols$ORF_category_Tx_compatible<-tx_category(sta,sto,ann_sta,ann_sto)
    }
    
    
    #annotate splice, genomic position and protein termini based on max cds and max pct
    
    orf_gen<-ORFs_gen[[orf_tx$ORF_id_tr]]
    
    #1 of multiple overlapping cds per gene
    
    if(length(annotated_cds)>1){
      moreg<-sapply(annotated_cds,FUN=function(x){
        orf_gen%over%x
      })
      if(length(orf_gen)==1){moreg<-t(data.frame(moreg,stringsAsFactors=F))}
      annotated_cds2<-reduce(unlist(annotated_cds[[colnames(moreg)[which.max(colSums(moreg))]]]))
      
    }
    
    #otherwise list with 1 gene
    if(length(annotated_cds)==1){annotated_cds2<-reduce(unlist(annotated_cds))}
    if(length(annotated_cds)==0){annotated_cds2<-annotated_cds}
    overl<-orf_gen%over%annotated_cds2
    
    
    if(sum(overl)==0){
      cols$ORF_category_Gen<-"novel"
      if(length(annotated_cds2)>0){
        nearest_cds<-annotated_cds[[names(unlist(annotated_cds))[nearest(orf_gen,unlist(annotated_cds))[1]]]]
        overl_whole<-orf_gen@ranges%over%IRanges(start=min(start(nearest_cds)),end=max(end(nearest_cds)))
        
        if(as.vector(strand(orf_gen[1]))=="+"){
          
          if(sum(overl_whole)==0){
            if(min(start(orf_gen))<min(start(nearest_cds))){
              cols$ORF_category_Gen<-"novel_Upstream"
            }
            if(min(start(orf_gen))>max(end(nearest_cds))){
              cols$ORF_category_Gen<-"novel_Downstream"
            }
          }
          if(sum(overl_whole)>0){
            cols$ORF_category_Gen<-"novel_Internal"
          }
          
        }
        if(as.vector(strand(orf_gen[1]))=="-"){
          if(sum(overl_whole)==0){
            if(max(end(orf_gen))>max(end(nearest_cds))){
              cols$ORF_category_Gen<-"novel_Upstream"
            }
            if(min(start(orf_gen))<min(start(nearest_cds))){
              cols$ORF_category_Gen<-"novel_Downstream"
            }
          }
          if(sum(overl_whole)>0){
            cols$ORF_category_Gen<-"novel_Internal"
          }                                        
        }
      }
    }
    
    #if overlaps CDS
    if(sum(overl)>0){
      
      #annotated NC wrt maxcds
      
      cols$ORF_category_Gen<-"overlaps_CDS"
      max_cdsok<-annotated_cds_tx[[max_cds[cols$gene_id]]]
      if(length(max_cdsok)==0){
        max_cdsok<-annotated_cds_tx[[which.max(sum(width(annotated_cds_tx)))]]
      }
      cols$ref_id<-max_cds[cols$gene_id]
      max_cds_seq<-unlist(getSeq(x=genome_sequence,max_cdsok))
      segm<-10
      if(length(max_cds_seq)<33 | nchar(orf_tx$Protein)<10){
        segm<-min(c(as.integer(length(max_cds_seq)/3)-1,nchar(orf_tx$Protein)))
      }
      N_pr<-AAString(unlist(orf_tx$Protein))[1:segm]
      C_pr<-AAString(unlist(orf_tx$Protein))[(nchar(unlist(orf_tx$Protein))-(segm-1)):nchar(unlist(orf_tx$Protein))]
      N_ann<-translate_solve(max_cds_seq[1:(segm*3)],genetic.code = genetic_code)
      C_ann<-translate_solve(head(tail(max_cds_seq,(segm*3+3)),(segm*3)),genetic.code = genetic_code,no.init.codon=segm<nchar(unlist(orf_tx$Protein)))
      
      
      cols$NC_protein_isoform<-"N_C"
      if(N_pr==N_ann){
        if(C_pr==C_ann){
          cols$NC_protein_isoform<-"same"
        } else {cols$NC_protein_isoform<-"C"}
      }
      if(C_pr==C_ann & N_pr!=N_ann){
        cols$NC_protein_isoform<-"N"                                        
      }
      
      #ann genomic
      
      if(as.vector(strand(orf_gen[1]))=="+"){
        
        gen_sta<-min(start(max_cdsok))
        gen_sto<-max(end(max_cdsok))
        
        sta_or<-min(start(orf_gen))
        sto_or<-max(end(orf_gen))+3
        if(sto_or==gen_sto){
          if(sta_or==gen_sta){cols$ORF_category_Gen<-"exact_start_stop"}
          if(sta_or<gen_sta){cols$ORF_category_Gen<-"Alt5_start"}
          if(sta_or>gen_sta){cols$ORF_category_Gen<-"Alt3_start"}
          
        }
        
        if(sto_or!=gen_sto){
          if(sta_or<gen_sta & sto_or<gen_sto){cols$ORF_category_Gen<-"Alt5_start_Alt5_stop"}
          if(sta_or<gen_sta & sto_or>gen_sto){cols$ORF_category_Gen<-"Alt5_start_Alt3_stop"}
          if(sta_or>gen_sta & sto_or>gen_sto){cols$ORF_category_Gen<-"Alt3_start_Alt3_stop"}
          if(sta_or>gen_sta & sto_or<gen_sto){cols$ORF_category_Gen<-"Alt3_start_Alt5_stop"}
          if(sta_or==gen_sta & sto_or<gen_sto){cols$ORF_category_Gen<-"Alt5_stop"}
          if(sta_or==gen_sta & sto_or>gen_sto){cols$ORF_category_Gen<-"Alt3_stop"}
          
        }
      }
      if(as.vector(strand(orf_gen[1]))=="-"){
        
        gen_sta<-max(end(max_cdsok))
        gen_sto<-min(start(max_cdsok))
        
        sta_or<-max(end(orf_gen))
        sto_or<-min(start(orf_gen))-3
        
        if(sto_or==gen_sto){
          if(sta_or==gen_sta){cols$ORF_category_Gen<-"exact_start_stop"}
          if(sta_or<gen_sta){cols$ORF_category_Gen<-"Alt3_start"}
          if(sta_or>gen_sta){cols$ORF_category_Gen<-"Alt5_start"}
          
        }
        
        if(sto_or!=gen_sto){
          if(sta_or>gen_sta & sto_or>gen_sto){cols$ORF_category_Gen<-"Alt5_start_Alt5_stop"}
          
          if(sta_or>gen_sta & sto_or<gen_sto){cols$ORF_category_Gen<-"Alt5_start_Alt3_stop"}
          if(sta_or<gen_sta & sto_or<gen_sto){cols$ORF_category_Gen<-"Alt3_start_Alt3_stop"}
          if(sta_or<gen_sta & sto_or>gen_sto){cols$ORF_category_Gen<-"Alt3_start_Alt5_stop"}
          if(sta_or==gen_sta & sto_or>gen_sto){cols$ORF_category_Gen<-"Alt5_stop"}
          if(sta_or==gen_sta & sto_or<gen_sto){cols$ORF_category_Gen<-"Alt3_stop"}
          
        }
        
        
        
      }
    }     
    
    ORFs_splice_feats[[orf_tx$ORF_id_tr]]<-annotate_splicing(orf_gen = orf_gen,ref_cds = max_cdsok)
    #to maxORF
    max_pct<-ORFs_gen[[maxORF_orf[cols$gene_id]]]
    cols$ref_id_maxORF<-maxORF_orf[cols$gene_id]
    ORFs_splice_feats_tomaxORF[[orf_tx$ORF_id_tr]]<-annotate_splicing(orf_gen = orf_gen,ref_cds = max_pct)
    mcols(ORFs_tx[[i]])<-cols
    
    
    
  }
  ORFs_splice_feats<-GRangesList(ORFs_splice_feats)
  ORFs_splice_feats_tomaxORF<-GRangesList(ORFs_splice_feats_tomaxORF)
  results_ORFs$ORFs_tx_position<-ORFs_tx
  list_spl_res<-list(ORFs_splice_feats,ORFs_splice_feats_tomaxORF)
  names(list_spl_res)<-c("annotation_wrt_longest","annotation_wrt_maxORF")
  results_ORFs$ORFs_splice_feats<-list_spl_res
  return(results_ORFs)
}
