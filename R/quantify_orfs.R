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

# The selection of the transcripts with Ribo-seq data, and the selection and
# quantification of ORFs (with the readthrough regions).


#' Select a subset of transcripts with Ribo-seq data
#'
#' This function flattens all annotated transcript structures and uses Ribo-seq to select a subset of transcripts.
#' @details Features (bins and junctions) are divided into shared and unique features, and into with support and without
#' support (with  or without reads mapping). A set of logical rules filters out transcripts with internal features with no support
#' and no unique features with reads. More specific details can be found in the ORFquant manuscript.
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param region genomic region being analyzed
#' @param annotation Rannot object containing annotation of CDS and transcript structures (see \code{prepare_annotation_files})
#' @param P_sites GRanges object with P_sites positions
#' @param P_sites_uniq GRanges object with uniquely mapping P_sites positions
#' @param junction_counts GRanges object containing Ribo-seq counts on the set of annotated junctions
#' @param uniq_signal Use only signal from uniquely mapping reads? Defaults to \code{FALSE}.
#' @return GRanges object with the set of counts on each exonic bin and junctions, together with the list
#' of selected transcripts
#' @seealso \code{\link{prepare_annotation_files}}
#' @export

select_txs<-function(region,annotation,P_sites,P_sites_uniq,junction_counts,uniq_signal=F){
  
  gene_feat <- junction_counts[junction_counts %over% region]
  
  nsns<-annotation$exons_bins
  sel_nsns<-nsns%over%region
  gen_nsns<-nsns[sel_nsns]
  genbin<-gen_nsns
  # columns are set on mcols() and put back once: a GRanges $<- also runs updateObject() (~20 ms)
  cols<-mcols(genbin)
  cols$exonic_part<-NULL
  cols$type<-"E"
  cols$reads<-0
  cols$unique_reads<-0
  mcols(genbin)<-cols
  
  d<-genbin$reads
  
  hts<-findOverlaps(genbin,P_sites,ignore.strand=F)
  hts<-cbind(queryHits(hts),P_sites[subjectHits(hts)]$score*width(P_sites[subjectHits(hts)]))
  if(length(hts)>0){
    hts<-aggregate(x = hts[,2],list(hts[,1]),FUN=sum)
    for(i in 1:dim(hts)[1]){
      d[hts[i,1]]<-hts[i,2]
    }
    mcols(genbin)$reads<-d
    
  }
  d<-genbin$unique_reads
  hts<-findOverlaps(genbin,P_sites_uniq,ignore.strand=F)
  hts<-cbind(queryHits(hts),P_sites_uniq[subjectHits(hts)]$score*width(P_sites_uniq[subjectHits(hts)]))
  #IMPORTANT, HERE THERE WAS A IF NO UNIQ RETURN GRANGESLIST()
  if(length(hts)>0){
    hts<-aggregate(x = hts[,2],list(hts[,1]),FUN=sum)
    for(i in 1:dim(hts)[1]){
      d[hts[i,1]]<-hts[i,2]
    }
    mcols(genbin)$unique_reads<-d
    
  }
  mcols(genbin)<-mcols(genbin)[,c("tx_name","gene_id","type","reads","unique_reads")]
  mcols(gene_feat)<-mcols(gene_feat)[,names(mcols(genbin))]
  
  gene_feat<-sort(c(gene_feat,genbin))
  cols<-mcols(gene_feat)
  cols$txs<-cols$tx_name
  cols$tx_name<-NULL
  mcols(gene_feat)<-cols
  
  
  a<-gene_feat$txs
  b<-gene_feat$gene_id
  c<-gene_feat$type
  
  #HERE uniq_flag
  if(uniq_signal){
    d<-gene_feat$unique_reads
  }
  if(!uniq_signal){
    d<-gene_feat$reads
  }
  
  
  if(sum(d)<4){
    return(GRangesList())
  }
  
  
  txs_gene<-unique(unlist(a))     
  
  
  gen_bins_junct<-gene_feat
  if(length(txs_gene)<2){
    cols<-mcols(gen_bins_junct)
    cols$genes<-b
    cols$txs<-a
    cols$genes_selected<-b
    cols$txs_selected<-a
    cols$use<-"unique"
    mcols(gen_bins_junct)<-cols
    final_ranges<-sort(gen_bins_junct)
    return(final_ranges)
  }
  
  
  
  #first round
  txs_sofar<-txs_gene
  
  mat<-matrix(data=0,nrow=length(d),ncol=length(txs_sofar))
  colnames(mat)<-txs_sofar
  for(i in 1:length(txs_sofar)){
    mat[,i]<-sapply(a,function(x){sum(x==txs_sofar[i])})
  }
  
  nest<-c()
  ident<-c()
  for(i in 1:dim(mat)[2]){
    yes<-which(mat[,i]==1)
    nesti<-c()
    for(j in (1:dim(mat)[2])[-i]){
      yesj<-which(mat[,j]==1)
      if(identical(yes,yesj)){ident<-c(ident,paste(sort(colnames(mat)[c(i,j)]),collapse=";"))}
      #added if length> otherwise txs with same structure are both deleted
      nesti<-c(nesti,sum(yesj%in%yes)==length(yes) & length(mat[,j])>length(yes))
    }
    
    if(sum(nesti)>0){nest<-c(nest,colnames(mat)[i])}
    
  }
  if(length(ident)>0){nest<-nest[!nest%in%unique(sapply(strsplit(ident,";"),"[[",1))]}
  txs_sofar<-txs_sofar[!txs_sofar%in%nest]
  
  change<-1
  
  while(change>0){
    
    mat<-matrix(data=0,nrow=length(d),ncol=length(txs_sofar))
    colnames(mat)<-txs_sofar
    for(i in 1:length(txs_sofar)){
      mat[,i]<-sapply(a,function(x){sum(x==txs_sofar[i])})
    }
    
    d_count<-paste(1:length(d),d,sep="_")
    good<-d_count[which(d>0)]
    bad<-d_count[which(d==0)]
    
    txs_good<-c()
    expl_good<-c()
    for(i in 1:dim(mat)[2]){
      expl_good_old<-expl_good
      #if new good feature
      tx<-d_count[which(mat[,i]>0)]
      tx_good<-tx[which(tx%in%good)]
      tx_bad<-tx[which(tx%in%bad)]
      
      if(length(tx_good)==0){next}
      
      if(sum(!tx_good%in%expl_good_old)>0){
        
        tx_torem<-c()
        tx_toscreen<-which(colnames(mat)%in%txs_good)
        if(length(tx_toscreen)>0){
          for(j in tx_toscreen){
            tx_contr<-d_count[which(mat[,j]>0)]
            tx_contr_good<-tx_contr[which(tx_contr%in%good)]
            tx_contr_bad<-tx_contr[which(tx_contr%in%bad)]
            if(sum(tx_contr_good%in%tx_good)==length(tx_contr_good)){
              tx_torem<-c(tx_torem,colnames(mat)[j])
              
            }
            
            if(length(tx_torem)>0){
              txs_good<-txs_good[!txs_good%in%tx_torem]
              
            }
          }
          
          
          
        }
        txs_good<-unique(c(txs_good,colnames(mat)[i]))
        expl_good<-unique(c(expl_good,tx_good))
      }
      #if same good feature, but fewer bad INTERNAL features than others.
      if(sum(!tx_good%in%expl_good_old)==0){
        tx_torem<-c()
        tx_toscreen<-which(colnames(mat)%in%txs_good)
        #here a counter when at least one good feature more than competing
        moref<-c()
        for(j in tx_toscreen){
          tx_contr<-d_count[which(mat[,j]>0)]
          tx_contr_good<-tx_contr[which(tx_contr%in%good)]
          tx_contr_bad<-tx_contr[which(tx_contr%in%bad)]
          moref<-c(moref,sum(!tx_good%in%tx_contr_good)>0)
          if(sum(tx_contr_good%in%tx_good)==length(tx_contr_good)){
            
            if(length(tx_good)>length(tx_contr_good)){
              tx_torem<-c(tx_torem,colnames(mat)[j])
            }
            #here
            if(length(tx_good)==length(tx_contr_good)){      
              fi<-which(tx==tx_good[1])
              la<-which(tx==tx_good[length(tx_good)])
              int_tx<-tx[fi:la]
              int_tx_bad<-tx_bad[tx_bad%in%int_tx]
              
              contr_fi<-which(tx_contr==tx_contr_good[1])
              contr_la<-which(tx_contr==tx_contr_good[length(tx_contr_good)])
              contr_int_tx<-tx_contr[contr_fi:contr_la]
              contr_int_tx_bad<-tx_contr_bad[tx_contr_bad%in%contr_int_tx]
              #SAME INTERNAL, TAKE
              if(length(int_tx_bad)<=length(contr_int_tx_bad)){
                
                txs_good<-unique(c(txs_good,colnames(mat)[i]))
                expl_good<-unique(c(expl_good,tx_good))
                #LESS INTERNAL, TAKE AND REMOVE OTHER
                if(length(int_tx_bad)<length(contr_int_tx_bad)){
                  
                  tx_torem<-c(tx_torem,colnames(mat)[j])}
              }
              
              
            }
          }
          
          
        }
        if(length(tx_torem)>0){
          txs_good<-txs_good[!txs_good%in%tx_torem]
          txs_good<-unique(c(txs_good,colnames(mat)[i]))
          expl_good<-unique(c(expl_good,tx_good))
        }
        
        if(length(tx_torem)==0 & sum(moref)==length(tx_toscreen)){
          txs_good<-unique(c(txs_good,colnames(mat)[i]))
          expl_good<-unique(c(expl_good,tx_good))
        }
        
        
      }
      
      
    }
    txs_sofar<-txs_good
    change<-abs(length(txs_good)-length(txs_sofar))
  }
  
  cols<-mcols(gen_bins_junct)
  cols$genes<-b
  cols$txs<-a
  a<-lapply(a,function(x){x[x%in%txs_sofar]})
  
  genes_sofar<-unique(subset(annotation$trann,annotation$trann$transcript_id%in%txs_sofar)$gene_id)
  
  b<-lapply(b,function(x){x[x%in%genes_sofar]})
  cols$genes_selected<-CharacterList(b)
  cols$txs_selected<-CharacterList(a)
  
  check<-sapply(cols$txs,FUN = length)
  use<-rep("shared",length(check))
  use[check==1]<-"unique"
  use[check==0]<-"absent"
  use[check>1]<-"shared"
  
  cols$use<-use
  mcols(gen_bins_junct)<-cols
  final_ranges<-sort(gen_bins_junct)
  return(final_ranges)
}


#' Extract possible readthrough sequences (beta) 
#'
#' This function extracts readthrough regions for subsequent analysis
#' @details The function looks for stop-stop pairs after the stop codon of the detected ORF
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param tx_name transcript_id
#' @param sequence DNAString object containing the sequence of the transcript
#' @param orf transcript-level ORF coordinates 
#' @param genetic_code GENETIC_CODE table to use
#' @return GRanges object with the set of possible readthrough sequences
#' @seealso \code{\link{detect_translated_orfs}}, \code{\link{select_quantify_ORFs}}
#' @export

get_reathr_seq<-function(tx_name,orf,sequence,genetic_code){
  length<-nchar(sequence)
  u=start(orf)%%3-1
  if(u==-1){u=2}
  pept<-NA
  pept<-unlist(strsplit(as.character(suppressWarnings(translate_solve(subseq(sequence,start=u+1),genetic.code = genetic_code))),split=""))
  
  
  starts<-pept=="M"
  
  stops<-pept=="*"
  
  start_pos<-((1:length(pept))[starts])*3
  if(length(start_pos)>0){
    start_pos<-start_pos+u-2
  } else {start_pos<-NA}
  
  stop_pos<-((1:length(pept))[stops])*3-3
  if(length(stop_pos)>0){
    #here was -2 at the end, I'll -1, as the GRanges puts 1 more ???
    stop_pos<-stop_pos+u
  } else {stop_pos<-NA}
  stop_pos<-stop_pos[stop_pos>=end(orf)]
  
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
      
      vector_starts[stoppo+1]<-start
      vector_stops[stoppo+1]<-end
      
    }
    
    
  }
  
  st_st<-data.frame(cbind(vector_starts,vector_stops))
  st_st<-st_st[!is.na(st_st[,"vector_starts"]),]
  st_st<-st_st[!is.na(st_st[,"vector_stops"]),]
  #st_st<-st_st[(st_st[,2]-st_st[,1]>11),]
  gra_par<-GRanges(seqnames = paste(tx_name,sep = "_"),strand = "+",ranges = IRanges(start = st_st[,1],end = st_st[,2]))
  gra_par<-gra_par[end(gra_par)>=end(orf)]
  return(gra_par)
}

#' Analyzed translation on possible readthrough regions (beta)
#'
#' This function uses the multitaper method to look for readthrough translation
#' @details The function looks for stop-stop pairs after the stop codon of the detected ORF
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param results_orf Full list of detected ORFs, from \code{select_quantify_ORFs} and \code{annotate_ORFs}
#' @param genome_sequence BSgenome object
#' @param annotation Rannot object containing annotation of CDS and transcript structures (see \code{prepare_annotation_files})
#' @param P_sites GRanges object with P_sites positions
#' @param P_sites_uniq GRanges object with uniquely mapping P_sites positions
#' @param P_sites_uniq_mm Rle signal of uniquely mapping P_sites with mismatches along the transcript
#' @param cutoff_fr_ave \code{cutoff} parameter for the  \code{calc_orf_pval} functions
#' @param genetic_code_table GENETIC_CODE table to use
#' @param uniq_signal Use only signal from uniquely mapping reads? Defaults to \code{FALSE}.
#' @return GRanges object with the set of translated readthrough regions
#' @seealso \code{\link{detect_translated_orfs}}, \code{\link{select_quantify_ORFs}}, \code{\link{annotate_ORFs}}, \code{\link{get_reathr_seq}}
#' @export

detect_readthrough<-function(results_orf,P_sites,P_sites_uniq,P_sites_uniq_mm,genome_sequence,annotation,genetic_code_table,cutoff_fr_ave=.5,uniq_signal=F){
  P_sites<-P_sites[!P_sites%over%unlist(results_orf$ORFs_genomic_position)]
  P_sites_uniq<-P_sites_uniq[!P_sites_uniq%over%unlist(results_orf$ORFs_genomic_position)]
  P_sites_uniq_mm<-P_sites_uniq_mm[!P_sites_uniq_mm%over%unlist(results_orf$ORFs_genomic_position)]
  readthroughs<-GRanges()
  if(length(P_sites)>4){
    for(i in 1:length(results_orf$ORFs_tx_position)){
      orf_tx<-results_orf$ORFs_tx_position[[i]]
      tx<-orf_tx$transcript_id
      ex_tx<-annotation$exons_txs[[tx]]
      
      #tx_gr<-transcripts_ranges[tx_ok]
      if(as.character(strand(ex_tx)[1])=="+"){
        covtx<-suppressWarnings(unlist(coverage(P_sites,weight = P_sites$score)[ex_tx]))
        if(length(P_sites_uniq)>0){
          covtx_uniq<-suppressWarnings(unlist(coverage(P_sites_uniq,weight = P_sites_uniq$score)[ex_tx]))
        }
        if(length(P_sites_uniq)==0){
          covtx_uniq<-suppressWarnings(unlist(coverage(P_sites_uniq)[ex_tx]))
        }
        if(length(P_sites_uniq_mm)>0){
          covtx_uniq_mm<-suppressWarnings(unlist(coverage(P_sites_uniq_mm,weight = P_sites_uniq_mm$score)[ex_tx]))
        }
        if(length(P_sites_uniq_mm)==0){
          covtx_uniq_mm<-suppressWarnings(unlist(coverage(P_sites_uniq_mm)[ex_tx]))
        }
        
      }
      if(as.character(strand(ex_tx)[1])=="-"){
        covtx<-suppressWarnings(unlist(unlist(RleList(lapply(coverage(P_sites,weight = P_sites$score)[ex_tx],FUN=rev)))))
        if(length(P_sites_uniq)>0){
          covtx_uniq<-suppressWarnings(unlist(unlist(RleList(lapply(coverage(P_sites_uniq,weight = P_sites_uniq$score)[ex_tx],FUN=rev)))))
        }
        if(length(P_sites_uniq)==0){
          covtx_uniq<-suppressWarnings(unlist(unlist(RleList(lapply(coverage(P_sites_uniq)[ex_tx],FUN=rev)))))
        }
        if(length(P_sites_uniq_mm)>0){
          covtx_uniq_mm<-suppressWarnings(unlist(unlist(RleList(lapply(coverage(P_sites_uniq_mm,weight = P_sites_uniq_mm$score)[ex_tx],FUN=rev)))))
        }
        if(length(P_sites_uniq_mm)==0){
          covtx_uniq_mm<-suppressWarnings(unlist(unlist(RleList(lapply(coverage(P_sites_uniq_mm)[ex_tx],FUN=rev)))))
        }
        
      }
      seq_tx<-unlist(getSeq(x=genome_sequence,ex_tx))
      
      orfs<-get_reathr_seq(tx_name = tx,orf = orf_tx,sequence = seq_tx,genetic_code = genetic_code_table)
      
      if(length(orfs)==0){next}
      vals1<-calc_orf_pval(ORFs = orfs[1],P_sites_rle = covtx,P_sites_uniq_rle = covtx_uniq,P_sites_uniq_mm_rle = covtx_uniq_mm,cutoff = cutoff_fr_ave)
      #uniq_flag
      if(!uniq_signal){
        if(is.na(vals1$pval) | vals1$pct_fr<.5 | vals1$pval>.05 ){next}
      }
      if(uniq_signal){
        if(is.na(vals1$pval_uniq) | vals1$pct_fr<.5 | vals1$pval_uniq>.05 ){next}
      }
      vals1$Protein<-AAStringSet(as.character(translate_solve(seq_tx[vals1@ranges],genetic.code = genetic_code_table)))
      vals1$ORF_orig_tr<-orf_tx$ORF_id_tr
      vals1$n_stops_readth<-1
      vals1$compatible_id<-CharacterList("")
      vals1$compatible_original_id<-CharacterList("")
      mcsva<-mcols(vals1)
      mcsva[,c("gene_id","gene_biotype","gene_name","transcript_id","transcript_biotype")]<-""
      mcols(vals1)<-mcsva
      vals1$gen_coords<-GRangesList(GRanges())
      
      if(length(orfs)>1){
        keep<-1
        while(keep>0){
          if(length(orfs)<(keep+1)){break}
          vals<-calc_orf_pval(ORFs = orfs[keep+1],P_sites_rle = covtx,P_sites_uniq_rle = covtx_uniq,P_sites_uniq_mm_rle = covtx_uniq_mm,cutoff = cutoff_fr_ave)
          keep<-keep+1
          #uniq_flag
          
          if(!uniq_signal){
            if(is.na(vals1$pval) | vals1$pct_fr<.5 | vals1$pval>.05 ){keep<-0}
            if(!is.na(vals$pval) & vals$pct_fr>.5 & vals$pval<.05 ){
              end(vals1)<-end(vals)
              vals1<-calc_orf_pval(ORFs = vals1,P_sites_rle = covtx,P_sites_uniq_rle = covtx_uniq,P_sites_uniq_mm_rle = covtx_uniq_mm,cutoff = cutoff_fr_ave)
              vals1$Protein<-AAStringSet(as.character(translate_solve(seq_tx[vals1@ranges],genetic.code = genetic_code_table)))
              vals1$ORF_orig_tr<-orf_tx$ORF_id_tr
              vals1$n_stops_readth<-keep
              vals1$compatible_id<-CharacterList("")
              vals1$compatible_original_id<-CharacterList("")
              mcsva<-mcols(vals1)
              mcsva[,c("gene_id","gene_biotype","gene_name","transcript_id","transcript_biotype")]<-""
              mcols(vals1)<-mcsva
              
            }
          }
          if(uniq_signal){
            if(is.na(vals1$pval_uniq) | vals1$pct_fr<.5 | vals1$pval_uniq>.05 ){keep<-0}
            if(!is.na(vals$pval_uniq) & vals$pct_fr>.5 & vals$pval_uniq<.05 ){
              end(vals1)<-end(vals)
              vals1<-calc_orf_pval(ORFs = vals1,P_sites_rle = covtx,P_sites_uniq_rle = covtx_uniq,P_sites_uniq_mm_rle = covtx_uniq_mm,cutoff = cutoff_fr_ave)
              vals1$Protein<-AAStringSet(as.character(translate_solve(seq_tx[vals1@ranges],genetic.code = genetic_code_table)))
              vals1$ORF_orig_tr<-orf_tx$ORF_id_tr
              vals1$n_stops_readth<-keep
              vals1$compatible_id<-CharacterList("")
              vals1$compatible_original_id<-CharacterList("")
              mcsva<-mcols(vals1)
              mcsva[,c("gene_id","gene_biotype","gene_name","transcript_id","transcript_biotype")]<-""
              mcols(vals1)<-mcsva
              
            }
          }
          
          
          
        }
      }
      #added for new R/Biocond compatibility (Dec 2021)
      if(length(readthroughs)==0){readthroughs<-vals1[-1]}
      readthroughs<-unique(suppressWarnings(c(readthroughs,vals1)))
      mcs<-mcols(readthroughs)
      mcs_orfs_orig<-mcs$ORF_orig_tr
      mcs_orfs<-mcs$ORF_id_tr
      mcs$ORF_orig_tr<-NULL
      mcs$ORF_id_tr<-NULL
      dups<-duplicated(mcs)
      readthroughs<-readthroughs[!dups]
      for(j in 1:length(readthroughs)){
        orig_orf<-readthroughs$ORF_orig_tr[j]
        mcs_orig<-mcols(results_orf$ORFs_tx_position[[which(names(results_orf$ORFs_tx_position)==orig_orf)]])[,c("gene_id","gene_biotype","gene_name","transcript_id","transcript_biotype")]
        mcols(readthroughs)[j,c("gene_id","gene_biotype","gene_name","transcript_id","transcript_biotype")]<-mcs_orig
        readthroughs$gen_coords[j]<-from_tx_togen(readthroughs[j],annotation$exons_txs[as.character(seqnames(readthroughs[j]))],annotation$introns_txs[[as.character(seqnames(readthroughs[j]))]])
        readthroughs$compatible_id[j]<-CharacterList(mcs_orfs[which(mcs$Protein==mcs[j,"Protein"])[-j]])
        readthroughs$compatible_original_id[j]<-CharacterList(mcs_orfs_orig[which(mcs$Protein==mcs[j,"Protein"])[-j]])
        
      }
      dups<-duplicated(readthroughs$Protein) & duplicated(readthroughs$P_sites_raw)
      
      if(sum(dups)>0){readthroughs<-readthroughs[!dups]}
    }
  }
  dups<-duplicated(readthroughs$Protein) & duplicated(readthroughs$P_sites_raw)
  
  if(sum(dups)>0){readthroughs<-readthroughs[!dups]}
  return(readthroughs)
}

#' Select and quantify ORF translation
#'
#' This function selects a subset of detected ORFs and quantifies their translation
#' @details ORFs are first selected using the same method as in the \code{select_txs} function, but 
#' using ORF features (ORF structures are treated as transcript structures).\cr
#' Ribo-seq coverage (reads/length) on bins and junctions (set to a length of 60) is used to derive a scaling factor (0-1) for each ORF,
#' which indicates how much of the ORF coverage can be assigned to such ORF (1 when no other ORF is present). 
#' When no unique features are present on an ORF, an adjusted scaling value is calculated subtracting coverage expected from a ORF with a unique feature. 
#' When no unique features are present on any ORF, scaling values are calculated assuming uniform coverage on each ORF.\cr
#' Scaling values are then further scaled to adjust for average coverage or for the total number of P_sites on the ORFs' exonic bins (recommended, the default).\cr
#' ORFs are then further filtered to exclude lowly translated ORFs and quantification/selection is re-iterated until no ORF is further filtered out. 
#' Percentage of total gene translation and length-adjusted quantification estimates are produced.
#' More details about the quantificatin procedure can be found in the ORFquant manuscript.\cr\cr
#' Additional columns are added to the ORFs_tx object:\cr
#' \code{P_sites}: P_sites_raw value from \code{detect_translated_orfs} (P_sites_raw_uniq with \code{uniq_signal = TRUE}) multiplied by the final ORF scaling factor.\cr
#' \code{ORF_pct_P_sites}: Percentage of gene translation output for the ORF, derived using P_sites values.\cr
#' \code{ORF_pct_P_sites_pN}: Percentage of gene translation ouptut (adjusted by length) for the ORF, derived using P_sites values.\cr
#' \code{unique_features_reads}:  initial number of reads on each unique ORF feature. \code{NA} when no unique feature is present.\cr
#' \code{adj_unique_features_reads}:  same as \code{unique_features_reads} for ORFs with unique features; for the others, number of reads on each feature shared only with ORFs that already have a scaling factor. \code{NA} when there is none.\cr
#' \code{scaling_factors}: Set of 3 scaling factors assigned to the ORF using intial unique ORF features, after adjusting for the presence of ORFs with no unique features, and final scaling factor after correcting for average Ribo-seq coverage (or total number of reads) on the ORFs.
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param results_ORFs Full list of detected ORFs, from \code{detect_translated_orfs}
#' @param P_sites GRanges object with P_sites positions
#' @param P_sites_uniq GRanges object with uniquely mapping P_sites positions
#' @param cutoff_cums cutoff to select ORFs until <x> percentage of total gene translation. Defaults to NA (not applied)
#' @param cutoff_pct minimum percentage of total gene translation for an ORF to be selected. Defaults to 2
#' @param cutoff_P_sites minimum number of P_sites assigned to the ORF to be selected. Defaults to NA (not applied)
#' @param optimiz (Beta) should numerical optimization (minimizing distance between observed coverage and expected coverage) 
#' be used to quantify ORF translation? Defaults to FALSE
#' @param scaling Additional scaling value taking into account average or total signal on the detected ORFs to adjust quantification estimates. 
#' Can be average_coverage or total_Psites. Defaults to total_Psites for consistency.
#' @param uniq_signal Use only signal from uniquely mapping reads? Defaults to \code{FALSE}.
#' @return modified \code{results_ORFs} object with the selected ORFs including quantification estimates.
#' An empty list when the cutoffs remove all ORFs.
#' @seealso \code{\link{detect_translated_orfs}}, \code{\link{select_txs}}, \code{\link{ORFquant_output}}
#' @export


select_quantify_ORFs<-function(results_ORFs,P_sites,P_sites_uniq,cutoff_cums=NA,cutoff_pct=2,cutoff_P_sites=NA,optimiz=FALSE,scaling="total_Psites",uniq_signal=F){
  
  if(!scaling%in%c("total_Psites","average_coverage")){stop(paste("scaling parameter must be either total_Psites (recommended) or average_coverage"),date())}
  select_feat <- results_ORFs[["ORFs_features"]]
  
  select_feat<-endoapply(select_feat,function(x){
    unqid<-paste(x@ranges,x$type,names(x),sep="_")
    x<-x[!duplicated(unqid)]
    x
  })
  
  select_feats<-unlist(select_feat)
  
  nmss<-c()
  for(nm in names(results_ORFs[["ORFs_features"]])){
    nmss<-c(nmss,rep(nm,length(select_feat[[nm]])))
  }
  names(select_feats)<-nmss
  select_feats_jun<-select_feats[select_feats$type=="J"]
  
  # ran_j<-IRanges()
  # allofthem<-unique(unlist(select_feats_jun$txs_orfs))
  # for(i in)
  if(length(select_feats_jun)>0){
    
    
    ran_j<-select_feats_jun@ranges
    #startend!
    df<-data.frame(stend=paste(ran_j@start,end(ran_j),sep="_"),orfs=names(ran_j),stringsAsFactors=F)
    tab_j<-as.matrix(table(df$stend,df$orfs))
    
    orfs_print<-colnames(tab_j)
    listorf_print<-list()
    for(junso in rownames(tab_j)){
      listorf_print[[junso]]<-orfs_print[which(tab_j[junso,]>0)]
    }
    
    orfs_print2<-listorf_print
    
    #names(orfs_print2)<-NULL
  }
  orfs <- results_ORFs[["ORFs_genomic_position"]]
  
  # the exonic parts of the ORFs and their tx_name, as exonicParts() gives them for a TxDb with one
  # transcript per ORF (exons in the order of orfs, on the chromosome and strand of the first ORF),
  # without building the TxDb (an SQLite database)
  ex<-unlist(orfs,use.names=FALSE)
  ex<-GRanges(as.character(unique(seqnames(orfs[[1]]))),ex@ranges,as.character(unique(strand(orfs[[1]]))))
  exbin<-disjoin(ex,with.revmap=TRUE)
  tx_name<-unique(IRanges::extractList(rep(as.character(names(orfs)),elementNROWS(orfs)),mcols(exbin)$revmap))
  mcols(exbin)<-DataFrame(tx_name=tx_name[!is.na(tx_name)])
  
  
  d<-rep(0,length(exbin))
  
  hts<-findOverlaps(exbin,P_sites,ignore.strand=F)
  hts<-cbind(queryHits(hts),P_sites[subjectHits(hts)]$score*width(P_sites[subjectHits(hts)]))
  if(length(hts)==0){
    return(GRangesList())
    
  }
  hts<-aggregate(x = hts[,2],list(hts[,1]),FUN=sum)
  for(i in 1:dim(hts)[1]){
    d[hts[i,1]]<-hts[i,2]
  }
  
  d2<-rep(0,length(exbin))
  
  hts<-findOverlaps(exbin,P_sites_uniq,ignore.strand=F)
  hts<-cbind(queryHits(hts),P_sites_uniq[subjectHits(hts)]$score*width(P_sites_uniq[subjectHits(hts)]))
  if(length(hts)>0){
    hts<-aggregate(x = hts[,2],list(hts[,1]),FUN=sum)
    for(i in 1:dim(hts)[1]){
      d2[hts[i,1]]<-hts[i,2]
    }
  }
  
  
  # columns are set on mcols() and put back once: a GRanges $<- also runs updateObject() (~20 ms)
  cols<-mcols(exbin)
  cols$reads<-d
  cols$unique_reads<-d2
  mcols(exbin)<-cols
  
  exbin<-GRanges(exbin)
  
  cols<-mcols(exbin)
  cols$gene_id<-NULL
  cols$exonic_part<-NULL
  cols$X<-NULL
  mcols(exbin)<-cols
  mcols(exbin)<-exbin[,c("reads","unique_reads","tx_name")]
  gene_feat<-exbin
  mcols(gene_feat)$type<-rep("E",length(exbin))
  
  if(length(select_feats_jun)>0){
    junc<-select_feats_jun
    names(junc)<-NULL
    reads_j<-junc$reads
    reads_j_uniq<-junc$unique_reads
    mcols(junc)<-NULL
    cols<-mcols(junc)
    cols$X<-junc
    cols$reads<-reads_j
    cols$unique_reads<-reads_j_uniq
    mcols(junc)<-cols
    
    orfs_jj<-orfs_print2[match(df$stend,names(orfs_print2))]
    #junc$tx_name<-as(orfs_jj,"CharacterList")
    mcols(junc)<-junc[,c("reads","unique_reads")]
    tx_ex<-exbin$tx_name
    mcols(exbin)<-exbin[,c("reads","unique_reads")]
    
    gene_feat<-c(exbin,junc)
    cols<-mcols(gene_feat)
    if(is(orfs_jj,"CompressedList")){
      cols$tx_name<-CharacterList(c(tx_ex,orfs_jj))
    }
    
    if(!is(orfs_jj,"CompressedList")){
      if(is.list(orfs_jj)){
        cols$tx_name<-CharacterList(c(tx_ex,CharacterList(orfs_jj)))
      }
      if(!is.list(orfs_jj)){
        cols$tx_name<-CharacterList(c(tx_ex,CharacterList(as.list(orfs_jj))))
      }
    }
    
    cols$type<-c(rep("E",length(exbin)),rep("J",length(junc)))
    mcols(gene_feat)<-cols
  }
  
  #first round
  gene_feat_cp<-gene_feat
  #change!
  gene_feat<-gene_feat_cp[!duplicated(paste(GRanges(gene_feat_cp),gene_feat_cp$type,sep="_"))]
  txs_gene<-unique(unlist(gene_feat$tx_name))
  txs_sofar<-txs_gene
  
  #uniq_flag
  
  if(uniq_signal){
    d<-gene_feat$unique_reads
  }
  if(!uniq_signal){
    d<-gene_feat$reads
  }
  
  a<-gene_feat$tx_name
  
  mat<-matrix(data=0,nrow=length(d),ncol=length(txs_sofar))
  colnames(mat)<-txs_sofar
  for(i in 1:length(txs_sofar)){
    mat[,i]<-sapply(a,function(x){sum(x==txs_sofar[i])})
  }
  
  #TAKE AWAY NESTED TXS
  nest<-c()
  ident<-c()
  for(i in 1:dim(mat)[2]){
    yes<-which(mat[,i]==1)
    nesti<-c()
    for(j in (1:dim(mat)[2])[-i]){
      yesj<-which(mat[,j]==1)
      if(identical(yes,yesj)){ident<-c(ident,paste(sort(colnames(mat)[c(i,j)]),collapse=";"))}
      #added if length> otherwise txs with same structure are both deleted
      nesti<-c(nesti,sum(yesj%in%yes)==length(yes) & length(mat[,j])>length(yes))
    }
    
    if(sum(nesti)>0){nest<-c(nest,colnames(mat)[i])}
    
  }
  if(length(ident)>0){nest<-nest[!nest%in%unique(sapply(strsplit(ident,";"),"[[",1))]}
  txs_sofar<-txs_sofar[!txs_sofar%in%nest]
  change<-1
  
  while(change>0){
    
    mat<-matrix(data=0,nrow=length(d),ncol=length(txs_sofar))
    colnames(mat)<-txs_sofar
    for(i in 1:length(txs_sofar)){
      mat[,i]<-sapply(a,function(x){sum(x==txs_sofar[i])})
    }
    
    d_count<-paste(1:length(d),d,sep="_")
    good<-d_count[which(d>0)]
    bad<-d_count[which(d==0)]
    
    txs_good<-c()
    expl_good<-c()
    for(i in 1:dim(mat)[2]){
      expl_good_old<-expl_good
      #if new good feature
      tx<-d_count[which(mat[,i]>0)]
      tx_good<-tx[which(tx%in%good)]
      tx_bad<-tx[which(tx%in%bad)]
      
      if(length(tx_good)==0){next}
      
      if(sum(!tx_good%in%expl_good_old)>0){
        
        tx_torem<-c()
        tx_toscreen<-which(colnames(mat)%in%txs_good)
        if(length(tx_toscreen)>0){
          for(j in tx_toscreen){
            tx_contr<-d_count[which(mat[,j]>0)]
            tx_contr_good<-tx_contr[which(tx_contr%in%good)]
            tx_contr_bad<-tx_contr[which(tx_contr%in%bad)]
            if(sum(tx_contr_good%in%tx_good)==length(tx_contr_good)){
              tx_torem<-c(tx_torem,colnames(mat)[j])
              
            }
            
            if(length(tx_torem)>0){
              txs_good<-txs_good[!txs_good%in%tx_torem]
              
            }
          }
          
          
          
        }
        txs_good<-unique(c(txs_good,colnames(mat)[i]))
        expl_good<-unique(c(expl_good,tx_good))
      }
      #if same good feature, but fewer bad INTERNAL features than others.
      if(sum(!tx_good%in%expl_good_old)==0){
        tx_torem<-c()
        tx_toscreen<-which(colnames(mat)%in%txs_good)
        #here a counter when at least one good feature more than competing
        moref<-c()
        for(j in tx_toscreen){
          tx_contr<-d_count[which(mat[,j]>0)]
          tx_contr_good<-tx_contr[which(tx_contr%in%good)]
          tx_contr_bad<-tx_contr[which(tx_contr%in%bad)]
          moref<-c(moref,sum(!tx_good%in%tx_contr_good)>0)
          #if same good features in competing tx, check bad internal ones
          if(sum(tx_contr_good%in%tx_good)==length(tx_contr_good)){
            
            if(length(tx_good)>length(tx_contr_good)){
              tx_torem<-c(tx_torem,colnames(mat)[j])
            }
            #here
            if(length(tx_good)==length(tx_contr_good)){      
              fi<-which(tx==tx_good[1])
              la<-which(tx==tx_good[length(tx_good)])
              int_tx<-tx[fi:la]
              int_tx_bad<-tx_bad[tx_bad%in%int_tx]
              
              contr_fi<-which(tx_contr==tx_contr_good[1])
              contr_la<-which(tx_contr==tx_contr_good[length(tx_contr_good)])
              contr_int_tx<-tx_contr[contr_fi:contr_la]
              contr_int_tx_bad<-tx_contr_bad[tx_contr_bad%in%contr_int_tx]
              #SAME INTERNAL, TAKE
              if(length(int_tx_bad)<=length(contr_int_tx_bad)){
                
                txs_good<-unique(c(txs_good,colnames(mat)[i]))
                expl_good<-unique(c(expl_good,tx_good))
                #LESS INTERNAL, TAKE AND REMOVE OTHER
                if(length(int_tx_bad)<length(contr_int_tx_bad)){
                  
                  tx_torem<-c(tx_torem,colnames(mat)[j])}
              }
              
              
            }
          }
          
          
        }
        if(length(tx_torem)>0){
          txs_good<-txs_good[!txs_good%in%tx_torem]
          txs_good<-unique(c(txs_good,colnames(mat)[i]))
          expl_good<-unique(c(expl_good,tx_good))
        }
        
        if(length(tx_torem)==0 & sum(moref)==length(tx_toscreen)){
          txs_good<-unique(c(txs_good,colnames(mat)[i]))
          expl_good<-unique(c(expl_good,tx_good))
        }
        
        
      }
      
      
    }
    txs_sofar<-txs_good
    change<-abs(length(txs_good)-length(txs_sofar))
  }
  
  cols<-mcols(gene_feat)
  cols$ORF_id_tr<-cols$tx_name
  cols$tx_name<-NULL
  a<-cols$ORF_id_tr
  a<-lapply(a,function(x){x[x%in%txs_good]})
  cols$ORF_id_tr_selected<-CharacterList(a)
  
  check<-sapply(cols$ORF_id_tr,FUN = length)
  use<-rep("shared",length(check))
  use[check==1]<-"unique"
  use[check==0]<-"absent"
  use[check>1]<-"shared"
  cols$use_ORF<-use
  
  check<-sapply(cols$ORF_id_tr_selected,FUN = length)
  use<-rep("shared",length(check))
  use[check==1]<-"unique"
  use[check==0]<-"absent"
  use[check>1]<-"shared"
  cols$use_ORF_selected<-use
  mcols(gene_feat)<-cols
  
  sel_feats<-list()
  for(i in txs_good){
    featexs<-gene_feat[gene_feat$type=="E"]
    featjuns<-gene_feat[gene_feat$type=="J"]
    
    
    ok<-featexs[featexs%over%orfs[[i]]]
    if(length(featjuns)>0){
      ok<-sort(c(ok,featjuns[which(featjuns%in%gaps(orfs[[i]]))]))
    }
    a<-sapply(ok$ORF_id_tr,function(x){length(x[x%in%i])})
    ok<-ok[a>0]
    sel_feats[[i]]<-unique(ok)
  }
  selected_ORFs<-lapply(results_ORFs,FUN=function(x){x[names(x)%in%txs_good]})
  selected_ORFs$selected_ORFs_features<-sel_feats
  
  #quantif 
  
  
  results_ORFs<-selected_ORFs
  feats<-results_ORFs$selected_ORFs_features
  orfs_tx<-results_ORFs$ORFs_tx_position
  orfs_tx<-lapply(orfs_tx,function(x){
    cols<-mcols(x)
    cols[,c("P_sites","ORF_pct_P_sites","ORF_pct_P_sites_pN")]<-NA
    cols[,"unique_features_reads"]<-NumericList("")
    cols[,"adj_unique_features_reads"]<-NumericList("")
    cols[,"scaling_factors"]<-NumericList("")
    mcols(x)<-cols
    x
  })
  
  #first round of unq
  
  orf_del<-c("")
  counter<-1
  cols<-mcols(gene_feat)
  cols$ORF_id_tr_selected_quant<- cols$ORF_id_tr_selected
  cols$use_ORF_selected_quant<-cols$use_ORF_selected
  mcols(gene_feat)<-cols
  
  while(length(orf_del)>0){
    
    feats<-feats[!names(feats)%in%orf_del]
    #the cutoffs removed all ORFs: the region has no ORFs, as when detect_translated_orfs finds none
    if(length(feats)==0){return(list())}
    nms<-sapply(feats,length)
    nms<-rep(names(nms),nms)
    feats<-unlist(GRangesList(unlist(feats)))
    names(feats)<-nms
    
    cols_f<-mcols(feats)
    cols_f$ORF_id_tr_selected<-CharacterList(lapply(cols_f$ORF_id_tr_selected,function(x){
      unique(x[!x%in%orf_del])
    }))
    
    cols<-mcols(gene_feat)
    cols$ORF_id_tr_selected_quant<-CharacterList(lapply(cols$ORF_id_tr_selected_quant,function(x){
      unique(x[!x%in%orf_del])
    }))
    lens<-sapply(cols_f$ORF_id_tr_selected,length)
    lens2<-sapply(cols$ORF_id_tr_selected_quant,length)
    cols$use_ORF_selected_quant<-"shared"
    cols$use_ORF_selected_quant[lens2==1]<-"unique"
    cols$use_ORF_selected_quant[lens2==0]<-"absent"
    mcols(gene_feat)<-cols
    
    cols_f$use_ORF_selected<-"shared"
    cols_f$use_ORF_selected[lens==1]<-"unique"
    mcols(feats)<-cols_f
    feats<-split(feats,names(feats))
    
    orfs_tx<-orfs_tx[!names(orfs_tx)%in%orf_del]
    
    
    unqs<-rep(NA,length(feats))
    names(unqs)<-names(feats)
    for(i in names(feats)){
      feat<-feats[[i]]
      orf_tx<-orfs_tx[[i]]
      
      #uniq_flag
      if(!uniq_signal){
        riz<-feat$reads
      }
      
      if(uniq_signal){
        riz<-feat$unique_reads
      }
      
      
      cov_feat<-riz/width(feat)
      js<-feat$type=="J"
      if(sum(js)>0){
        cov_feat[js]<-riz[js]/60
      }
      unq<-feat$use_ORF_selected=="unique"
      if(sum(unq)>0){
        unq_rat<-mean(cov_feat[unq])/mean(cov_feat)
        if(unq_rat>1){unq_rat<-1}
        cols<-mcols(orfs_tx[[i]])
        cols$unique_features_reads<-NumericList(riz[unq])
        cols$adj_unique_features_reads<-NumericList(riz[unq])
        mcols(orfs_tx[[i]])<-cols
        
        
      }
      if(sum(unq)==0){
        mcols(orfs_tx[[i]])$unique_features_reads<-NumericList(NA)
        unq_rat<-NA
      }
      unqs[i]<-unq_rat
    }
    unqs_adj<-unqs
    unqnas<-unqs_adj[is.na(unqs_adj)]
    unqs_okk<-unqs_adj[!is.na(unqs_adj)]
    #NA_adjustment
    
    #if all NAs
    if(length(unqs_okk)==0){
      unqs_na<-names(unqs_adj[is.na(unqs_adj)])
      for(i in unqs_na){
        feat<-feats[[i]]
        
        riz<-feat$reads
        
        if(uniq_signal){
          riz<-feat$unique_reads
        }
        
        cov_feat<-riz/width(feat)
        
        js<-feat$type=="J"
        if(sum(js)>0){
          cov_feat[js]<-riz[js]/60
        }
        cov_adj<-cov_feat
        for(j in 1:length(feat)){
          fea<-feat[j]
          txs_fea<-unlist(fea$ORF_id_tr_selected)
          nass<-txs_fea%in%unqs_na
          txs_fea<-txs_fea[!nass]
          adj<-cov_feat[j]-(cov_feat[j]*sum(unqs_adj[txs_fea]))
          if(sum(nass)>1){
            adj<-(cov_feat[j]-(cov_feat[j]*sum(unqs_adj[txs_fea])))/sum(nass)
          }
          if(length(txs_fea)==0){
            adj<-cov_feat[j]/(sum(nass))
          }
          if(adj<0){adj=0}
          cov_adj[j]<-adj
        }
        unq<-mean(cov_adj)/mean(cov_feat)
        if(unq>1){unq<-1}
        unqs_adj[i]<-unq
      }
    }
    unqnas<-unqs_adj[is.na(unqs_adj)]
    unqs_okk<-unqs_adj[!is.na(unqs_adj)]
    
    #if all 0
    
    if(sum(unqs_okk==0)==length(unqs_adj)){
      unqs_zero<-names(unqs_okk)
      for(i in unqs_zero){
        feat<-feats[[i]]
        
        riz<-feat$reads
        
        if(uniq_signal){
          riz<-feat$unique_reads
        }
        
        cov_feat<-riz/width(feat)
        
        js<-feat$type=="J"
        if(sum(js)>0){
          cov_feat[js]<-riz[js]/60
        }
        cov_adj<-cov_feat
        for(j in 1:length(feat)){
          fea<-feat[j]
          txs_fea<-unlist(fea$ORF_id_tr_selected)
          nass<-txs_fea%in%unqs_zero
          txs_fea<-txs_fea[!nass]
          adj<-cov_feat[j]-(cov_feat[j]*sum(unqs_okk[txs_fea]))
          if(sum(nass)>1){
            adj<-(cov_feat[j]-(cov_feat[j]*sum(unqs_okk[txs_fea])))/sum(nass)
          }
          if(length(txs_fea)==0){
            adj<-cov_feat[j]/(sum(nass))
          }
          if(adj<0){adj=0}
          cov_adj[j]<-adj
        }
        unq<-mean(cov_adj)/mean(cov_feat)
        if(unq>1){unq<-1}
        unqs_okk[i]<-unq
        unqs_adj[i]<-unq
        
      }
    }
    unqnas<-unqs_adj[is.na(unqs_adj)]
    unqs_okk<-unqs_adj[!is.na(unqs_adj)]
    
    if(length(unqnas)>0){
      nonas<--1
      
      while(nonas<0){
        
        prev_nas<-length(unqnas)
        if(length(unqs_okk)>0){
          for(i in names(unqnas)){
            unqs_okk_noi<-unqs_okk[names(unqs_okk)!=i]
            feat<-feats[[i]]
            orf_tx<-orfs_tx[[i]]
            
            riz<-feat$reads
            
            if(uniq_signal){
              riz<-feat$unique_reads
            }
            
            cov_feat<-riz/width(feat)
            
            js<-feat$type=="J"
            if(sum(js)>0){
              cov_feat[js]<-riz[js]/60
            }
            adj_use<-feat$use_ORF_selected
            adj_cov_feat<-cov_feat
            for(j in 1:length(feat)){
              txs_fea<-feat$ORF_id_tr_selected[[j]]
              unqs_okk_noi_feat<-unqs_okk_noi[names(unqs_okk_noi)%in%txs_fea]
              txs_fea_adj<-txs_fea[!txs_fea%in%names(unqs_okk_noi)]
              adj<-cov_feat[j]
              usef<-"shared"
              if(length(txs_fea_adj)==1){usef<-"unique"}
              
              if(length(unqs_okk_noi_feat[!is.na(unqs_okk_noi_feat)])>0){
                adj<-adj-(adj*sum(unqs_okk_noi_feat[!is.na(unqs_okk_noi_feat)]))     
                if(adj<0){adj=0}
              }
              adj_use[j]<-usef
              adj_cov_feat[j]<-adj
            }
            
            unq<-adj_use=="unique"
            
            if(sum(unq)>0){
              unq_rat<-mean(adj_cov_feat[adj_use=="unique"])/mean(adj_cov_feat)
              if(mean(adj_cov_feat)==0){unq_rat<-0}
              if(unq_rat>1){unq_rat<-1}
              mcols(orfs_tx[[i]])$adj_unique_features_reads<-NumericList(riz[unq])
              
            }
            if(sum(unq)==0){
              mcols(orfs_tx[[i]])$adj_unique_features_reads<-NumericList(NA)
              unq_rat<-NA
            }
            unqs_adj[i]<-unq_rat
            
            
          }
          unqnas<-unqs_adj[is.na(unqs_adj)]
          unqs_okk<-unqs_adj[!is.na(unqs_adj)]
          nonas<-length(unqnas)-prev_nas
          
        }
        
        
      }
      #if still NAs
      if(sum(is.na(unqs_adj))>0){
        unqs_na<-names(unqs_adj[is.na(unqs_adj)])
        for(i in unqs_na){
          feat<-feats[[i]]
          
          riz<-feat$reads
          
          if(uniq_signal){
            riz<-feat$unique_reads
          }
          
          cov_feat<-riz/width(feat)
          
          
          js<-feat$type=="J"
          if(sum(js)>0){
            cov_feat[js]<-riz[js]/60
          }
          cov_adj<-cov_feat
          for(j in 1:length(feat)){
            fea<-feat[j]
            txs_fea<-unlist(fea$ORF_id_tr_selected)
            nass<-txs_fea%in%unqs_na
            txs_fea<-txs_fea[!nass]
            adj<-cov_feat[j]-(cov_feat[j]*sum(unqs_adj[txs_fea]))
            if(sum(nass)>1){
              adj<-(cov_feat[j]-(cov_feat[j]*sum(unqs_adj[txs_fea])))/sum(nass)
            }
            if(length(txs_fea)==0){
              adj<-cov_feat[j]/(sum(nass))
            }
            if(adj<0){adj=0}
            cov_adj[j]<-adj
          }
          unq<-mean(cov_adj)/mean(cov_feat)
          if(unq>1){unq<-1}
          unqs_adj[i]<-unq
        }
        
      }
      
      
    }
    
    
    if(sum(unqs_adj>0)==0){
      nmmm<-names(unqs_adj)
      unqs_adj<-rep(1/length(nmmm),length(nmmm))
      names(unqs_adj)<-nmmm
    }
    
    orftxs<-unlist(results_ORFs[["ORFs_tx_position"]])
    unqs_optim<-unqs_adj
    
    if(!is(feats,"GRangesList")){feats<-GRangesList(feats)}
    featsall<-sort(unlist(feats))
    featsall<-featsall[!duplicated(paste(GRanges(featsall),featsall$type,sep="_"))]
    
    #optimization option
    
    if(optimiz==TRUE){
      unqs_optim_noopt<-unqs_optim
      #uniq_flag
      
      if(!uniq_signal){
        cov_ps_unqs<-orftxs$P_sites_raw/width(orftxs)
      }
      
      if(uniq_signal){
        cov_ps_unqs<-orftxs$P_sites_raw_uniq/width(orftxs)
      }
      
      
      names(cov_ps_unqs)<-orftxs$ORF_id_tr
      cov_ps_unqs<-cov_ps_unqs[names(unqs_adj)]
      #uniq_flag
      
      if(!uniq_signal){
        featsall$coverage<-featsall$reads/width(featsall)
        jxs<-featsall$type=="J"
        featsall$coverage[jxs]<-featsall$reads[jxs]/60
      }
      
      if(uniq_signal){
        featsall$coverage<-featsall$unique_reads/width(featsall)
        jxs<-featsall$type=="J"
        featsall$coverage[jxs]<-featsall$unique_reads[jxs]/60
      }
      
      vals<-cov_ps_unqs*unqs_optim
      expect<-unlist(lapply(featsall$ORF_id_tr_selected,function(x){sum(vals[x])}))
      maxdist<-sum(abs(expect-featsall$coverage))
      
      calc_dist_exp<-function(unqs_optim,maxval=maxdist){
        vals<-cov_ps_unqs*unqs_optim
        #sizzs<-width(featsall)
        #sizzs[featsall$type=="J"]<-60
        
        expect<-unlist(lapply(featsall$ORF_id_tr_selected,function(x){sum(vals[x])}))
        if(sum(expect)>0){
          expect<-expect/sum(expect)
        }
        trucov<-featsall$coverage
        if(sum(expect)>0){
          trucov<-trucov/sum(trucov)
        }
        #or do something different, like minimize % error per each feature after adding pseudocount
        if(sum(trucov>0 & expect==0)>0){return(maxval+1)}
        
        #sum(abs(expect-featsall$coverage)*log(sizzs+1))
        #-cor.test(trucov,expect,method = "p")$estimate
        mean(abs(expect-trucov))
      }
      #inspired by Alpine: https://github.com/mikelove/alpine/blob/master/R/estimate_abundance.R
      optimm<-optim(par = unqs_optim,fn = calc_dist_exp,lower = rep(0,length(unqs_optim)),upper = rep(1,length(unqs_optim)),method = "L-BFGS-B")
      unqs_optim<-optimm$par
      if(sum(unqs_optim>0)==0){unqs_optim<-unqs_optim_noopt}
    }
    
    if(scaling=="average_coverage"){
      
      #scale to adjust coverage?
      
      #uniq_flag
      
      if(!uniq_signal){
        cov_ps_unqs<-orftxs$P_sites_raw/width(orftxs)
      }
      
      if(uniq_signal){
        cov_ps_unqs<-orftxs$P_sites_raw_uniq/width(orftxs)
      }
      
      names(cov_ps_unqs)<-orftxs$ORF_id_tr
      cov_ps_unqs<-cov_ps_unqs[names(unqs_optim)]
      
      #uniq_flag
      if(uniq_signal){
        covo<-featsall$unique_reads/width(featsall)
      }
      
      if(!uniq_signal){
        covo<-featsall$reads/width(featsall)
      }
      featsall$coverage<-covo
      jxs<-featsall$type=="J"
      featsall$coverage[jxs]<-covo[jxs]/60
      
      cov_ps_unqs_adj<-cov_ps_unqs*unqs_optim
      cov_unpcts<-NumericList(lapply(featsall$ORF_id_tr_selected,function(x){cov_ps_unqs_adj[x]}))
      cov_unpcts<-sum(cov_unpcts)
      
      scale_f<-sum(featsall$coverage)/sum(cov_unpcts)
      if(sum(scale_f>1)>0){scale_f<-scale_f/max(scale_f)}
      unqs_optim<-unqs_optim*scale_f
    }
    
    if(scaling=="total_Psites"){
      
      unqs_tott<-unqs_optim
      orfgrps<-list()
      for(grpo in 1:length(featsall)){
        fezzo=featsall$ORF_id_tr_selected[[grpo]]
        presalready<-fezzo%in%unlist(orfgrps)
        if(sum(presalready)==0){orfgrps<-append(orfgrps,fezzo);next()}
        if(sum(presalready)>0){
          presalready<-which(orfgrps%in%fezzo)
          for(pres in presalready){
            orfgrps[[pres]]<-sort(unique(c(orfgrps[[pres]],fezzo)))
          }
        }
        
      }
      
      orfgrps<-unique(orfgrps)
      if(length(orfgrps)>1){
        orfgrpsdello<-rep(FALSE,length(orfgrps))
        for(grpo in 1:length(orfgrps)){
          orfes<-orfgrps[[grpo]]
          otherse<-orfgrps[-grpo]
          remo<-FALSE
          for(garpo in 1:length(otherse)){
            contest<-otherse[[garpo]]
            if(sum(orfes%in%contest)==length(orfes)){remo=T;break}
          }
          orfgrpsdello[grpo]=remo
        }
        
        orfgrps<-orfgrps[!orfgrpsdello]
      }
      
      ps_feats_tott<-c()
      fetse<-featsall[featsall$type=="E"]
      rizz<-fetse$reads
      ps_tott<-orftxs$P_sites_raw
      names(ps_tott)<-orftxs$ORF_id_tr
      ps_tott<-ps_tott[names(unqs_tott)]
      
      if(uniq_signal){
        ps_tott<-orftxs$P_sites_raw_uniq
        names(ps_tott)<-orftxs$ORF_id_tr
        ps_tott<-ps_tott[names(unqs_tott)]
        rizz<-fetse$unique_reads
      }
      
      for(grpo in 1:length(orfgrps)){
        ps_feats_tott<-c(ps_feats_tott,sum(rizz[sum(fetse$ORF_id_tr_selected%in%orfgrps[[grpo]])>0]))
      }
      
      adj_psts<-unqs_tott*ps_tott
      final_scls<-unqs_tott
      
      for(grpo in 1:length(orfgrps)){
        adj_grp<-adj_psts[orfgrps[[grpo]]]
        final_scls[orfgrps[[grpo]]]<-ps_feats_tott[grpo]/sum(adj_grp)
      }
      final_scls[is.infinite(final_scls)]<-0
      
      # to be added?
      # final_scls[is.infinite(final_scls)]<-0
      # final_scls[is.na(final_scls)]<-0
      # final_scls[is.nan(final_scls)]<-0
      # unqs_optim<-final_scls*unqs_tott
      # 
      # unqs_optim[is.infinite(unqs_optim)]<-0
      # unqs_optim[is.na(unqs_optim)]<-0
      # unqs_optim[is.nan(unqs_optim)]<-0
      
      
      unqs_optim<-final_scls*unqs_tott
      
      unqs_optim[unqs_optim>1]<-1
      
      
    }
    #apply quantification factor
    
    for(i in names(feats)){
      feat<-feats[[i]]
      orf_tx<-orfs_tx[[i]]
      
      #uniq_flag
      if(!uniq_signal){
        ps<-orf_tx$P_sites_raw
      }
      
      if(uniq_signal){
        ps<-orf_tx$P_sites_raw_uniq
      }
      unq_rat<-unqs_optim[i]
      
      ps_norm<-ps*unq_rat
      
      cols<-mcols(orfs_tx[[i]])
      cols$P_sites<-ps_norm
      scalss<-c(unqs[i],unqs_adj[i],unqs_optim[i])
      names(scalss)<-c("unq_feats","adj_feats","optim_feats")
      cols$scaling_factors<-NumericList(round(scalss,digits = 4))
      mcols(orfs_tx[[i]])<-cols
      
    }
    
    orfs_genes<-split(unlist(GRangesList(orfs_tx)),f=unlist(GRangesList(orfs_tx))$gene_id)
    orfs_genes<-GRangesList(lapply(orfs_genes,FUN=function(x){
      xnot<-x[is.na(x$P_sites)]
      x<-x[!is.na(x$P_sites)]
      cols<-mcols(x)
      cols$ORF_pct_P_sites<-cols$P_sites*100/sum(cols$P_sites)
      cols$ORF_pct_P_sites_pN<-(cols$P_sites/width(x))*100/sum(cols$P_sites/width(x))
      #added this when genes, mostly overlapping ones, get no reads
      
      if(sum(cols$P_sites)==0){cols$ORF_pct_P_sites<-0;cols$ORF_pct_P_sites_pN<-0}
      mcols(x)<-cols
      
      c(x[order(x$ORF_pct_P_sites,decreasing=T)],xnot)
      
      
    }))
    
    orf_del_cums<-c()
    orf_del_iso<-c()
    orf_del_ps<-c()
    
    cums_list<-list()
    for(h in names(orfs_genes)){
      x<-orfs_genes[[h]]
      cums<-cumsum(x$ORF_pct_P_sites)
      names(cums)<-x$ORF_id_tr
      cums_list[[h]]<-cums
    }
    if(is.numeric(cutoff_cums)){
      orf_del_cums<-lapply(cums_list,function(x){
        del<-which(x>cutoff_cums)
        delna<-which(is.na(x))
        if(length(del)==1){del<-c()}
        if(length(del)>1){del<-del[-1]}
        if(length(delna)==1){del<-c(del,delna)}
        del
        
      })
      names(orf_del_cums)<-NULL
      orf_del_cums<-names(unlist(orf_del_cums))
    }
    
    if(is.numeric(cutoff_pct)){
      orfgrlun<-unlist(orfs_genes)
      orf_del_iso<-unique(c(orf_del_iso,orfgrlun$ORF_id_tr[which(orfgrlun$ORF_pct_P_sites<cutoff_pct)]))
    }
    
    if(is.numeric(cutoff_P_sites)){
      orfgrlun<-unlist(orfs_genes)
      orf_del_ps<-orfgrlun$ORF_id_tr[which(orfgrlun$P_sites<cutoff_P_sites)]
      
    }
    
    orf_del<-unique(c(orf_del_cums,orf_del_iso,orf_del_ps))
    
    cols<-mcols(gene_feat)
    cols$ORF_id_tr_selected_quant<-CharacterList(lapply(cols$ORF_id_tr_selected_quant,function(x){
      unique(x[!x%in%orf_del])
    }))
    lens2<-sapply(cols$ORF_id_tr_selected_quant,length)
    cols$use_ORF_selected_quant<-"shared"
    cols$use_ORF_selected_quant[lens2==1]<-"unique"
    cols$use_ORF_selected_quant[lens2==0]<-"absent"
    mcols(gene_feat)<-cols
    
    fs<-unlist(orfs_genes)
    #put isovalues
    for(g in names(orfs_tx)){
      cols<-mcols(orfs_tx[[g]])
      cols$ORF_pct_P_sites<-round(fs$ORF_pct_P_sites[fs$ORF_id_tr==g],digits = 4)
      cols$ORF_pct_P_sites_pN<-round(fs$ORF_pct_P_sites_pN[fs$ORF_id_tr==g],digits = 4)
      mcols(orfs_tx[[g]])<-cols
    }
    counter<-counter+1
  }
  
  results_ORFs$ORFs_tx_position<-orfs_tx
  results_ORFs$ORFs_genomic_position<-results_ORFs$ORFs_genomic_position[names(results_ORFs$ORFs_tx_position)]
  results_ORFs$ORFs_features<-results_ORFs$ORFs_features[names(results_ORFs$ORFs_tx_position)]
  results_ORFs$tx_annotated_ORFs<-results_ORFs$tx_annotated_ORFs[names(results_ORFs$ORFs_tx_position)]
  #adjust feats!
  
  results_ORFs$selected_ORFs_features<-gene_feat
  return(results_ORFs)
  
}
