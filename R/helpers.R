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

# Helpers used by several steps: faster versions of Bioconductor functions, the parallel
# helpers (mc_*), from_tx_togen() and the circular chromosomes.


setMethods(f = "unlist","GRanges",function(x){return(x)})


DEFAULT_CIRC_SEQS <- unique(c("chrM","MT","MtDNA","mit","Mito","mitochondrion",
                              "dmel_mitochondrion_genome","Pltd","ChrC","Pt","chloroplast",
                              "Chloro","2micron","2-micron","2uM",
                              "Mt", "NC_001879.2", "NC_006581.1","ChrM"))

# translate(x, ..., if.fuzzy.codon = "solve"), faster. "solve" costs a fixed
# ~80 ms per call, to build the lookup of all fuzzy codons; a sequence with
# only A/C/G/T has none and gives the same protein with the default ("error").
translate_solve<-function(x,...){
  translate(x,...,if.fuzzy.codon=if(hasOnlyBaseLetters(x)) "error" else "solve")
}

# GRangesList(x) for a list x of GRanges, faster. GRangesList() binds the
# elements with c(), which makes a few S4 calls per element and mcols column,
# and merges their Seqinfo one element at a time: minutes for the thousands of
# regions or ORFs that run_ORFquant exports. bind_GRanges does the same bind on
# the slots; where it can't (see there), GRangesList() is called.
GRangesList_fast<-function(x){
  #GRangesList() suppresses the warnings of its bind too
  unl<-suppressWarnings(bind_GRanges(unname(x)))
  if(is.null(unl)){
    return(GRangesList(x))
  }
  relist(unl,IRanges::PartitioningByEnd(x))
}

# do.call(c,x) for a list x of GRanges, or NULL where this can't give the same
# result: x has other classes or only empty elements; an element has only some
# of the mcols columns (empty elements may have none); a column is Rle, factor
# or list, or its class differs between elements (unless all are plain
# vectors); or the Seqinfo don't merge. As in c(), each column is bound with
# bindROWS() after dropping empty elements; GRanges and GRangesList columns
# with bind_GRanges where it can. Seqinfo are merged pairwise in a tree, which
# gives the same seqlevels in the same order as merging them one by one.
bind_GRanges<-function(x){
  if(length(x)==0 || !all(vapply(x,function(g){class(g)[1]},"")=="GRanges")){
    return(NULL)
  }
  bind_rows<-function(v){
    v<-v[!vapply(v,is.null,TRUE)]
    if((length(unique(lapply(v,class)))>1 && !all(vapply(v,function(i){is.atomic(i) && !is.factor(i)},TRUE))) ||
       is.list(v[[1]]) || is(v[[1]],"Rle") || is.factor(v[[1]])){
      return(NULL)
    }
    v_ne<-v[vapply(v,NROW,1L)>0]
    if(length(v_ne)==0){
      return(v[[1]])
    }
    if(length(v_ne)==1){
      return(v_ne[[1]])
    }
    if(is(v_ne[[1]],"GRanges")){
      unl<-bind_GRanges(v_ne)
      if(!is.null(unl)){
        return(unl)
      }
    }
    #a CompressedGRangesList binds its unlistData with bindROWS() (empty elements kept) and has no other
    #parallel slots than mcols, here without columns
    if(is(v_ne[[1]],"CompressedGRangesList") && all(vapply(v_ne,function(i){class(i@elementMetadata)[1]=="DFrame" && length(i@elementMetadata)==0 && is.null(i@elementMetadata@rownames)},TRUE))){
      unl<-bind_GRanges(lapply(v_ne,function(i){i@unlistData}))
      if(!is.null(unl)){
        ans<-v_ne[[1]]
        ans@unlistData<-unl
        ans@partitioning<-IRanges::PartitioningByEnd(cumsum(unlist(lapply(v_ne,S4Vectors::elementNROWS))))
        ans@elementMetadata@nrows<-length(ans@partitioning)
        return(ans)
      }
    }
    S4Vectors::bindROWS(v_ne[[1]],v_ne[-1])
  }
  len<-vapply(x,length,1L)
  mc<-lapply(x,function(g){g@elementMetadata})
  nc<-vapply(mc,function(m){length(m@listData)},1L)
  cn<-names(mc[[which.max(nc>0)]]@listData)
  ok<-vapply(seq_along(mc),function(i){
    class(mc[[i]])[1]=="DFrame" && is.null(mc[[i]]@rownames) &&
      ((nc[i]==length(cn) && setequal(names(mc[[i]]@listData),cn)) || (nc[i]==0 && len[i]==0))
  },TRUE)
  if(sum(len)==0 || anyDuplicated(cn) || !all(ok)){
    return(NULL)
  }
  cols<-lapply(setNames(cn,cn),function(j){bind_rows(lapply(mc,function(m){m@listData[[j]]}))})
  rngs<-bind_rows(lapply(x,function(g){g@ranges}))
  if(any(vapply(cols,is.null,TRUE)) || is.null(rngs)){
    return(NULL)
  }
  #merging is needed (it rebuilds the Seqinfo) from 2 elements on, also if theirs are identical
  si<-lapply(x,function(g){g@seqinfo})
  si<-si[c(TRUE,!vapply(seq_along(si)[-1],function(i){identical(si[[i]],si[[i-1]])},TRUE))]
  if(length(si)==1 && length(x)>1){
    si<-c(si,si)
  }
  si<-tryCatch({
    while(length(si)>1){
      i<-seq(1,length(si)-1,by=2)
      si<-c(Map(Seqinfo::merge,si[i],si[i+1]),si[-c(i,i+1)])
    }
    si[[1]]
  },error=function(e){NULL})
  if(is.null(si)){
    return(NULL)
  }
  ans<-x[[1]]
  ans@seqnames<-Rle(factor(unlist(lapply(x,function(g){as.character(g@seqnames@values)})),levels=seqlevels(si)),
                    unlist(lapply(x,function(g){g@seqnames@lengths})))
  ans@ranges<-rngs
  ans@strand<-Rle(unlist(lapply(x,function(g){g@strand@values})),unlist(lapply(x,function(g){g@strand@lengths})))
  ans@seqinfo<-si
  mcs<-mc[[1]]
  mcs@listData<-cols
  mcs@nrows<-sum(len)
  ans@elementMetadata<-mcs
  ans
}

# Forked processes for steps that only read their data. mc_job(expr,n_cores) evaluates expr in a
# forked process (parallel::mcparallel) when n_cores>1, here otherwise; mc_value(job) waits for
# it and gives the value of expr, or NULL with wait=FALSE if it isn't done yet. mc_lapply(fs,n_cores)
# is lapply(fs,function(f){f()}), with each f() in a forked process, at most n_cores at a time
# (parallel::mclapply). The warnings of a forked process are given again here (they would be
# lost), and its error stops here.
mc_job<-function(expr,n_cores){
  if(n_cores>1){
    return(parallel::mcparallel(mc_capture(function(){expr})))
  }
  list(value=expr)
}
mc_value<-function(job,wait=TRUE){
  if(!inherits(job,"parallelJob")){
    return(job$value)
  }
  res<-parallel::mccollect(job,wait=wait)
  if(is.null(res)){
    return(NULL)
  }
  mc_release(res[[1]])
}
mc_lapply<-function(fs,n_cores){
  if(n_cores>1){
    return(lapply(parallel::mclapply(fs,mc_capture,mc.cores=n_cores,mc.preschedule=FALSE),mc_release))
  }
  lapply(fs,function(f){f()})
}
mc_capture<-function(f){
  warns<-list()
  value<-withCallingHandlers(f(),warning=function(w){
    warns[[length(warns)+1]]<<-w
    invokeRestart("muffleWarning")
  })
  list(value=value,warnings=warns)
}
mc_release<-function(res){
  if(inherits(res,"try-error")){
    stop(attr(res,"condition"))
  }
  if(is.null(res)){
    stop("A forked process ended without a result (out of memory?)")
  }
  for(w in res$warnings){
    warning(w)
  }
  res$value
}


#' Map transcript coordinates to genomic coordinates
#'
#' This function uses the \code{mapFromTranscripts} function to switch between transcript
#' and genomic coordinates
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param ORFs Set of detected ORFs from the \code{calc_orf_pval} function
#' @param exons exonic regions of the analyzed transcripts, as a GRangesList object
#' @param introns intronic regions of the analyzed transcripts, as a GRangesList object
#' @return exonic coordinates for each ORF.
#' @seealso \code{\link[GenomicFeatures]{mapFromTranscripts}}
#' @export

from_tx_togen<-function(ORFs,exons,introns){
  strand(ORFs)<-rep("*",length(ORFs))
  orfs_gen<-mapFromTranscripts(x = ORFs,transcripts = exons,ignore.strand=F)
  strand(orfs_gen)<-strand(exons[[1]][1])
  # setdiff() of each ORF and the introns, for all ORFs at once: per ORF, setdiff() and [[<- on
  # the GRangesList cost several ms
  list_ma<-GenomicRanges::psetdiff(orfs_gen[1:length(ORFs)],rep(GRangesList(introns),length(ORFs)))
  names(list_ma)<-ORFs$ORF_id_tr
  return(list_ma)
}
