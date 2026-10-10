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

# Step 3: ORFquant() runs detection, quantification and annotation on one genomic
# region; run_ORFquant() runs it on all regions of a sample.


#' Detection, quantification and annotation of translated ORFs in a genomic region
#'
#' This function detects, quantifies and annotates actively translated ORF in a genomic region
#' @details A set of transcripts, together with genome sequence and Ribo-signal are analyzed to extract translated ORFs
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param region GRanges object with genomic coordinates of the genomic region analyzed
#' @param for_ORFquant "for_ORFquant" Robject containing P_sites positions and junction reads
#' @param genetic_code_region GENETIC_CODE table to use
#' @param orf_find.all_starts \code{get_all_starts} parameter for the \code{detect_translated_orfs} function
#' @param orf_find.nostarts \code{nostarts} parameter for the \code{detect_translated_orfs} function
#' @param orf_find.start_sel_cutoff \code{cutoff} parameter for the \code{detect_translated_orfs} function
#' @param orf_find.start_sel_cutoff_ave \code{cutoff_ave} parameter for the \code{detect_translated_orfs} function
#' @param orf_find.cutoff_fr_ave \code{cutoff} parameter for the \code{detect_translated_orfs} function
#' @param orf_quant.cutoff_cums \code{cutoff_cums} parameter for the \code{select_quantify_ORFs} function
#' @param orf_quant.cutoff_pct \code{cutoff_pct} parameter for the \code{select_quantify_ORFs} function
#' @param orf_quant.cutoff_P_sites \code{cutoff_P_sites} parameter for the \code{select_quantify_ORFs} function
#' @param unique_reads Use only signal from uniquely mapping reads? Defaults to \code{FALSE}.
#' @param orf_quant.scaling \code{scaling} parameter for the \code{select_quantify_ORFs} function. Defaults to total_Psites 
#' @return A list containing transcript coordinates, exonic coordinates and annotation for each ORF.\cr\cr
#' The description for each list object is as follows:\cr\cr
#' \code{ORFs_tx}: transcript coordinates of the detected ORFs, with the columns described in \code{\link{ORFquant_output}}.\cr
#' \code{ORFs_gen}: genomic (exon) coordinates of the detected ORFs.\cr
#' \code{ORFs_feat}: list of ORF features together with mapping reads and uniqueness.\cr
#' \code{ORFs_txs_feats}: list of transcript features present in the genomic region, together with mapping reads and uniqueness.\cr
#' \code{ORFs_spl_feat_longest}: splicing annotation for each ORF exon, with respect to the longest annotated coding transcript for each gene.\cr
#' \code{ORFs_spl_feat_maxORF}: splicing annotation for each ORF exon, with respect to the most translated ORF in each gene.\cr
#' \code{selected_txs}: character vector containing the transcript ids of the selected transcripts.\cr
#' \code{ORFs_readthroughs}: (Beta) transcript coordinates of the detected ORFs readthroughs.\cr
#' @seealso \code{\link{select_txs}}, \code{\link{detect_translated_orfs}}, \code{\link{select_quantify_ORFs}}, \code{\link{annotate_ORFs}}, \code{\link{detect_readthrough}}, \code{\link{ORFquant_output}}
#' @export

ORFquant<-function(region,for_ORFquant,genetic_code_region,
                   orf_find.all_starts=T,orf_find.nostarts=F,orf_find.start_sel_cutoff = NA,orf_find.start_sel_cutoff_ave = .5,
                   orf_find.cutoff_fr_ave=.5,orf_quant.cutoff_cums = NA,orf_quant.cutoff_pct = 2,orf_quant.cutoff_P_sites=NA,unique_reads=F,orf_quant.scaling="total_Psites"){
  if(!orf_quant.scaling%in%c("total_Psites","average_coverage")){stop(paste("orf_quant.scaling parameter must be either total_Psites (recommended) or average_coverage"),date())}
  
  P_sites_region<-for_ORFquant$P_sites_all[for_ORFquant$P_sites_all%over%region]
  P_sites_uniq_region<-for_ORFquant$P_sites_uniq[for_ORFquant$P_sites_uniq%over%region]
  P_sites_uniq_mm_region<-for_ORFquant$P_sites_uniq_mm[for_ORFquant$P_sites_uniq_mm%over%region]
  
  res_orfs<-list()
  minimum_reads<-length(P_sites_region)>4
  if(unique_reads){minimum_reads<-length(P_sites_uniq_region)>4}
  
  if(minimum_reads){
    
    selected_transcripts<-select_txs(region = region,P_sites = P_sites_region,P_sites_uniq = P_sites_uniq_region,annotation = GTF_annotation,junction_counts=for_ORFquant$junctions,uniq_signal = unique_reads)
    if(length(selected_transcripts)>0){
      res_orfs<-suppressWarnings(detect_translated_orfs(selected_txs = selected_transcripts,genome_sequence = genome_seq,annotation = GTF_annotation,
                                                        P_sites = P_sites_region,P_sites_uniq = P_sites_uniq_region,P_sites_uniq_mm = P_sites_uniq_mm_region,
                                                        genomic_region=region,genetic_code=genetic_code_region,
                                                        all_starts=orf_find.all_starts,nostarts=orf_find.nostarts,
                                                        start_sel_cutoff=orf_find.start_sel_cutoff,
                                                        start_sel_cutoff_ave=orf_find.start_sel_cutoff_ave,
                                                        cutoff_fr_ave=orf_find.cutoff_fr_ave,uniq_signal = unique_reads))
    }
    if(length(res_orfs)>0){
      res_orfs<-select_quantify_ORFs(results_ORFs=res_orfs,P_sites = P_sites_region,P_sites_uniq = P_sites_uniq_region,
                                     cutoff_cums = orf_quant.cutoff_cums,cutoff_pct = orf_quant.cutoff_pct,cutoff_P_sites=orf_quant.cutoff_P_sites,uniq_signal = unique_reads,scaling = orf_quant.scaling)
    }
    if(length(res_orfs)>0){
      
      res_orfs<-annotate_ORFs(results_ORFs=res_orfs,Annotation=GTF_annotation,genome_sequence = genome_seq,region=region,genetic_code=genetic_code_region)
      res_orfs[["readthrough"]]<-suppressWarnings(detect_readthrough(results_orf = res_orfs,P_sites = P_sites_region,P_sites_uniq = P_sites_uniq_region,P_sites_uniq_mm = P_sites_uniq_mm_region,
                                                                     genome_sequence=genome_seq, annotation = GTF_annotation,genetic_code_table=genetic_code_region,cutoff_fr_ave=orf_find.cutoff_fr_ave,uniq_signal = unique_reads))
      
      
    }
    res_orfs[["genomic_features"]]<-selected_transcripts
    
  }
  
  return(res_orfs)
}

# ORFs_tx as a data.frame: one row per ORF, multiple values joined with ","
ORFs_tx_as_table<-function(ORFs_tx){
  df<-data.frame(seqnames=as.character(seqnames(ORFs_tx)),start=start(ORFs_tx),end=end(ORFs_tx),
                 width=width(ORFs_tx),strand=as.character(strand(ORFs_tx)),stringsAsFactors=F)
  for(n in names(mcols(ORFs_tx))){
    v<-mcols(ORFs_tx)[[n]]
    if(is(v,"GRanges") || is(v,"XStringSet") || is.factor(v)){
      v<-as.character(v)
    }else if(is(v,"List") || is.list(v)){
      v<-vapply(as.list(v),function(x){paste(as.character(x),collapse=",")},"",USE.NAMES=F)
    }
    df[[n]]<-v
  }
  df
}

#' Run the ORFquant pipeline
#'
#' This wrapper function runs the entire ORFquant pipeline
#' @details A set of transcripts, together with genome sequence and Ribo-signal are analyzed to extract translated ORFs
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param annotation_file REQUIRED - path to the *Rannot R file in the annotation directory used in the \code{prepare_annotation_files function}
#' @param for_ORFquant_file REQUIRED - path to the "for_ORFquant" file containing P_sites positions and junction reads
#' @param n_cores REQUIRED - number of cores to use
#' @param prefix prefix to use for the output files. Defaults to same as \code{for_ORFquant_file} (appends to its filename)
#' @param gene_name \code{character} vector of gene names to analyze.
#' @param gene_id \code{character} vector of gene ids to analyze
#' @param genomic_region \code{GRanges} object with genomic regions to analyze
#' @param write_temp_files write temporary files. Defaults to \code{TRUE}
#' @param write_GTF_file write a GTF files with the ORF coordinates. Defaults to \code{TRUE}
#' @param write_protein_fasta write a protein fasta file. Defaults to \code{TRUE}
#' @param interactive should put R object in global environment? Defaults to \code{TRUE}
#' @param stn.orf_find.all_starts \code{orf_find.all_starts} parameter for the \code{ORFquant} function
#' @param stn.orf_find.nostarts \code{orf_find.nostarts} parameter for the \code{ORFquant} function: also find ORFs without a start codon (see \code{Stop_Stop} in \code{\link{get_orfs}}, and \code{\link{ORFquant_output}})? Defaults to \code{FALSE}
#' @param stn.orf_find.start_sel_cutoff \code{orf_find.start_sel_cutoff} parameter for the \code{ORFquant} function
#' @param stn.orf_find.start_sel_cutoff_ave \code{orf_find.start_sel_cutoff_ave} parameter for the \code{ORFquant} functio
#' @param stn.orf_find.cutoff_fr_ave \code{orf_find.cutoff_fr_ave} parameter for the \code{ORFquant} function
#' @param stn.orf_quant.cutoff_cums \code{orf_quant.cutoff_cums} parameter for the \code{ORFquant} function
#' @param stn.orf_quant.cutoff_pct \code{orf_quant.cutoff_pct} parameter for the \code{ORFquant} function
#' @param stn.orf_quant.cutoff_P_sites \code{orf_quant.cutoff_P_sites} parameter for the \code{ORFquant} function
#' @param unique_reads_only Use only signal from uniquely mapping reads? Defaults to \code{FALSE}.
#' @param stn.orf_quant.scaling \code{orf_quant.scaling} parameter for the \code{ORFquant} function. Defaults to total_Psites
#' @param canonical_start_only Use only the canonical start codon (no alternative initiation codons)? Defaults to \code{TRUE}.
#' @param write_TSV_file write a tab-separated file with one row per ORF, the \code{ORFs_tx} table. Defaults to \code{TRUE}
#' @return A set of output files containing transcript coordinates, exonic coordinates and annotation for each ORF, including optional GTF, TSV and protein fasta files.\cr\cr
#' The description for each list object is as follows:\cr\cr
#' \code{tmp_ORFquant_results}: (Optional) RData object file containing the entire set of results for each genomic region.\cr
#' \code{final_ORFquant_results}: RData object file containing the final ORFquant results, see \code{ORFquant}.\cr
#' \code{Protein_sequences.fasta}: (Optional) Fasta file containing the set of translated proteins .\cr
#' \code{Detected_ORFs.gtf}: GTF file with the exons of the selected transcripts (\code{exon} lines) and the exons of the detected ORFs, without the stop codon (\code{CDS} lines). A \code{CDS} line has the \code{transcript_id} of the ORF's transcript, and the \code{ORF_id} (\code{ORF_id_tr}), \code{P_sites}, \code{ORF_pct_P_sites}, \code{ORF_pct_P_sites_pN} and \code{ORFs_pM} of the ORF. A transcript can have more than one ORF, for example a uORF and the annotated CDS, so group the \code{CDS} lines by \code{ORF_id}: tools that take one CDS for each \code{transcript_id}, such as a conversion to genePred or \code{txdbmaker::makeTxDbFromGFF()}, join these ORFs into one CDS or keep only one of them.\cr
#' \code{Detected_ORFs.tsv}: (Optional) Tab-separated file with one row per ORF: the columns of \code{as.data.frame(ORFs_tx)}, transcript coordinates first. Columns with several values per ORF have them separated by commas, and ranges are written as \code{seqname:start-end:strand}. \code{\link{ORFquant_output}} describes its columns.\cr\cr
#' In addition, new columns are added in the ORFs_tx file:\cr\cr
#' \code{ORFs_pM}: number of P_sites for each ORF, divided by ORF length and summing up to a million (akin to TPM).\cr
#' @seealso \code{\link{prepare_annotation_files}}, \code{\link{load_annotation}}, \code{\link{ORFquant}}, \code{\link{ORFquant_output}}
#' @export

run_ORFquant<-function(for_ORFquant_file,annotation_file,n_cores,prefix=for_ORFquant_file,gene_name=NA,gene_id=NA,genomic_region=NA,
                       write_temp_files=T,write_GTF_file=T,write_protein_fasta=T,interactive=T,
                       stn.orf_find.all_starts=T,stn.orf_find.nostarts=F,stn.orf_find.start_sel_cutoff = NA,
                       stn.orf_find.start_sel_cutoff_ave = .5,stn.orf_find.cutoff_fr_ave=.5,
                       stn.orf_quant.cutoff_cums = NA,stn.orf_quant.cutoff_pct = 2,stn.orf_quant.cutoff_P_sites=NA,unique_reads_only=F,canonical_start_only=T,stn.orf_quant.scaling="total_Psites",write_TSV_file=T){    
  
  if(!stn.orf_quant.scaling%in%c("total_Psites","average_coverage")){stop(paste("stn.orf_quant.scaling parameter must be either total_Psites (recommended) or average_coverage! ",date(),sep=""))}
  
  check_n_cores(n_cores)
  
  for (f in c(for_ORFquant_file,annotation_file)){
    if(file.access(f, 0)==-1) {
      stop("The following files don't exist:\n",
           f, "\n")
    }
  }
  
  check_output_dir(prefix)
  
  message("Loading annotation and Ribo-seq signal ... ",date())
  
  load_annotation(annotation_file)
  for_ORFquant_data<-get(load(for_ORFquant_file))
  
  #the P-sites must come from the same genome as the annotation (e.g. not from a for_ORFquant_file made with another annotation)
  lengths_psites<-seqlengths(for_ORFquant_data$P_sites_all)
  if(length(lengths_psites)>0){
    lengths_annot<-seqlengths(GTF_annotation$seqinfo)
    common_chroms<-intersect(names(lengths_psites),names(lengths_annot))
    if(length(common_chroms)==0 || any(lengths_psites[common_chroms]!=lengths_annot[common_chroms],na.rm = TRUE)){
      stop("The P-sites of ",for_ORFquant_file," are not on the genome of the annotation ",annotation_file,": ",
           "their chromosomes (e.g. ",paste(head(names(lengths_psites),3),collapse=", "),") or chromosome lengths are different. ",
           "Run prepare_for_ORFquant with this annotation.")
    }
  }
  if(isTRUE(is.na(GTF_annotation$stop_in_gtf))){
    message("Most annotated CDS regions end before their stop codon (stop_in_gtf is NA): the ORF categories use the codon after the CDS as the annotated stop codon")
  }
  
  genes_red<-reduce(unlist(GTF_annotation$txs_gene))
  
  if(!is.na(gene_name[1])){
    gnid<-unique(GTF_annotation$trann$gene_id[GTF_annotation$trann$gene_name%in%gene_name])
    #the genes' ranges from their transcripts, as genes() makes them, also for genes on several chromosomes
    #or strands (e.g. UCSC's PAR genes), which aren't in GTF_annotation$genes
    genes_red<-genes_red[genes_red%over%unlist(range(GTF_annotation$txs_gene[gnid]))]
  }
  
  if(!is.na(gene_id[1])){
    genes_red<-genes_red[genes_red%over%unlist(range(GTF_annotation$txs_gene[gene_id]))]
  }
  
  if(!is.na(genomic_region[1])){
    genes_red<-genes_red[genes_red%over%genomic_region]
  }
  
  if(length(genes_red)==0){
    stop(paste("Incorrect gene_id | gene_name | genomic_region! ",date(),sep = ""))
  }
  
  message("Loading annotation and Ribo-seq signal --- Done! ",date())
  
  ovs_genesred<-summarizeOverlaps(genes_red,reads = for_ORFquant_data$P_sites_all)
  
  genes_red<-genes_red[which(assay(ovs_genesred)>4)]
  
  if(length(genes_red)==0){
    stop(paste("Not enough P_sites signal over genomic regions! ",date(),sep = ""))
  }
  
  message("Summoning ORFquant with ", sum(for_ORFquant_data$P_sites_all%over%genes_red)," P_sites positions over ", length(genes_red), " genomic regions using ",n_cores," processor(s) ... ",date())
  lengg<-length(genes_red)
  pcts_leng<-as.integer(seq(1,lengg,length.out = 11))
  labs_top<-paste(c(0,seq(10,100,by = 10)),"% completed",sep = "")
  
  #each region's P-sites and junctions, found once: ORFquant() takes them from the genome-wide
  #objects with a %over% per region, 10-30 ms each, so it gets only its region's part (as is, in order)
  regions_idx<-lapply(for_ORFquant_data[c("P_sites_all","P_sites_uniq","P_sites_uniq_mm","junctions")],function(x){
    hts<-findOverlaps(x,genes_red)
    split(queryHits(hts),factor(subjectHits(hts),levels=seq_along(genes_red)))
  })
  
  if(n_cores>1){
    #the regions go to the workers one at a time, in order, each to the next free worker: region costs
    #vary a lot, and equal shares fixed at the start left most workers idle at the end of the loop.
    #the workers are forked once, after run_region is set, so they have it with its data, and a task sends
    #only the region's number (a task sends its function's environment too, unless it is a namespace).
    #outfile="": the workers print the progress lines
    .ORFquant_loop$run_region<-function(g){
      
      if(g%in%pcts_leng){
        message(labs_top[pcts_leng==g])
      }
      
      gen_region<-genes_red[g]
      genetcd<-GTF_annotation$genetic_codes$genetic_code[rownames(GTF_annotation$genetic_codes)==as.character(seqnames(gen_region))]
      genetcd<-getGeneticCode(genetcd)
      if(canonical_start_only){
        attributes(genetcd)$alt_init_codons<-names(which(genetcd=="M"))
      }
      
      for_ORFquant_region<-for_ORFquant_data
      for(n in names(regions_idx)){
        for_ORFquant_region[[n]]<-for_ORFquant_data[[n]][regions_idx[[n]][[g]]]
      }
      
      ORFquant(region=gen_region,for_ORFquant=for_ORFquant_region,genetic_code_region=genetcd,
               orf_find.all_starts=stn.orf_find.all_starts,orf_find.nostarts=stn.orf_find.nostarts,
               orf_find.start_sel_cutoff = stn.orf_find.start_sel_cutoff,orf_find.start_sel_cutoff_ave = stn.orf_find.start_sel_cutoff_ave,
               orf_find.cutoff_fr_ave=stn.orf_find.cutoff_fr_ave,orf_quant.cutoff_cums = stn.orf_quant.cutoff_cums,
               orf_quant.cutoff_pct = stn.orf_quant.cutoff_pct,orf_quant.cutoff_P_sites=stn.orf_quant.cutoff_P_sites,unique_reads = unique_reads_only,orf_quant.scaling = stn.orf_quant.scaling)
      
    }
    cl<-NULL
    ORFs_found<-tryCatch({
      cl<-parallel::makeForkCluster(n_cores,outfile="")
      parallel::clusterApplyLB(cl,seq_along(genes_red),ORFquant_loop_task)
    },finally={
      if(!is.null(cl)) parallel::stopCluster(cl)
      rm("run_region",envir=.ORFquant_loop)
    })
    
    
  }
  if(n_cores==1){
    ORFs_found<-list()
    for(g in 1:length(genes_red)){
      
      if(g%in%pcts_leng){
        message(labs_top[pcts_leng==g])
      }
      
      gen_region<-genes_red[g]
      genetcd<-GTF_annotation$genetic_codes$genetic_code[rownames(GTF_annotation$genetic_codes)==as.character(seqnames(gen_region))]
      genetcd<-getGeneticCode(genetcd)
      if(canonical_start_only){
        attributes(genetcd)$alt_init_codons<-names(which(genetcd=="M"))
      }
      
      for_ORFquant_region<-for_ORFquant_data
      for(n in names(regions_idx)){
        for_ORFquant_region[[n]]<-for_ORFquant_data[[n]][regions_idx[[n]][[g]]]
      }
      
      ORFs_found[[g]]<-ORFquant(region=gen_region,for_ORFquant=for_ORFquant_region,genetic_code_region=genetcd,
                                orf_find.all_starts=stn.orf_find.all_starts,orf_find.nostarts=stn.orf_find.nostarts,
                                orf_find.start_sel_cutoff = stn.orf_find.start_sel_cutoff,orf_find.start_sel_cutoff_ave = stn.orf_find.start_sel_cutoff_ave,
                                orf_find.cutoff_fr_ave=stn.orf_find.cutoff_fr_ave,orf_quant.cutoff_cums = stn.orf_quant.cutoff_cums,
                                orf_quant.cutoff_pct = stn.orf_quant.cutoff_pct,orf_quant.cutoff_P_sites=stn.orf_quant.cutoff_P_sites,unique_reads = unique_reads_only,orf_quant.scaling = stn.orf_quant.scaling)
      
      
    }
    
  }
  message("Summoning ORFquant --- Done! ",date())
  
  message("Exporting ORFquant results ... ",date())
  
  #with more than one core, a forked process saves the tmp file while the tables are built
  tmp_save<-NULL
  if(write_temp_files){
    tmp_save<-mc_job({save(ORFs_found,file=paste(prefix,"tmp_ORFquant_results",sep="_"));TRUE},n_cores)
  }
  
  lens<-elementNROWS(ORFs_found)
  
  ORFs_found<-ORFs_found[lens>0]
  if(length(ORFs_found)==0){mc_value(tmp_save);stop(paste("No ORFs found! Please check sub-codon resolution of Ribo-seq reads and ensure the annotation is correct --- ",date(),"\n"))}
  
  ORFs_found_feats<-ORFs_found
  
  lens<-elementNROWS(ORFs_found)
  
  ORFs_found<-ORFs_found[lens>1]
  if(length(ORFs_found)==0){mc_value(tmp_save);stop(paste("No ORFs found! Please check sub-codon of Ribo-seq reads or that the annotation is correct --- ",date(),"\n"))}
  
  #the tables, each from ORFs_found alone: with more than one core, each one is built in a forked process,
  #at most 3 at a time. Each one can copy GBs of the memory of this process, and more don't make it
  #faster: ORFs_tx alone takes about as long as the other tables on 2 cores
  tables<-mc_lapply(list(ORFs_txs_feats=function(){
  ORFs_txs_feats<-unlist(GRangesList_fast(lapply(ORFs_found_feats,function(x){unlist(x$genomic_features)})))
  ORFs_txs_feats<-ORFs_txs_feats[!duplicated(mcols(ORFs_txs_feats)) | !duplicated(ORFs_txs_feats)]
  },ORFs_tx=function(){

  #ORFs_tx<-unlist(GRangesList(unlist(sapply(ORFs_found,function(x){unlist(x$ORFs_tx_position)}))))
  
  chunks<-seq(1,length(ORFs_found),by = 1000)
  if(length(chunks)==1 || chunks[length(chunks)]<length(ORFs_found)){chunks<-c(chunks,length(ORFs_found))}
  ORFs_tx<-GRangesList()
  for(i in 1:(length(chunks)-1)){
  
    if(i!=(length(chunks)-1)){
    ao<-unlist(GRangesList_fast(unlist(lapply(ORFs_found[chunks[i]:(chunks[i+1]-1)],function(x){unlist(x$ORFs_tx_position)}))))
    ORFs_tx[[i]]<-ao
  }
  if(i==(length(chunks)-1)){
    ORFs_tx[[i]]<-unlist(GRangesList_fast(unlist(lapply(ORFs_found[chunks[i]:(chunks[i+1])],function(x){unlist(x$ORFs_tx_position)}))))
  }
  }

  ORFs_tx<-unlist(ORFs_tx)
  },ORFs_feat=function(){

  ORFs_feat<-unlist(sapply(ORFs_found,function(x){unlist(x$selected_ORFs_features)}))
  ORFs_feat<-GRangesList_fast(sapply(ORFs_feat,function(x){
    #x$X$tx_name<-NULL, on the slots: the $<- calls check the objects, 12 ms per region
    x@elementMetadata@listData$X@elementMetadata@listData$tx_name<-NULL
    return(x)
  }))
  },ORFs_gen=function(){
  ORFs_gen<-unlist(GRangesList_fast(sapply(ORFs_found,function(x){unlist(x$ORFs_genomic_position)})))
  },ORFs_spl_feat_longest=function(){
  ORFs_spl_feat_longest<-unlist(GRangesList_fast(sapply(ORFs_found,function(x){unlist(x$ORFs_splice_feats$annotation_wrt_longest)})))
  },ORFs_spl_feat_maxORF=function(){
  ORFs_spl_feat_maxORF<-unlist(GRangesList_fast(sapply(ORFs_found,function(x){unlist(x$ORFs_splice_feats$annotation_wrt_maxORF)})))
  },ORFs_readthroughs=function(){
  ORFs_readthroughs<-unlist(GRangesList_fast(unlist(sapply(ORFs_found,function(x){unlist(x$readthrough)}))))
  if(length(ORFs_readthroughs)>0){
    ORFs_readthroughs<-ORFs_readthroughs[order(ORFs_readthroughs$P_sites_raw,decreasing = T)]
  }
  ORFs_readthroughs
  }),min(n_cores,3))
  ORFs_txs_feats<-tables$ORFs_txs_feats
  ORFs_tx<-tables$ORFs_tx
  ORFs_feat<-tables$ORFs_feat
  ORFs_gen<-tables$ORFs_gen
  ORFs_spl_feat_longest<-tables$ORFs_spl_feat_longest
  ORFs_spl_feat_maxORF<-tables$ORFs_spl_feat_maxORF
  ORFs_readthroughs<-tables$ORFs_readthroughs
  rm(tables)
  
  selected_txs<-sort(unique(unlist(ORFs_txs_feats$txs_selected)))
  
  
  #toAdd - use the features unique to ORFs as confidence score
  
  na_ps<-is.na(ORFs_tx$P_sites)
  ORFs_tx$P_sites_pN<-NA
  ORFs_tx$P_sites_pN[!na_ps]<-ORFs_tx$P_sites[!na_ps]/(width(ORFs_tx)[!na_ps])
  ORFs_tx$ORFs_pM<-NA
  ORFs_tx$ORFs_pM[!na_ps]<-ORFs_tx$P_sites_pN[!na_ps]*(1000000/(sum(ORFs_tx$P_sites_pN[!na_ps])))
  ORFs_tx$P_sites_pN<-NULL
  
  
  ORFquant_results<-list(ORFs_tx,ORFs_gen,ORFs_feat,ORFs_spl_feat_longest,ORFs_spl_feat_maxORF,ORFs_readthroughs,ORFs_txs_feats,selected_txs)
  names(ORFquant_results)<-c("ORFs_tx","ORFs_gen","ORFs_feat","ORFs_spl_feat_longest","ORFs_spl_feat_maxORF","ORFs_readthroughs","ORFs_txs_feats","selected_txs")
  
  #added for seqinfo problem
  
  
  x<-ORFquant_results$ORFs_readthroughs
  seqf<-Seqinfo(seqnames = names(GTF_annotation$exons_txs),seqlengths = sum(width(GTF_annotation$exons_txs)),isCircular = NA,genome = NA)
  x@seqnames<-Rle(factor(as.character(x@seqnames),levels = seqlevels(seqf)))
  x@seqinfo<-seqf
  ORFquant_results$ORFs_readthroughs<-x
  
  x<-ORFquant_results$ORFs_tx
  seqf<-Seqinfo(seqnames = names(GTF_annotation$exons_txs),seqlengths = sum(width(GTF_annotation$exons_txs)),isCircular = NA,genome = NA)
  x@seqnames<-Rle(factor(as.character(x@seqnames),levels = seqlevels(seqf)))
  x@seqinfo<-seqf
  aa<-suppressWarnings(GRanges(x$longest_ORF))
  aa@seqnames<-Rle(factor(as.character(x$longest_ORF@seqnames),levels = seqlevels(seqf)))
  aa@seqinfo<-seqf
  mcols(x)$longest_ORF<-NULL
  x$longest_ORF<-aa
  
  
  ORFquant_results$ORFs_tx<-x
  
  ORFquant_results$annotation_file<-annotation_file
  ORFquant_results$psite_data_file <- for_ORFquant_file
  
  #with more than one core, a forked process saves the results while the files below are written
  final_save<-mc_job({save(ORFquant_results ,file = paste(prefix,"final_ORFquant_results",sep="_"));TRUE},n_cores)
  
  if(write_TSV_file){
    write.table(ORFs_tx_as_table(ORFquant_results$ORFs_tx),file=paste(prefix,"Detected_ORFs.tsv",sep="_"),sep="\t",quote=F,row.names=F)
  }
  
  if(write_GTF_file){
    map_tx_genes<-mcols(ORFs_tx)[,c("ORF_id_tr","gene_id","gene_biotype","gene_name","transcript_id","transcript_biotype","P_sites","ORF_pct_P_sites","ORF_pct_P_sites_pN","ORFs_pM")]
    
    match_ORF<-match(names(ORFs_gen),map_tx_genes$ORF_id_tr)
    
    ORFs_gen$transcript_id<-map_tx_genes[match_ORF,"transcript_id"]
    
    match_tx<-match(ORFs_gen$transcript_id,map_tx_genes$transcript_id)
    
    ORFs_gen$transcript_biotype<-map_tx_genes[match_tx,"transcript_biotype"]
    ORFs_gen$gene_id<-map_tx_genes[match_tx,"gene_id"]
    ORFs_gen$gene_biotype<-map_tx_genes[match_tx,"gene_biotype"]
    ORFs_gen$gene_name<-map_tx_genes[match_tx,"gene_name"]
    ORFs_gen$ORF_id<-map_tx_genes[match_ORF,"ORF_id_tr"]
    
    ORFs_gen$P_sites<-round(map_tx_genes[match_ORF,"P_sites"],digits=4)
    ORFs_gen$ORF_pct_P_sites<-round(map_tx_genes[match_ORF,"ORF_pct_P_sites"],digits=4)
    ORFs_gen$ORF_pct_P_sites_pN<-round(map_tx_genes[match_ORF,"ORF_pct_P_sites_pN"],digits=4)
    ORFs_gen$ORFs_pM<-round(map_tx_genes[match_ORF,"ORFs_pM"],digits=4)
    ORFs_readthroughs<-ORFquant_results$ORFs_readthroughs
    if(is(ORFs_readthroughs$Protein,"list")){
      proteins_readthrough<-AAStringSet(lapply(ORFs_readthroughs$Protein,"[[",1))
    }else{proteins_readthrough<-AAStringSet(ORFs_readthroughs$Protein)}
    if(length(proteins_readthrough)>0){
      names(proteins_readthrough)<-paste(ORFs_readthroughs$ORF_id_tr,ORFs_readthroughs$gene_biotype,ORFs_readthroughs$gene_id,"readthrough","readthrough",sep="|")
      proteins_readthrough<-narrow(proteins_readthrough,start = start(proteins_readthrough)[1]+1)
      proteins_readthrough<-AAStringSet(gsub(proteins_readthrough,pattern = "[*]",replacement = "X"))
    }
    
    proteins<-AAStringSet(ORFs_tx$Protein)
    names(proteins)<-paste(ORFs_tx$ORF_id_tr,ORFs_tx$gene_biotype,ORFs_tx$gene_id,ORFs_tx$ORF_category_Gen,ORFs_tx$ORF_category_Tx_compatible,sep="|")
    proteins<-c(proteins,proteins_readthrough)
    if(write_protein_fasta){
      writeXStringSet(proteins,filepath= paste(prefix,"Protein_sequences.fasta",sep="_"))
    }
    
    map_tx_genes<-GTF_annotation$trann
    #ORFs_gen$transcript_id<-names(ORFs_gen)
    ORFs_gen$type="CDS"
    exs_gtf<-unlist(GTF_annotation$exons_txs[selected_txs])
    mcols(exs_gtf)<-NULL
    exs_gtf$transcript_id<-names(exs_gtf)
    exs_gtf$transcript_biotype<-map_tx_genes[match(exs_gtf$transcript_id,map_tx_genes$transcript_id),"transcript_biotype"]
    exs_gtf$gene_id<-map_tx_genes[match(exs_gtf$transcript_id,map_tx_genes$transcript_id),"gene_id"]
    exs_gtf$gene_biotype<-map_tx_genes[match(exs_gtf$transcript_id,map_tx_genes$transcript_id),"gene_biotype"]
    exs_gtf$gene_name<-map_tx_genes[match(exs_gtf$transcript_id,map_tx_genes$transcript_id),"gene_name"]
    mcols(exs_gtf)[,names(mcols(ORFs_gen))[!names(mcols(ORFs_gen))%in%names(mcols(exs_gtf))]]<-NA
    exs_gtf$type<-"exon"
    mcols(ORFs_gen)<-mcols(ORFs_gen)[,names(mcols(exs_gtf))]
    all<-sort(c(exs_gtf,ORFs_gen))
    all$`source`="ORFquant"
    names(all)<-NULL
    suppressWarnings(export.gff2(object=all,con=paste(prefix,"Detected_ORFs.gtf",sep="_")))
  }
  mc_value(final_save)
  mc_value(tmp_save)
  
  if(interactive){
    for_ORFquant<<-for_ORFquant_data
    ORFquant_results<<-ORFquant_results
  }
  message("Exporting ORFquant results --- Done! ",date())
  invisible(ORFquant_results)
}

#the task of run_ORFquant's workers (n_cores>1): its environment is the namespace, so a task doesn't send
#run_region's data, and the workers find run_region in their copy of .ORFquant_loop
.ORFquant_loop<-new.env()
ORFquant_loop_task<-function(g) .ORFquant_loop$run_region(g)
