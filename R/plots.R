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

# The plots of the results.


#' Plot general statistics about ORFquant results
#'
#' This function produces a series of plots and statistics about the set ORFs called by ORFquant compared to the annotation.
#' IMPORTANT: Use only on transcriptome-wide ORFquant results. See \code{run_ORFquant}
#' @keywords ORFquant, Ribo-seQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param for_ORFquant_file path to the "for_ORFquant" file containing P_sites positions and junction reads
#' @param ORFquant_output_file Full path to the "_final_ORFquant_results" RData object output by ORFquant. See \code{run_ORFquant}
#' @param annotation_file Full path to the *Rannot R file in the annotation directory used in the \code{prepare_annotation_files function}
#' @param coverage_file_plus Optional. Full path to a Ribo-seq coverage (no P-sites but read coverage) bigwig file (plus strand), as the ones created by \code{RiboseQC}
#' @param coverage_file_minus Optional. Full path to a Ribo-seq coverage (no P-sites but read coverage) bigwig file (minus strand), as the ones created by \code{RiboseQC}
#' @param output_plots_path Full path to the directory where plots in .pdf format are stored.
#' @param prefix prefix appended to output filenames
#' @return the function exports a RData object (*ORFquant_plots_RData) containing data to produce all plots, and produces different QC plots in .pdf format.
#' The plots created are as follows:\cr\cr
#' \code{ORFs_found}: Number of ORF categories detected per gene biotype.\cr
#' \code{ORFs_found_pct_tr}: Distribution of ORF_pct_P_sites (% of gene translation) for different ORF categories and gene biotypes.\cr
#' \code{ORFs_found_ORFs_pM}: Distribution of ORFs_pM (ORFs per Million, similar to TPM) for different ORF categories and gene biotypes.\cr
#' \code{ORFs_found_len}: Distribution of ORF length for different ORF categories and gene biotypes.\cr
#' \code{ORFs_genes}: Number of detected ORFs per gene.\cr
#' \code{ORFs_genes_tpm}: Gene level TPM values, plotted by number of ORFs detected.\cr
#' \code{ORFs_maxiso}:  Number of genes plotted against the percentages of gene translation of their most translated ORF.\cr
#' \code{ORFs_maxiso_tpm}: Gene level TPM values, plotted against the percentages of gene translation of their most translated ORF.\cr
#' \code{Sel_txs_genes}: Number of genes plotted against the number of selected transcripts.\cr
#' \code{Sel_txs_genes_tpm}: Gene level TPM values, plotted against the number of selected transcripts.\cr
#' \code{Sel_txs_genes_pct}: Percentages of annotated trascripts per gene, plotted against the number of selected transcripts.\cr
#' \code{Sel_txs_bins_juns}: Percentages of covered exonic bins or junctions, using all annotated transcripts, coding transcripts only, or the set of selected transcripts.\cr
#' \code{Meta_splicing_coverage}: Aggregate signal of Ribo-seq coverage and normalized ORF coverage across different splice sites combinations, with different mixtures of translated overlapping ORFs.
#' @seealso \code{\link{run_ORFquant}}
#' @export



plot_ORFquant_results<-function(for_ORFquant_file,ORFquant_output_file,annotation_file,coverage_file_plus=NA,coverage_file_minus=NA,output_plots_path=NA,prefix=NA){
  
  
  message("Plotting ORFquant results for ",ORFquant_output_file," ... ",date())
  
  if(is.na(prefix)){prefix<-gsub(gsub(basename(ORFquant_output_file),pattern = "_final_ORFquant_results",replacement = ""),pattern = "final_ORFquant_results",replacement = "")}
  
  for (f in c(for_ORFquant_file,ORFquant_output_file,annotation_file)){
    if(file.access(f, 0)==-1) {
      stop("The following files don't exist:\n",
           f, "\n")
    }
  }
  
  if(is.na(output_plots_path)){output_plots_path<-paste0(ORFquant_output_file, "_plots")}
  dir.create(output_plots_path,recursive=TRUE, showWarnings=FALSE)
  
  load(annotation_file)
  load(ORFquant_output_file)
  ORFs_tx<-ORFquant_results$ORFs_tx
  ORFs_gen<-ORFquant_results$ORFs_gen
  selected_txs<-ORFquant_results$selected_txs
  ORFs_txs_feats<-ORFquant_results$ORFs_txs_feats
  ORFs_feat<-ORFquant_results$ORFs_feat
  
  for_ORFquant<-get(load(for_ORFquant_file))
  
  list_ORFquant_plots<-list()
  
  
  gids<-ORFs_tx$gene_id
  cat_biot<-ORFs_tx$gene_biotype
  cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
  cat_biot[!cat_biot%in%c("protein_coding","pseudogene")]<-"non-coding RNA"
  tbid<-table(gids)
  tbid[tbid>3]<-">3"
  #boxplot(log(tpms)~tbid)
  tb_bio<-cat_biot[match(names(tbid),gids)]
  
  
  cat_tx<-ORFs_tx$ORF_category_Tx_compatible
  cat_tx[cat_tx=="novel"]<-"not_annotated"
  cat_tx[cat_tx=="overl_uORF"]<-"uORF"
  cat_tx[cat_tx=="overl_dORF"]<-"other"
  cat_tx[cat_tx=="dORF"]<-"other"
  
  cat_gen<-ORFs_tx$ORF_category_Gen
  cat_biot<-ORFs_tx$compatible_biotype
  cat_biot_gene<-ORFs_tx$gene_biotype
  cat_biot[cat_biot!="protein_coding" & cat_biot_gene=="protein_coding"]<-"non-coding isoform"
  cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
  cat_biot[!cat_biot%in%c("protein_coding","pseudogene","non-coding isoform")]<-"non-coding RNA"
  
  #cat_tx[intersect(grep(ORFs_tx$gene_biotype,pattern = "pseudogene"),which(cat_tx=="novel"))]<-"novel_pseudogene"
  cat_tx[cat_tx%in%c("C_extension","C_truncation","NC_extension","N_extension","nested_ORF","overl_dORF")]<-"other"
  
  cat_gen[grep(cat_gen,pattern = "novel")]<-"non-overlapping\nCDS regions"
  cat_gen[grep(cat_gen,pattern = "non-overlapping\nCDS regions",invert = T)]<-"overlapping\nCDS regions"
  levvs<-c("ORF_annotated","N_truncation","uORF","other","not_annotated")
  
  df<-melt(table(cat_tx,cat_biot))
  
  df$value[df$value==0]<-NA
  df$cat_tx<-factor(df$cat_tx,levels=levvs)
  df$cat_biot<-factor(df$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA"))
  
  
  colli<-c("red","orange","dark blue","cornflowerblue")
  
  
  dfiso<-data.frame(cat_tx,cat_biot,ORFs_tx$ORF_pct_P_sites)
  dfiso$cat_tx<-factor(dfiso$cat_tx,levels=levvs)
  dfiso$cat_biot<-factor(dfiso$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA"))
  
  dfpnpm<-data.frame(cat_tx,cat_biot,ORFs_tx$ORFs_pM)
  #dfpnpm$value[dfpnpm$value==0]<-NA
  dfpnpm$cat_tx<-factor(dfpnpm$cat_tx,levels=levvs)
  dfpnpm$cat_biot<-factor(dfpnpm$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA"))
  
  dfwid<-data.frame(cat_tx,cat_biot,width(ORFs_tx))
  dfwid$cat_tx<-factor(dfwid$cat_tx,levels=levvs)
  dfwid$cat_biot<-factor(dfwid$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA"))
  
  df<-df[!is.na(df$value),]
  
  
  a<-ggplot(df,aes(x=cat_biot,y=value,fill=cat_biot))
  a<-a + geom_bar(stat="identity",position = "dodge",colour="black")
  a<-a + facet_grid(. ~ cat_tx ,drop = T,scales = "free_x")
  a<-a + ylim(0,max(df$value,na.rm=T)*1.1)
  a<-a + geom_text(aes(x=cat_biot, y=value, hjust=.5,vjust=-.2,label=value),check_overlap = TRUE,colour="black",position = position_dodge(width = .9),size=3.5)
  a<-a + theme_bw()
  a<-a + ylab("n of ORFs")
  a<-a + xlab("")
  a<-a + scale_fill_manual(values = colli,"biotype")
  a<-a + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
  a<-a + theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=45, vjust=0.5, size=15))
  orfs_found_n<-a + theme(strip.text.x = element_text(size=14, face="bold"),strip.text.y = element_text(size=9),strip.background = element_rect(colour="black", fill=c(rep("darkkhaki",8),rep("red",8))))
  
  list_ORFquant_plots[["ORFs_found"]]<-orfs_found_n
  list_ORFquant_plots[["ORFs_found"]][["pars"]]<-c(14,4)
  
  rownames(dfiso)<-NULL
  b<-ggplot(dfiso,aes(x=cat_biot,y=ORFs_tx.ORF_pct_P_sites,fill=cat_biot))
  b<-b + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
  b<-b + theme_bw()
  b<-b + ylab("ORF_pct_P-sites")
  b<-b + xlab("")
  b<-b + facet_grid(. ~  cat_tx ,drop = T,scales = "free_x")
  #b<-b + theme(legend.position="none")
  b<-b + scale_fill_manual(values = colli,"biotype")
  b<-b + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
  b<-b + theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=45, vjust=0.5, size=15))
  orfs_found_pct<-b + theme(strip.text.x = element_text(size=14, face="bold"),strip.text.y = element_text(size=9),strip.background = element_rect(colour="black", fill=c(rep("darkkhaki",8),rep("red",8))))
  
  list_ORFquant_plots[["ORFs_found_pct_tr"]]<-orfs_found_pct
  list_ORFquant_plots[["ORFs_found_pct_tr"]][["pars"]]<-c(14,4)
  
  
  c<-ggplot(dfpnpm,aes(x=cat_biot,y=ORFs_tx.ORFs_pM+1,fill=cat_biot))
  c<-c + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
  #c<-c + geom_jitter(aes(x=cat_tx,y=ORFs_tx.ORFs_pM),position = position_jitterdodge(jitter.width = .5), alpha = 0.2)
  c<-c + scale_y_log10(breaks=c(1,11,101,1001),limits=c(1,max(dfpnpm$ORFs_tx.ORFs_pM)*1.1),labels=c(1,11,101,1001)-1)
  #c<-c + geom_text(aes(x=cat_tx, y=ORFs_tx.ORFs_pM, hjust="top",label=ORFs_tx.ORFs_pM),colour="black",position = position_dodge(width = 1),size=5)
  c<-c + facet_grid(. ~ cat_tx ,drop = T,scales = "free_x")
  c<-c + theme_bw()
  c<-c + xlab("")
  c<-c + ylab("P-sites_pNpM")
  c<-c + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
  #c<-c + theme(legend.position="none")
  c<-c + scale_fill_manual(values = colli,"biotype")
  c<-c + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
  c<-c + theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=45, vjust=0.5, size=16))
  c<-c + theme(strip.text.x = element_text(size=16, face="bold"),strip.text.y = element_text(size=20),strip.background = element_rect(colour="black", fill=c("darkkhaki")))
  orfs_found_pspn<-c + theme(strip.text.x = element_text(size=14, face="bold"),strip.text.y = element_text(size=9),strip.background = element_rect(colour="black", fill=c(rep("darkkhaki",8),rep("red",8))))
  
  list_ORFquant_plots[["ORFs_found_ORFs_pM"]]<-orfs_found_pspn
  list_ORFquant_plots[["ORFs_found_ORFs_pM"]][["pars"]]<-c(14,4)
  
  d<-ggplot(dfwid,aes(x=cat_biot,y=width.ORFs_tx.,fill=cat_biot))
  d<-d + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
  #c<-c + geom_jitter(aes(x=cat_tx,y=ORFs_tx.ORFs_pM),position = position_jitterdodge(jitter.width = .5), alpha = 0.2)
  d<-d + scale_y_log10(breaks=c(10,100,1000,10000),labels=c(10,100,1000,10000),limits=c(min(dfwid$width.ORFs_tx.),max(dfwid$width.ORFs_tx.)*1.1))
  #c<-c + geom_text(aes(x=cat_tx, y=ORFs_tx.ORFs_pM, hjust="top",label=ORFs_tx.ORFs_pM),colour="black",position = position_dodge(width = 1),size=5)
  d<-d + theme_bw()
  d<-d + facet_grid(. ~ cat_tx ,drop = T,scales = "free_x")
  d<-d + xlab("")
  d<-d + ylab("ORF length (nt)")
  d<-d + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
  #c<-c + theme(legend.position="none")
  d<-d + scale_fill_manual(values = colli,"biotype")
  d<-d + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
  d<-d + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=16))
  orfs_found_len<-d + theme(strip.text.x = element_text(size=14, face="bold"),strip.text.y = element_text(size=9),strip.background = element_rect(colour="black", fill=c(rep("darkkhaki",8),rep("red",8))))
  
  list_ORFquant_plots[["ORFs_found_len"]]<-orfs_found_len
  list_ORFquant_plots[["ORFs_found_len"]][["pars"]]<-c(14,4)
  
  
  red_ex <- unlist(GTF_annotation$exons_txs)
  red_ex$gene_id<-GTF_annotation$trann$gene_id[match(names(red_ex),GTF_annotation$trann$transcript_id)]
  red_ex<-reduce(split(red_ex,red_ex$gene_id))
  
  ovv<-findOverlaps(red_ex,for_ORFquant$P_sites_all)
  aggo<-aggregate(for_ORFquant$P_sites_all$score[ovv@to],by=list(ovv@from),sum)
  cnts<-aggo[,2]
  names(cnts)<-names(red_ex)[aggo[,1]]
  cnts0<-rep(0,sum(!names(red_ex)%in%names(cnts)))
  names(cnts0)<-names(red_ex)[!names(red_ex)%in%names(cnts)]
  cnts<-c(cnts,cnts0)[names(red_ex)]
  
  cnts<-DataFrame(gene_id=names(cnts),P_sites=as.numeric(cnts))
  cnts<-cnts[match(cnts$gene_id,names(red_ex)),]
  cnts$RPKM<-cnts$P_sites/(sum(cnts$P_sites)/1e06)
  cnts$RPKM<-round(cnts$RPKM/(sum(width(red_ex))/1000),digits = 4)
  cnts$TPM<-cnts$P_sites/(sum(width(red_ex))/1000)
  cnts$TPM<-round(cnts$TPM/(sum(cnts$TPM)/1e06),digits = 4)
  
  
  tpms<-cnts[match(names(tbid),cnts$gene_id),"TPM"]
  gids<-ORFs_tx$gene_id
  cat_biot<-ORFs_tx$gene_biotype
  cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
  cat_biot[!cat_biot%in%c("protein_coding","pseudogene")]<-"non-coding RNA"
  tbid<-table(gids)
  tbid[tbid>3]<-">3"
  #boxplot(log(tpms)~tbid)
  tb_bio<-cat_biot[match(names(tbid),gids)]
  max_ORF<-aggregate(ORFs_tx$ORF_pct_P_sites,by=list(ORFs_tx$gene_id),max)
  colnames(max_ORF)<-c("gene_id","ORF_pct_P_sites")
  
  multg<-names(which(tbid!="1"))
  #mult_max_ORF<-max_ORF[max_ORF$gene_id%in%multg,]
  mult_max_ORF<-max_ORF
  
  maxiso<-cut(mult_max_ORF$ORF_pct_P_sites,breaks = seq(0,100,by = 10),include.lowest = T,right = T)
  maxiso[is.na(maxiso)]<-"(90,100]"
  maxiso<-gsub(maxiso,pattern = ",",replacement = "-")
  maxiso<-gsub(maxiso,pattern = "\\[",replacement = "")
  maxiso<-gsub(maxiso,pattern = "]",replacement = "")
  maxiso<-gsub(maxiso,pattern = "\\(",replacement = "")
  tpms_mult<-tpms[mult_max_ORF$gene_id]
  #boxplot(log(tpms_mult+1)~maxiso)
  df<-data.frame(tbid,tb_bio,tpms)
  
  df$Freq<-factor(df$Freq,levels=c("1","2","3",">3"))
  
  df$tb_bio<-factor(df$tb_bio,levels=c("protein_coding","pseudogene","non-coding RNA"))
  df<-table(df$Freq)
  df<-melt(df)
  
  df$Var1<-factor(df$Var1,levels=c("1","2","3",">3"))
  
  colli<-c("red","dark blue","cornflowerblue")
  a<-ggplot(df,aes(x=Var1,y=value,fill="dark grey"))
  a<-a + geom_bar(stat="identity",position = "stack",colour="black")
  #a<-a + scale_y_log10(breaks=c(1,10,100,1000,10000),limits=c(1,15000))
  a<-a + geom_text(aes(x=Var1, y=value, hjust=.5,vjust=-.8,label=value),check_overlap = TRUE,colour="black",size=4)
  a<-a + theme_bw()
  a<-a + theme(legend.position="none")
  a<-a  + ylim(0,max(df$value)*1.1)
  a<-a + ylab("n of genes")
  a<-a + xlab("n of ORFs")
  a<-a + scale_fill_manual(values = "dark grey","biotype")
  a<-a + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))
  orfsn_genes<-a + theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(angle=45, vjust=0.5, size=18))
  
  list_ORFquant_plots[["ORFs_genes"]]<-orfsn_genes
  list_ORFquant_plots[["ORFs_genes"]][["pars"]]<-c(7,6)
  
  
  df<-data.frame(tbid,tb_bio,tpms)
  df$tb_bio<-factor(df$tb_bio,levels=names(sort(table(df$tb_bio),decreasing = T)))
  df$Freq<-factor(df$Freq,levels=c("1","2","3",">3"))
  #df<-df[df$Freq!="1",]
  rownames(df)<-NULL
  df$gids<-NULL
  b<-ggplot(df,aes(x=Freq,y=tpms+1,fill="dark grey"))
  b<-b + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
  b<-b + scale_y_log10(breaks=c(1,11,101,1001),limits=c(1,max(df$tpms+1)*1.1),labels=c(0,10,100,1000))
  #b<-b + geom_text(aes(x=cat_tx, y=value, hjust="top",label=value),colour="black",position = position_dodge(width = 1),size=5)
  b<-b + theme_bw()
  b<-b + xlab("n of ORFs")
  b<-b + theme(legend.position="none")
  
  b<-b + ylab("TPM (gene level)")
  b<-b + scale_fill_manual(values ="dark grey")
  b<-b + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))
  orfsn_genes_tpm<-b + theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(angle=45, vjust=0.5, size=20))
  
  list_ORFquant_plots[["ORFs_genes_tpm"]]<-orfsn_genes_tpm
  list_ORFquant_plots[["ORFs_genes_tpm"]][["pars"]]<-c(7,6)
  
  
  df<-melt(table(maxiso))
  df$maxiso<-factor(df$maxiso,levels = rev(df$maxiso))
  c<-ggplot(df,aes(x=maxiso,y=value,fill="dark grey"))
  c<-c + geom_bar(stat="identity",position = "stack",colour="black")
  #a<-a + scale_y_log10(breaks=c(1,10,100,1000,10000),limits=c(1,15000))
  c<-c + geom_text(aes(x=maxiso, y=value, hjust=.5,vjust=-.8,label=value),check_overlap = TRUE,colour="black",size=4)
  c<-c + theme_bw()
  c<-c + theme(legend.position="none")
  c<-c + ylim(0,max(df$value)*1.1)
  c<-c + ylab("n of genes")
  c<-c + xlab("ORF_pct_P_sites (max ORF)")
  c<-c + scale_fill_manual(values = "dark grey")
  c<-c + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))
  orfsn_genes_maxiso<-c + theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(angle=45, vjust=0.5, size=20))
  
  list_ORFquant_plots[["ORFs_maxiso"]]<-orfsn_genes_maxiso
  list_ORFquant_plots[["ORFs_maxiso"]][["pars"]]<-c(7,6)
  
  
  df<-data.frame(tpms_mult,maxiso)
  df$maxiso<-factor(df$maxiso,levels = rev(levels(df$maxiso)))
  rownames(df)<-NULL
  d<-ggplot(df,aes(x=maxiso,y=tpms_mult+1,fill="dark grey"))
  d<-d + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
  d<-d + scale_y_log10(breaks=c(1,11,101,1001),limits=c(1,max(tpms_mult+1)*1.1),labels=c(0,10,100,1000))
  #b<-b + geom_text(aes(x=cat_tx, y=value, hjust="top",label=value),colour="black",position = position_dodge(width = 1),size=5)
  d<-d + theme_bw()
  d<-d + theme(legend.position="none")
  
  d<-d + xlab("ORF_pct_P_sites (max ORF)")
  d<-d + ylab("TPM (gene level)")
  d<-d + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
  d<-d + scale_fill_manual(values ="dark grey")
  
  d<-d + theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(angle=45, vjust=0.5, size=20))
  orfsn_genes_maxiso_tpm<-d + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))
  
  list_ORFquant_plots[["ORFs_maxiso_tpm"]]<-orfsn_genes_maxiso_tpm
  list_ORFquant_plots[["ORFs_maxiso_tpm"]][["pars"]]<-c(7,6)
  
  
  
  ch_txs_sel<-CharacterList(split(selected_txs,GTF_annotation$trann$gene_id[match(selected_txs,GTF_annotation$trann$transcript_id)]))
  tot_n_tx<-table(GTF_annotation$trann$gene_id)[names(ch_txs_sel)]
  genes_m<-names(which(tot_n_tx>1))
  genes_m<-genes_m[genes_m%in%ORFs_tx$gene_id]
  #genes_m<-unique(ORFs_tx$gene_id)
  tot_n_tx<-tot_n_tx[genes_m]
  ch_txs_sel<-ch_txs_sel[genes_m]
  
  tpms<-cnts[match(names(ch_txs_sel),cnts$gene_id),"TPM"]
  pct_sel<-(elementNROWS(ch_txs_sel))/tot_n_tx*100
  nsels<-elementNROWS(ch_txs_sel)
  qnt<-cut(nsels,breaks = unique(c(0,3,6,9,max(nsels))),include.lowest = T)
  qnt<-gsub(qnt,pattern = ",",replacement = "-")
  qnt<-gsub(qnt,pattern = "\\[",replacement = "")
  qnt<-gsub(qnt,pattern = "]",replacement = "")
  qnt<-gsub(qnt,pattern = "\\(",replacement = "")
  qnt[qnt=="0-3"]<-"1-3"
  qnt[qnt=="3-6"]<-"4-6"
  qnt[qnt=="6-9"]<-"7-9"
  qnt[qnt=="9-49"]<-"10-49"
  qnt<-factor(qnt,levels = names(sort(table(qnt),decreasing = T)))
  df<-melt(table(qnt))
  
  a<-ggplot(df,aes(x=qnt,y=value,fill="dark grey"))
  a<-a + geom_bar(stat="identity",position = "stack",colour="black")
  #a<-a + scale_y_log10(breaks=c(1,10,100,1000,10000),limits=c(1,15000))
  a<-a + geom_text(aes(x=qnt, y=value, hjust=.5,vjust=-.8,label=value),check_overlap = TRUE,colour="black",size=4)
  a<-a + theme_bw()
  a<- a + ylim(0,max(df$value)*1.1)
  a<-a + theme(legend.position="none")
  a<-a + ylab("n of genes")
  a<-a + xlab("selected Txs per gene")
  a<-a + scale_fill_manual(values = "dark grey","biotype")
  a<-a + theme(axis.title.y = element_text(size=16),axis.text.y  = element_text(angle=45, vjust=0.5, size=16))
  sel_txs_genes<-a + theme(axis.title.x = element_text(size=16),axis.text.x  = element_text(angle=45, vjust=0.5, size=16))
  
  list_ORFquant_plots[["Sel_txs_genes"]]<-sel_txs_genes
  list_ORFquant_plots[["Sel_txs_genes"]][["pars"]]<-c(7,6)
  
  df<-data.frame(tpm=tpms[names(ch_txs_sel)],elementNROWS(ch_txs_sel))
  df$qnt<-qnt
  rownames(df)<-NULL
  b<-ggplot(df,aes(x=qnt,y=tpms+1,fill="dark grey"))
  b<-b + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
  b<-b + scale_y_log10(breaks=c(1,11,101,1001,10001),limits=c(1,max(df$tpm+1)*1.1),labels=c(0,10,100,1000,10000))
  #b<-b + geom_text(aes(x=cat_tx, y=value, hjust="top",label=value),colour="black",position = position_dodge(width = 1),size=5)
  b<-b + theme_bw()
  b<-b + xlab("selected Txs per gene")
  b<-b + theme(legend.position="none")
  b<-b + ylab("TPM (gene)")
  b<-b + scale_fill_manual(values ="dark grey")
  b<-b + theme(axis.title.y = element_text(size=16),axis.text.y  = element_text(angle=45, vjust=0.5, size=16))
  sel_txs_genes_tpm<-b + theme(axis.title.x = element_text(size=16),axis.text.x  = element_text(angle=45, vjust=0.5, size=16))
  
  list_ORFquant_plots[["Sel_txs_genes_tpm"]]<-sel_txs_genes_tpm
  list_ORFquant_plots[["Sel_txs_genes_tpm"]][["pars"]]<-c(7,6)
  
  
  df<-data.frame(pct_sel=as.numeric(pct_sel),elementNROWS(ch_txs_sel))
  df$qnt<-qnt
  rownames(df)<-NULL
  c<-ggplot(df,aes(x=qnt,y=pct_sel,fill="dark grey"))
  c<-c + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
  #b<-b + scale_y_log10(breaks=c(1,11,101,1001,10001),limits=c(1,10001),labels=c(0,10,100,1000,10000))
  #b<-b + geom_text(aes(x=cat_tx, y=value, hjust="top",label=value),colour="black",position = position_dodge(width = 1),size=5)
  c<-c + theme_bw()
  c<-c + xlab("selected Txs per gene")
  c<-c + theme(legend.position="none")
  c<-c + ylab("% annotated Txs")
  c<-c + scale_fill_manual(values ="dark grey")
  c<-c + theme(axis.title.y = element_text(size=16),axis.text.y  = element_text(angle=45, vjust=0.5, size=16))
  sel_txs_genes_pct<-c + theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=45, vjust=0.5, size=16))
  
  list_ORFquant_plots[["Sel_txs_genes_pct"]]<-sel_txs_genes_pct
  list_ORFquant_plots[["Sel_txs_genes_pct"]][["pars"]]<-c(7,6)
  
  gens_sel<-unique(GTF_annotation$trann$gene_id[GTF_annotation$trann$transcript_id%in%selected_txs])
  tx_all<-GTF_annotation$trann$transcript_id[GTF_annotation$trann$gene_id%in%gens_sel]
  tx_sel<-selected_txs
  if(sum(!tx_sel%in%names(GTF_annotation$exons_txs))>0){warning(paste("some selected Txs not present in annotation!, like:",head(tx_sel[!tx_sel%in%names(GTF_annotation$exons_txs)],1)))}
  tx_sel<-tx_sel[tx_sel%in%names(GTF_annotation$exons_txs)]
  
  tx_cds<-names(GTF_annotation$cds_txs)
  #NOCDS, USE ALL
  b<-disjoin(unlist(reduce(GTF_annotation$exons_txs[GTF_annotation$trann$transcript_id[GTF_annotation$trann$gene_id%in%gens_sel]])),with.revmap=T)
  b$ps<-assay(summarizeOverlaps(b,reads = for_ORFquant$P_sites_all,ignore.strand=F,mode="Union",inter.feature=FALSE))
  
  #b<-b[b%over%GTF_annotation$exons_txs[tx_sel]]
  
  b<-b[as.vector(b$ps>0)]
  
  ov<-findOverlaps(b,GTF_annotation$exons_txs[tx_sel])
  ov<-split(subjectHits(ov),queryHits(ov))
  b$tx_sel<-0
  b$tx_sel[as.numeric(names(ov))]<-elementNROWS(ov)
  
  ov<-findOverlaps(b,GTF_annotation$exons_txs[tx_all])
  ov<-split(subjectHits(ov),queryHits(ov))
  b$tx_all<-0
  b$tx_all[as.numeric(names(ov))]<-elementNROWS(ov)
  
  ov<-findOverlaps(b,GTF_annotation$exons_txs[tx_cds])
  ov<-split(subjectHits(ov),queryHits(ov))
  b$tx_cds<-0
  b$tx_cds[as.numeric(names(ov))]<-elementNROWS(ov)
  
  nts_all<-c()
  nts_cds<-c()
  nts_sel<-c()
  
  
  for(i in c(0,1,2,3,4,5,"more")){
    if(i!="more"){
      nts_all<-c(nts_all,length(b[(b$tx_all)==i]))
      nts_cds<-c(nts_cds,length(b[(b$tx_cds)==i]))
      nts_sel<-c(nts_sel,length(b[(b$tx_sel)==i]))
    }
    if(i=="more"){
      nts_all<-c(nts_all,length(b[(b$tx_all)>5]))
      nts_cds<-c(nts_cds,length(b[(b$tx_cds)>5]))
      nts_sel<-c(nts_sel,length(b[(b$tx_sel)>5]))
    }
  }
  
  juns<-for_ORFquant$junctions
  juns<-juns[juns$reads>0]
  juns$tx_cds<-match(juns$tx_name,tx_cds)
  juns$tx_sel<-match(juns$tx_name,tx_sel)
  juns$tx_cds<-CharacterList(lapply(juns$tx_cds,function(x){x[!is.na(x)]}))
  juns$tx_sel<-CharacterList(lapply(juns$tx_sel,function(x){x[!is.na(x)]}))
  #juns_ok<-juns[elementNROWS(juns$tx_sel)>0 & elementNROWS(juns$tx_cds)>0]
  juns_ok<-juns[sum(juns$gene_id%in%gens_sel)>0]
  
  juns_all<-c()
  juns_cds<-c()
  juns_sel<-c()
  for(i in c(0,1,2,3,4,5,"more")){
    if(i!="more"){
      juns_all<-c(juns_all,sum((elementNROWS(juns_ok$tx_name)==i)))
      juns_cds<-c(juns_cds,sum((elementNROWS(juns_ok$tx_cds)==i)))
      juns_sel<-c(juns_sel,sum((elementNROWS(juns_ok$tx_sel)==i)))
    }
    if(i=="more"){
      juns_all<-c(juns_all,sum((elementNROWS(juns_ok$tx_name)>5)))
      juns_cds<-c(juns_cds,sum((elementNROWS(juns_ok$tx_cds)>5)))
      juns_sel<-c(juns_sel,sum((elementNROWS(juns_ok$tx_sel)>5)))     
    }
  }
  
  juns_all<-round(juns_all/sum(juns_all)*100,digits = 2)
  juns_cds<-round(juns_cds/sum(juns_cds)*100,digits = 2)
  juns_sel<-round(juns_sel/sum(juns_sel)*100,digits = 2)
  
  nts_all<-round(nts_all/sum(nts_all)*100,digits = 2)
  nts_cds<-round(nts_cds/sum(nts_cds)*100,digits = 2)
  nts_sel<-round(nts_sel/sum(nts_sel)*100,digits = 2)
  
  
  
  juns_nt_mapping<-data.frame(t(rbind(juns_all,juns_cds,juns_sel)))
  juns_nt_mapping<-suppressMessages(melt(juns_nt_mapping))
  juns_nt_mapping$type="covered junctions"
  juns_nt_mapping$n_tx<-rep(c(0:5,">5"),3)
  juns_nt_mapping$ref<-c(rep("All Txs",7),rep("Coding Txs",7),rep("Selected Txs",7))
  cds_nt_mapping<-data.frame(t(rbind(nts_all,nts_cds,nts_sel)))
  cds_nt_mapping<-suppressMessages(melt(cds_nt_mapping))
  cds_nt_mapping$type="covered exon bins"
  cds_nt_mapping$n_tx<-rep(c(0:5,">5"),3)
  cds_nt_mapping$ref<-c(rep("All Txs",7),rep("Coding Txs",7),rep("Selected Txs",7))
  
  
  
  df<-rbind(cds_nt_mapping,juns_nt_mapping)
  df$n_tx<-factor(df$n_tx,levels=c("0",">5","5","4","3","2","1"))
  colss<-rev(c(colorRampPalette(colors = c("black","white"))(6),"red"))
  a<-ggplot(df,aes(x=ref,y=value,fill=n_tx))
  a<-a + geom_bar(stat="identity",position = "stack",colour="black")
  #a<-a + scale_y_log10(breaks=c(1,10,100,1000,10000),limits=c(1,15000))
  a<-a + theme_bw()
  a<-a + ylab("percentage")
  a<-a + xlab("")
  a<-a + facet_grid(type~.,scales="free",)
  a<-a + scale_fill_manual(values = colss,"mapping Txs")
  a<-a + theme(axis.title.y = element_text(size=16),axis.text.y  = element_text(angle=45, vjust=0.5, size=16),strip.text.y =  element_text(size=12))
  selection_bins_juns<-a + theme(axis.title.x = element_text(size=12),axis.text.x  = element_text(angle=45, vjust=0.5, size=16))
  
  list_ORFquant_plots[["Sel_txs_bins_juns"]]<-selection_bins_juns
  list_ORFquant_plots[["Sel_txs_bins_juns"]][["pars"]]<-c(5,4.5)
  
  
  
  if(!is.na(coverage_file_plus) & !is.na(coverage_file_minus)){
    
    message("Plotting alternative splice sites profiles ... ",date())
    
    aggregate_regions<-function(range,coverage_plus,coverage_minus,norm_x=NA,norm_y=NA,nozero=F){
      pl<-range[strand(range)=="+"]
      min<-range[strand(range)=="-"]
      mat_pl<-t(sapply(coverage_plus[pl],as.vector))
      if(nozero){
        mat_pl<-mat_pl[rowSums(mat_pl)>0,]
      }
      if(!is.na(norm_x)){
        mat_pl<-mat_pl/(rowSums(mat_pl)/norm_x)
      }
      if(!is.na(norm_y)){
        mat_pl<-mat_pl/(colSums(mat_pl)/norm_y)
      }
      mat_min<-t(sapply(coverage_minus[min],as.vector))
      if(nozero){
        mat_min<-mat_min[rowSums(mat_min)>0,]
      }
      if(!is.na(norm_x)){
        mat_min<-mat_min/(rowSums(mat_min)/norm_x)
      }
      if(!is.na(norm_y)){
        mat_min<-mat_min/(colSums(mat_min)/norm_y)
      }
      
      mat<-rbind(mat_pl,mat_min[,dim(mat_min)[2]:1])
      rownames(mat)<-c(names(pl),names(min))
      mat<-mat[match(names(range),rownames(mat)),]
      mat
    }
    
    ORFs_gen<-ORFquant_results$ORFs_gen
    ORFs_spl_feat_maxORF<-ORFquant_results$ORFs_spl_feat_maxORF
    
    seqlevels(ORFs_gen)<-seqlevels(for_ORFquant$P_sites_all)
    seqlengths(ORFs_gen)<-seqlengths(for_ORFquant$P_sites_all)
    #nope, use coverage
    
    isovals<-cbind.data.frame(ORFs_tx$ORF_id_tr,ORFs_tx$ORF_pct_P_sites_pN)
    ORFs_gen$ORF_pct_P_sites_pN<-isovals[match(names(ORFs_gen),rownames(isovals)),2]
    ORFs_gen$ORF_pct_P_sites_pN[which(is.na(ORFs_gen$ORF_pct_P_sites_pN))]<-0
    ORFs_gen$ORF_pct_P_sites_pN<-round(ORFs_gen$ORF_pct_P_sites_pN,digits = 4)
    covisopl<-coverage(ORFs_gen[strand(ORFs_gen)=="+"],weight = ORFs_gen[strand(ORFs_gen)=="+"]$ORF_pct_P_sites_pN)
    covisomn<-coverage(ORFs_gen[strand(ORFs_gen)=="-"],weight = ORFs_gen[strand(ORFs_gen)=="-"]$ORF_pct_P_sites_pN)
    
    
    covpspl<-import(coverage_file_plus)
    strand(covpspl)<-"+"
    covpspl<-coverage(x = covpspl,weight = as.numeric(covpspl$score))
    covpsmn<-import(coverage_file_minus)
    strand(covpsmn)<-"-"
    covpsmn<-coverage(x = covpsmn,weight = as.numeric(covpsmn$score))
    
    
    #fix the previous annotation, then re-run
    listvals<-list()
    listiso<-list()
    listrang<-list()
    #CHECK DISTS
    types<-c("up_5ss","down_5ss","same_5ss","up_3ss","down_3ss","same_3ss","CDS_spanning")
    types_alt<-types[grep(types,pattern = "same",invert = T)]
    orfs_alt<-ORFs_spl_feat_maxORF[sapply(strsplit(ORFs_spl_feat_maxORF$spl_type,";"),function(x){sum(x%in%types_alt)>0})]
    ORFs_ok<-names(orfs_alt)
    genes_ok<-unique(ORFs_tx$gene_id[ORFs_tx$ORF_id_tr%in%ORFs_ok])
    ORFs_ok<-ORFs_tx$ORF_id_tr[ORFs_tx$gene_id%in%genes_ok]
    ORFs_ok<-ORFs_spl_feat_maxORF[names(ORFs_spl_feat_maxORF)%in%ORFs_ok]
    
    
    for(i in types){
      region<-ORFs_ok[grep(i,ORFs_ok$spl_type)]
      region<-region[!duplicated(GRanges(region))]
      if(i%in%c("up_3ss","up_lastCDS")){
        region<-region[elementNROWS(region$ref)==1]
        set.seed(666)
        if(length(region)>1000){region<-region[sample(x = 1:length(region),size = 1000,replace = F)]}
        region$dist<-end(region$ref)-end(region)
        region$dist[as.vector(strand(region)=="-")]<-start(region[as.vector(strand(region)=="-")])-start(region[as.vector(strand(region)=="-")]$ref)
        region$dist<-as.numeric(region$dist)
        dst<-region$dist
        nms<-names(region)
        region1<-promoters(resize(region,width = 1,fix = "end"),upstream = 25,downstream = 0)
        mcols(region1)<-mcols(region)
        
        region2<-promoters(resize(region,width = 1,fix = "end"),upstream = 0,downstream = 25)
        mcols(region2)<-mcols(region)
        
        region<-promoters(resize(region,width = 1,fix = "end"),upstream = 25,downstream = 25)
        mcols(region)<-mcols(region1)
        
        region$dist<-dst
        names(region)<-nms
      }
      
      if(i%in%c("down_3ss","down_lastCDS","same_3ss","same_lastCDS")){
        region<-region[elementNROWS(region$ref)==1]
        set.seed(666)
        if(length(region)>1000){region<-region[sample(x = 1:length(region),size = 1000,replace = F)]}
        region$dist<-end(region)-end(region$ref)
        region$dist[as.vector(strand(region)=="-")]<-start(region[as.vector(strand(region)=="-")]$ref)-start(region[as.vector(strand(region)=="-")])
        region$dist<-as.numeric(region$dist)
        dst<-region$dist
        nms<-names(region)
        region1<-promoters(resize(unlist(region$ref),width = 1,fix = "end"),upstream = 25,downstream = 0)
        mcols(region1)<-mcols(region)
        
        region2<-promoters(resize(unlist(region$ref),width = 1,fix = "end"),upstream = 0,downstream = 25)
        mcols(region2)<-mcols(region)
        
        region<-promoters(resize(unlist(region$ref),width = 1,fix = "end"),upstream = 25,downstream = 25)
        mcols(region)<-mcols(region2)
        
        region$dist<-dst
        names(region)<-nms
      }
      
      
      if(i%in%c("down_5ss","down_firstCDS")){
        region<-region[elementNROWS(region$ref)==1]
        set.seed(666)
        if(length(region)>1000){region<-region[sample(x = 1:length(region),size = 1000,replace = F)]}
        region$dist<-start(region)-start(region$ref)
        region$dist[as.vector(strand(region)=="-")]<-end(region[as.vector(strand(region)=="-")]$ref)-end(region[as.vector(strand(region)=="-")])
        region$dist<-as.numeric(region$dist)
        #invert numbers for donwstr-upstr
        dst<-region$dist
        nms<-names(region)
        region2<-promoters(resize(region,width = 1,fix = "start"),upstream = 25,downstream = 0)
        mcols(region2)<-mcols(region)
        
        region1<-promoters(resize(region,width = 1,fix = "start"),upstream = 0,downstream = 25)
        mcols(region1)<-mcols(region)
        
        region<-promoters(resize(region,width = 1,fix = "start"),upstream = 25,downstream = 25)
        mcols(region)<-mcols(region1)
        
        region$dist<-dst
        names(region)<-nms
      }
      
      if(i%in%c("up_5ss","up_firstCDS","same_5ss","same_firstCDS")){
        region<-region[elementNROWS(region$ref)==1]
        set.seed(666)
        if(length(region)>1000){region<-region[sample(x = 1:length(region),size = 1000,replace = F)]}
        region$dist<-start(region$ref)-start(region)
        region$dist[as.vector(strand(region)=="-")]<-end(region[as.vector(strand(region)=="-")])-end(region[as.vector(strand(region)=="-")]$ref)
        region$dist<-as.numeric(region$dist)
        dst<-region$dist
        nms<-names(region)
        region2<-promoters(resize(unlist(region$ref),width = 1,fix = "start"),upstream = 25,downstream = 0)
        mcols(region2)<-mcols(region)
        region1<-promoters(resize(unlist(region$ref),width = 1,fix = "start"),upstream = 0,downstream = 25)
        mcols(region1)<-mcols(region)
        region<-promoters(resize(unlist(region$ref),width = 1,fix = "start"),upstream = 25,downstream = 25)
        mcols(region)<-mcols(region1)
        region$dist<-dst
        names(region)<-nms
      }
      
      if(!i%in%c("CDS_spanning")){
        region1$iso<-0
        region1$iso[as.vector(strand(region1)=="+")]<-mean(covisopl[region1[as.vector(strand(region1)=="+")]])
        region1$iso[as.vector(strand(region1)=="-")]<-mean(covisomn[region1[as.vector(strand(region1)=="-")]])
        region2$iso<-0
        region2$iso[as.vector(strand(region2)=="+")]<-mean(covisopl[region2[as.vector(strand(region2)=="+")]])
        region2$iso[as.vector(strand(region2)=="-")]<-mean(covisomn[region2[as.vector(strand(region2)=="-")]])
        region$iso_1<-region1$iso
        region$iso_2<-region2$iso
        region$delta_iso<-region$iso_1-region$iso_2
        #check for those depending on region type
        #region<-region[region$iso_1<101 & region$iso_2<101]
        #region<-region[region$delta_iso<1]
        listrang[[i]]<-region
        listvals[[i]]<-aggregate_regions(range=region,coverage_plus=covpspl,coverage_minus=covpsmn)
        listiso[[i]]<-aggregate_regions(range=region,coverage_plus=covisopl,coverage_minus=covisomn)
        
        
      }
      
      if(i%in%c("CDS_spanning")){
        region<-region[elementNROWS(region$ref)==2]
        ress<-unlist(region$ref)
        start(ress[seq(1,length(region)*2,by = 2)])<-end(ress[seq(1,length(region)*2,by = 2)])
        
        region1_a<-promoters(ress[seq(1,length(region)*2,by = 2)],upstream = 25,downstream = 0)
        region2_a<-promoters(ress[seq(1,length(region)*2,by = 2)],upstream = 0,downstream = 25)
        region2_b<-promoters(ress[seq(2,length(region)*2,by = 2)],upstream = 25,downstream = 0)
        region1_b<-promoters(ress[seq(2,length(region)*2,by = 2)],upstream = 0,downstream = 25)
        region_a<-promoters(ress[seq(1,length(region)*2,by = 2)],upstream = 25,downstream = 25)
        names(region_a)<-names(region)
        region_b<-promoters(ress[seq(2,length(region)*2,by = 2)],upstream = 25,downstream = 25)
        names(region_b)<-names(region)
        region_a$dist<-sapply(region$ref,function(x){min(width(gaps(x)))})
        region_b$dist<-sapply(region$ref,function(x){min(width(gaps(x)))})
        
        #check for upstr and dowstr, or maybe main and secondary
        
        region1_a$iso<-0
        region1_a$iso[as.vector(strand(region1_a)=="+")]<-mean(covisopl[region1_a[as.vector(strand(region1_a)=="+")]])
        region1_a$iso[as.vector(strand(region1_a)=="-")]<-mean(covisomn[region1_a[as.vector(strand(region1_a)=="-")]])
        region1_b$iso<-0
        region1_b$iso[as.vector(strand(region1_b)=="+")]<-mean(covisopl[region1_b[as.vector(strand(region1_b)=="+")]])
        region1_b$iso[as.vector(strand(region1_b)=="-")]<-mean(covisomn[region1_b[as.vector(strand(region1_b)=="-")]])
        region2_a$iso<-0
        region2_a$iso[as.vector(strand(region2_a)=="+")]<-mean(covisopl[region2_a[as.vector(strand(region2_a)=="+")]])
        region2_a$iso[as.vector(strand(region2_a)=="-")]<-mean(covisomn[region2_a[as.vector(strand(region2_a)=="-")]])
        region2_b$iso<-0
        region2_b$iso[as.vector(strand(region2_b)=="+")]<-mean(covisopl[region2_b[as.vector(strand(region2_b)=="+")]])
        region2_b$iso[as.vector(strand(region2_b)=="-")]<-mean(covisomn[region2_b[as.vector(strand(region2_b)=="-")]])
        
        
        
        region_a$iso_1<-region1_a$iso
        region_a$iso_2<-region2_a$iso
        region_b$iso_1<-region1_b$iso
        region_b$iso_2<-region2_b$iso
        
        region_a$delta_iso<-region_a$iso_1-region_a$iso_2
        region_b$delta_iso<-region_b$iso_1-region_b$iso_2
        names(region_a)<-paste(names(region),"left",sep = "_")
        names(region_b)<-paste(names(region),"right",sep = "_")
        region<-c(region_a,region_b)
        
        listrang[[i]]<-region
        listvals[[i]]<-aggregate_regions(range=region,coverage_plus=covpspl,coverage_minus=covpsmn)
        listiso[[i]]<-aggregate_regions(range=region,coverage_plus=covisopl,coverage_minus=covisomn)
        
      }
    }
    
    same_all<-ORFs_ok[grep("same_5ss;same_3ss",ORFs_ok$spl_type)]
    same_all<-same_all[width(same_all)>110]
    same_all$dist<-50
    same_all$iso_1<-50
    same_all$iso_2<-50
    same_all$delta_iso<-50
    set.seed(666)
    if(length(same_all)>1000){same_all<-same_all[sample(x = 1:length(same_all),size = 1000,replace = F)]}
    same_all_5<-promoters(resize(same_all,width = 1,fix = "center"),upstream = 50,downstream = 0)
    same_all_3<-promoters(resize(same_all,width = 1,fix = "center"),upstream = 0,downstream = 50)
    
    listrang[["sameall_5"]]<-same_all_5
    listvals[["sameall_5"]]<-aggregate_regions(range=same_all_5,coverage_plus=covpspl,coverage_minus=covpsmn)
    listiso[["sameall_5"]]<-aggregate_regions(range=same_all_5,coverage_plus=covisopl,coverage_minus=covisomn)
    
    listrang[["sameall_3"]]<-same_all_3
    listvals[["sameall_3"]]<-aggregate_regions(range=same_all_3,coverage_plus=covpspl,coverage_minus=covpsmn)
    listiso[["sameall_3"]]<-aggregate_regions(range=same_all_3,coverage_plus=covisopl,coverage_minus=covisomn)
    
    
    topl_ok<-c("up_5ss","down_3ss","down_5ss","up_3ss","CDS_spanning")
    topl_cons<-c("same_5ss","same_3ss","sameall_5","sameall_3","CDS_spanning")
    
    plots_ribo<-list()
    plots_iso<-list()
    plots_sketch<-list()
    
    for(i in 1:length(topl_ok)){
      rans<-listrang[[topl_ok[i]]]
      vals<-listvals[[topl_ok[i]]]
      vals_contr<-listvals[[topl_cons[i]]]
      isos_contr<-listiso[[topl_cons[i]]]
      rans_contr<-listrang[[topl_cons[i]]]
      if(topl_ok[i]=="CDS_spanning"){
        vals_contr<-rbind(listvals[["same_3ss"]],listvals[["same_5ss"]])
        isos_contr<-rbind(listiso[["same_3ss"]],listiso[["same_5ss"]])
        rans_contr<-c(listrang[["same_3ss"]],listrang[["same_5ss"]])
        names(rans_contr)<-paste(names(rans_contr),c(rep("left",length(listrang[["same_3ss"]])),rep("right",length(listrang[["same_3ss"]]))),sep = "_")
        
      }
      
      isos<-listiso[[topl_ok[i]]]
      rans$id<-paste(names(rans),GRanges(rans),sep = "_")
      #ok<-intersect(intersect(which(abs(rans$dist)>25),which(rans$iso_1<100 & rans$iso_2<100)),which(rans$delta_iso<0))
      ok<-intersect(which(abs(rans$dist)>25),which(rans$iso_1<101 & rans$iso_2<101))
      ok<-intersect(ok,which(rans$delta_iso>0))
      ok<-intersect(ok,which(rowSums(vals)>.1))
      
      if(topl_ok[i]=="down_5ss"){
        blocks=data.frame(x1=c(.5,0,0,.14,.28,.42), x2=c(1,1,.09,.23,.37,.5), y1=c(.5,.05,.675,.675,.675,.675), y2=c(.85,.4,.685,.685,.685,.685),coll=c("red","blue","red","red","red","red"),stringsAsFactors = F)
      }
      
      if(topl_ok[i]=="down_firstCDS"){
        blocks=data.frame(x1=c(.5,0,0), x2=c(1,1,1), y1=c(.5,.05,.575), y2=c(.85,.4,.785),coll=c("red","blue","red"),stringsAsFactors = F)
      }
      
      if(topl_ok[i]=="up_5ss"){
        blocks=data.frame(x1=c(0,.5,0,.14,.28,.42), x2=c(1,1,.09,.23,.37,.5), y1=c(.5,.05,.225,.225,.225,.225), y2=c(.85,.4,.235,.235,.235,.235),coll=c("red","blue","blue","blue","blue","blue"),stringsAsFactors = F)
      }
      
      if(topl_ok[i]=="up_firstCDS"){
        blocks=data.frame(x1=c(0,.5,0), x2=c(1,1,1), y1=c(.5,.05,.125), y2=c(.85,.4,.330),coll=c("red","blue","blue"),stringsAsFactors = F)
      }
      
      if(topl_ok[i]=="down_3ss"){
        blocks=data.frame(x1=c(0,0,.5,.64,.78,.92), x2=c(1,.5,.59,.73,.87,1), y1=c(.5,.05,.225,.225,.225,.225), y2=c(.85,.4,.235,.235,.235,.235),coll=c("red","blue","blue","blue","blue","blue"),stringsAsFactors = F)
      }
      
      if(topl_ok[i]=="down_lastCDS"){
        blocks=data.frame(x1=c(0,0,.5), x2=c(1,.5,1), y1=c(.5,.05,.125), y2=c(.85,.4,.330),coll=c("red","blue","blue"),stringsAsFactors = F)
      }
      
      
      if(topl_ok[i]=="up_3ss"){
        blocks=data.frame(x1=c(0,0,.5,.64,.78,.92), x2=c(.5,1,.59,.73,.87,1), y1=c(.5,.05,.675,.675,.675,.675), y2=c(.85,.4,.685,.685,.685,.685),coll=c("red","blue","red","red","red","red"),stringsAsFactors = F)
      }
      
      if(topl_ok[i]=="up_lastCDS"){
        blocks=data.frame(x1=c(0,0,.5), x2=c(.5,1,1), y1=c(.5,.05,.575), y2=c(.85,.4,.785),coll=c("red","blue","red"),stringsAsFactors = F)
      }
      
      
      if(topl_ok[i]=="CDS_spanning"){
        blocks=data.frame(x1=c(0,0,0.25,0.39,0.53,0.67,.75), x2=c(1,.25,0.34,0.48,0.62,0.75,1), y1=c(.5,.05,.225,.225,.225,.225,.05), y2=c(.85,.4,.235,.235,.235,.235,.4),coll=c("red","blue","blue","blue","blue","blue","blue"),stringsAsFactors = F)
      }
      
      
      colss<-blocks$coll
      blocks$coll<-factor(blocks$coll,levels=c("red","blue","white"))
      df<-blocks
      bl<-ggplot() + geom_rect(data=df, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2,fill=coll,colour=coll))+
        scale_fill_manual(values = colss) +
        scale_color_manual(values = colss) + 
        theme_nothing() + ylim(c(0,1))
      
      
      if(length(ok)==0){
        a<-ggplot() + theme_void() + ggtitle(paste("Event =",topl_ok[i],"\nNot enough distinct regions"))
        b<-a
        plots_ribo[[i]]<-a
        plots_iso[[i]]<-b
        plots_sketch[[i]]<-bl
        next
      }
      
      rans_contr<-rans_contr[rowSums(vals_contr)>.1]
      isos_contr<-isos_contr[rowSums(vals_contr)>.1,]
      vals_contr<-vals_contr[rowSums(vals_contr)>.1,]
      vals<-vals[ok,]
      isos<-isos[ok,]
      rans<-rans[ok]
      rownames(vals)<-rans$id
      rownames(isos)<-rans$id
      rans$delta_iso<-round(as.numeric(abs(rans$delta_iso)),digits=2)
      rans$delta_iso2<-abs(rans$iso_1-rans$delta_iso)
      if(topl_ok[i]%in%c("up_3ss","down_5ss","down_firstCDS","up_lastCDS")){
        rans$delta_iso2<-abs(rans$iso_1-rans$iso_2)
      }
      rans$delta_iso<-rans$delta_iso2
      qnt<-unique(quantile(rans$delta_iso,probs=seq(0,1,length.out = 4),include.lowest = T))
      
      rans$group<-as.character(cut(rans$delta_iso,breaks = qnt,include.lowest = T))
      
      tb<-aggregate(rans$delta_iso,by=list(rans$group),mean)
      rans$group<-round(tb$x[match(rans$group,tb$Group.1)],digits = 2)
      ok<-which(!is.na(rans$group))
      vals<-vals[ok,]
      isos<-isos[ok,]
      rans<-rans[ok]
      rans_contr$id<-paste(names(rans_contr),GRanges(rans_contr),sep = "_")
      rans_contr$group<-"No_mixture"
      rownames(isos_contr)<-rans_contr$id
      rownames(vals_contr)<-rans_contr$id
      levss<-sort(unique(rans$group),decreasing = F)
      levss<-c("No_mixture",levss)
      if(topl_ok[i]=="CDS_spanning"){
        rans$spanning<-"downstream"
        rans$spanning[grep(names(rans),pattern = "left")]<-"upstream"
        rans_contr$spanning<-"downstream"
        rans_contr$spanning[grep(names(rans_contr),pattern = "left")]<-"upstream"
      }
      notcol<-names(mcols(rans))[!names(mcols(rans))%in%names(mcols(rans_contr))]
      rans$ref<-NULL
      mcols(rans_contr)[,notcol]<-""
      mcols(rans_contr)<-mcols(rans_contr)[,names(mcols(rans))]
      rans<-c(rans,rans_contr)
      rans$group<-factor(rans$group,levels=levss)
      
      
      #HERE IT WAS THE PROBLEM
      #levss<-levels(rans$group)[order(as.numeric(sapply(strsplit(unique(rans$group),"-"),"[[",1)),decreasing = F)]
      #levels(rans$group)<-levss
      
      
      #vals<-t(apply(vals,1,scale))
      #vals<-vals/apply(vals,1,max)
      
      vals<-rbind(vals,vals_contr)
      isos<-rbind(isos,isos_contr)
      
      vals<-vals/rowSums(vals)
      #vals<-scale(vals)
      colnames(vals)<--25:24
      vls<-DataFrame(melt(vals))
      vls$delta_iso<-rans$delta_iso[match(vls$Var1,rans$id)]
      vls$group<-rans$group[match(vls$Var1,rans$id)]
      #isos<-t(apply(isos,1,scale))
      #isos<-isos/apply(isos,1,max)
      
      isos<-isos/rowSums(isos)
      #isos<-scale(isos)
      colnames(isos)<--25:24
      
      vls2<-DataFrame(melt(isos))
      vls2$delta_iso<-rans$delta_iso[match(vls2$Var1,rans$id)]
      
      vls2$group<-rans$group[match(vls2$Var1,rans$id)]
      vls$type<-"Ribo-seq_coverage"
      vls2$type<-"Iso_values"
      if(topl_ok[i]=="CDS_spanning"){
        vls2$spanning<-rans$spanning[match(vls2$Var1,rans$id)]
        vls$spanning<-rans$spanning[match(vls$Var1,rans$id)]
      }
      
      vlsall<-rbind(vls,vls2)
      vlsall[,"Position_from_splice_site"]<-vlsall$Var2
      vlsall$type<-factor(vlsall$type)
      
      
      dat<-as.data.frame(vlsall)
      cols<- colorRampPalette(c("dark blue","blue","red"))(8)[c(1,3,6,7)]
      cols<-cols[(4-length(qnt)+1):4]
      
      means<-aggregate(dat$value,by=list(dat$group,dat$Position_from_splice_site,dat$type),mean)
      colnames(means)<-c("group","Position_from_splice_site","type","value")
      
      if(topl_ok[i]=="CDS_spanning"){
        means<-aggregate(dat$value,by=list(dat$group,dat$Position_from_splice_site,dat$type,dat$spanning),mean)
        colnames(means)<-c("group","Position_from_splice_site","type","spanning","value")
        ok<-ok[1:(length(ok)/2)]
        means$spanning<-factor(means$spanning,levels=c("upstream","downstream"))
      }
      
      
      means$groups_all<-paste(means$group,means$type)
      means$shape<-means$type
      
      means$color<-cols[as.numeric(means$group)]
      means_ribo<-means[means$type!="Iso_values",]
      means_Iso<-means[means$type=="Iso_values",]
      #ss<-aes(aes(ymax = value + sd, ymin = value - sd,group = groups_all,color=groups_all,shape=type,linetype=type))
      nm_titl<-topl_ok[i]
      
      a<-ggplot(means_ribo, aes(x = Position_from_splice_site,y =  value, group = group,color=group)) +
        geom_line(size=1.8) +
        geom_vline(xintercept = 0,col="black",lty=2) +
        theme_classic() +
        #geom_ribbon(aes(ymin=value-sds, ymax=value+sds, x = Position_from_splice_site, fill = group), alpha = 0.01) +
        ylab("Normalized coverage \nRibo-seq") +
        theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(angle=45, vjust=0.5, size=15)) +
        theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=15))+
        scale_color_manual(values=cols,name="Avg.\ncontribution\nadditional ORF(s)") +
        ggtitle(paste("Event =",nm_titl,"\nn =",length(ok)))
      
      b<-ggplot(means_Iso, aes(x = Position_from_splice_site,y =  value, group = group,color=group)) +
        geom_line(size=1.8) +
        geom_vline(xintercept = 0,col="black",lty=2) +
        theme_classic() +
        #geom_ribbon(aes(ymin=value-sds, ymax=value+sds, x = Position_from_splice_site, fill = group), alpha = 0.01) +
        theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(angle=45, vjust=0.5, size=15)) +
        theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=15))+
        scale_color_manual(values=cols,name="Avg.\ncontribution\nadditional ORF(s)") +
        ggtitle(paste("Event =",nm_titl,"\nn =",length(ok)))
      if(topl_ok[i]=="CDS_spanning"){
        a<-a+facet_grid(facets = .~spanning) + xlab("Position from Splice Sites") +  ylab("Normalized coverage\nRibo-seq")
        b<-b+facet_grid(facets = .~spanning) + xlab("Position from Splice Sites") +  ylab("Normalized coverage\nORFquant ORFs")
      }
      if(topl_ok[i]%in%c("up_5ss")){
        a<-a + xlab("Position from Splice Site") +  ylab("Normalized coverage\nRibo-seq")
        b<-b + xlab("Position from Splice Site") +  ylab("Normalized coverage\nORFquant ORFs")
      }
      if(topl_ok[i]%in%c("up_firstCDS")){
        a<-a + xlab("Position from Start codon") +  ylab("Normalized coverage\nRibo-seq")
        b<-b + xlab("Position from Start codon") +  ylab("Normalized coverage\nORFquant ORFs")
      }
      
      if(topl_ok[i]%in%c("up_lastCDS")){
        a<-a + xlab("Position from Stop codon") +  ylab("")
        b<-b + xlab("Position from Stop codon") +  ylab("")
      }
      
      
      if(!topl_ok[i]%in%c("up_lastCDS","up_firstCDS","up_5ss","CDS_spanning")){
        a<-a + xlab("") +  ylab("")
        b<-b + xlab("") +  ylab("")
      }
      
      
      
      plots_ribo[[i]]<-a
      plots_iso[[i]]<-b
      plots_sketch[[i]]<-bl
      
    }
    
    x1s<-c(0,.25,.5,.75,.3)
    x1s_iso<-x1s
    x1s_sk<-x1s
    
    wids<-c(rep(.25,4),.5)
    wids_iso<-wids
    wids_sk<-wids
    
    y2s<-c(rep(.8,4),.3)
    y2s_iso<-y2s-.2
    y2s_sk<-round(as.numeric(y2s_iso-.1),digits = 3)
    
    heis<-c(rep(.2,4),.2)
    heis_iso<-heis
    heis_sk<-rep(.1,10)
    
    all<-ggdraw()
    for(i in 1:length(plots_ribo)){
      all<-all+draw_plot(plots_ribo[[i]], x = x1s[i], y = y2s[i], width = wids[i], height = heis[i])  +
        draw_plot(plots_iso[[i]],  x = x1s_iso[i], y = y2s_iso[i], width = wids_iso[i], height = heis_iso[i])  +
        draw_plot(plots_sketch[[i]],  x = x1s_sk[i], y = y2s_sk[i], width = wids_sk[i], height = heis_sk[i])
      
    }
    
    all_metapl<-all
    list_ORFquant_plots[["Meta_splicing_coverage"]]<-all_metapl
    list_ORFquant_plots[["Meta_splicing_coverage"]][["pars"]]<-c(22,16)
  }
  
  save(list_ORFquant_plots,file = paste0(output_plots_path,"/",prefix,"_","ORFquant_plots_RData"))
  for(i in names(list_ORFquant_plots)){
    pdf(file = paste0(output_plots_path,"/",prefix,"_",i,".pdf"),width = list_ORFquant_plots[[i]][["pars"]][1] ,height = list_ORFquant_plots[[i]][["pars"]][2])
    print(list_ORFquant_plots[[i]])
    dev.off()
  }
  message("Plotting ORFquant results for ",ORFquant_output_file,"  --- Done! ",date())
  
}


#' Create a plot of the ORFquant results at a locus
#'
#' Create a plot of the ORFquant results at a locus. Uses the info form the orfquant results about where the psite data is located.
#' 
#' @keywords ORFquant
#' @author Dermot Harnett, \email{dermot.p.harnett@@gmail.com}
#' 
#' @param locus String; a gene name, must be present in names(orfquant_results$ORFs_gen)
#' @param orfquant_results A list containing processed output from ORFquant
#' @param bam_files bam files, (or pre-processed bam data from RiboseQC) to be plotted
#' @param plotfile the file into which the plot will be saved as a pdf
#' @param col not used (the P-site and junction tracks are drawn in forestgreen)
#' @return returns the value of plotfile if successfull.
#' @export
#' 
plot_orfquant_locus<-function(locus,orfquant_results,bam_files, plotfile='locusplot.pdf', col ='green' ){
  if (!all(vapply(c("Gviz",'lemon','dplyr','GenomeInfoDb'), requireNamespace, logical(1), quietly = TRUE))) {
    stop("Packages \"Gviz\",\"dplyr\",\"lemon\",\"GenomeInfoDb\" needed for this function to work. Please install them.",
         call. = FALSE)
  }
  options(ucscChromosomeNames=FALSE)
  `%>%` <- dplyr::`%>%`
  
  if(is.null(orfquant_results$psite_data_file)){stop("this object looks like it's form an old version of ORFquant, it doesn't list the psite data file")}
  if(length(orfquant_results$psite_data_file)>1){stop("locus plots aren't supported for multiple psite tracks - either unify the psite tracks or modify the input object to have only one track")}
  riboseqcoutput<-get(load(orfquant_results$psite_data_file))
  
  anno <- GTF_annotation
  
  
  selgene <- locus
  seltxs <- orfquant_results$ORFs_tx[orfquant_results$ORFs_tx$gene_id==selgene]$transcript_id
  selorfs <- orfquant_results$ORFs_tx[orfquant_results$ORFs_tx$gene_id==selgene]$ORF_id_tr
  stopifnot(length(locus)==1)
  stopifnot(any(orfquant_results$ORFs_tx$gene_id==selgene))
  orfs_quantified_gen <-  orfquant_results$ORFs_gen[selorfs]
  
  mcols(orfs_quantified_gen) <- mcols(orfquant_results$ORFs_tx)[match(names(orfs_quantified_gen),orfquant_results$ORFs_tx$ORF_id_tr),]
  seltxs <- orfs_quantified_gen$transcript_id%>%unique
  orfs_quantified_gen$feature <- 'CDS'
  orfs_quantified_gen$transcript=orfs_quantified_gen$transcript_id
  anno <- GTF_annotation
  orfs_quantified_tr <- anno$exons_txs
  orfs_quantified_tr <- orfs_quantified_tr[unique(orfs_quantified_gen$transcript_id)]
  #get the orfs, then get the negative coverage
  #add transcript info to the ORFs_tx object
  seqinf <- Seqinfo(names(anno$exons_txs),anno$exons_txs%>%width%>%sum)
  #Now get the negatives for each ORF
  
  utrs <- orfquant_results$ORFs_tx%>%subset(gene_id==selgene)%>%GenomeInfoDb::keepSeqlevels(seltxs)%>%{seqinfo(.)<-seqinf[seltxs];.}%>%
    coverage%>%
    as('GRanges')%>%subset(score==0)%>%mapFromTranscripts(anno$exons_txs)%>%
    {.$transcript <- names(anno$exons_txs)[.$transcriptsHits];.}%>%
    {.$feature='utr';.}
  
  orfquantgr <- c(
    orfs_quantified_gen[,c('feature','transcript')]
    ,utrs[,c('feature','transcript')]
  )
  orfquantgr$feature[orfquantgr$feature=='CDS'] <- names(orfquantgr)[orfquantgr$feature=='CDS']
  #get correct col name
  metacols = colnames(mcols(orfs_quantified_gen))
  quantcol = 'ORFs_pM'
  
  ufeats <- orfquantgr$feature%>%{.=.[.!='utr'];.}%>%unique
  orfscores<- ufeats%>%setNames(match(ufeats,orfs_quantified_gen$ORF_id_tr)%>%{mcols(orfs_quantified_gen[.])[[quantcol]]},.)
  
  
  orfcols <- orfscores%>%{./max(na.omit(.))}%>%
    c(0,.)
  
  orfquantgrscores = mcols(orfs_quantified_gen)[[quantcol]][match(names(orfquantgr),orfs_quantified_gen$ORF_id_tr)]
  
  orfcols <- orfcols%>% vapply(function(.)tryCatch({rgb(0,.,0)},error=function(e){'white'}),'foo')%>%setNames(c('0',ufeats))
  
  
  ###Define non selected
  disctxs<-anno$txs_gene[selgene]%>%unlist%>%.$tx_name%>%unique%>%setdiff(seltxs)
  disctxsint <- disctxs%>%intersect(seqnames(anno$cds_txs_coords))%>%as.character 
  disc_orfquantgr <- if(length(disctxsint)==0) GRanges(transcript=character(0),feature=character(0)) else anno$cds_txs_coords%>%
    GenomeInfoDb::keepSeqlevels(disctxsint,'coarse')%>%
    {seqinfo(.)<-seqinf[disctxsint];.}%>%  coverage%>%
    as('GRanges')%>%
    subset(.$score==0)%>%
    {   
      txgr = .
      out = mapFromTranscripts(txgr,anno$exons_txs)
      out$score = txgr$score[out$xHits]
      out
    }%>%
    {.$transcript <- names(anno$exons_txs)[.$transcriptsHits];.}%>%
    {.$feature=ifelse(.$score==0,'utr','CDS');.}
  disc_orfquantgr <- disc_orfquantgr%>% c(.,anno$cds_txs[disctxsint]%>%unlist%>%{.$feature=rep('CDS',length(.));.$transcript=names(.);.})
  discORFnames<-paste0(disctxsint,'_',start(anno$cds_txs_coords[disctxsint]),'_',end(anno$cds_txs_coords[disctxsint]))%>%setNames(disctxsint)
  disc_orfquantgr$symbol = discORFnames[disc_orfquantgr$transcript]
  fakejreads <- riboseqcoutput$junctions%>%subset(any(gene_id==selgene))%>%resize(width(.)+2,'center')%>%
    {.$cigar <- sprintf('1M%dN1M',width(.)-2);.}
  fakejreads <- fakejreads[mapply(seq_along(fakejreads[]),fakejreads$reads,FUN=rep)%>%unlist]
  ncols <- 2
  nrows <- 1
  orfcols <- orfcols[order(-orfscores[names(orfcols)])]
  orfquantgr_sorted <- orfquantgr[order(orfscores[names(orfquantgr)])]
  fix_utrs <- function(orfquantgr_sorted){
    orfquantgr_sorted$symbol = names(orfquantgr_sorted)
    #add utrs for each ORF
    #for each selected ORF
    orftrpairs<-orfquantgr_sorted%>%subset(feature!='utr')%>%mcols%>%as.data.frame%>%{dplyr::distinct(.)}
    orfutrs <- orfquantgr_sorted%>%subset(feature=='utr')
    orfutrs<-lapply(1:nrow(orftrpairs),function(i){
      orfutrs <- orfutrs%>%subset(transcript==orftrpairs$transcript[i])
      orfutrs$symbol = orftrpairs$feature[i]
      names(orfutrs) = orfutrs$symbol 
      orfutrs
    })%>%GRangesList%>%unlist
    orfquantgr_sorted<-orfquantgr_sorted%>%subset(feature!='utr')%>%c(.,orfutrs)
    orfquantgr_sorted
  }
  orfquantgr_sorted<-fix_utrs(orfquantgr_sorted)
  # disc_orfquantgrfix<-fix_utrs(disc_orfquantgr)
  disc_orfquantgrfix<-(disc_orfquantgr)
  #dimenions, extent of the plot
  #genes on several chromosomes or strands (e.g. UCSC's PAR genes) aren't in anno$genes: then the range of the
  #gene's transcripts on the chromosome and strand of the first selected transcript
  selgenerange <-  if(selgene%in%names(anno$genes)) anno$genes[selgene] else
    range(anno$txs_gene[[selgene]])%>%subsetByOverlaps(anno$exons_txs[[seltxs[1]]])
  plotstart = start(selgenerange) - (0.2 * (end(selgenerange)-start(selgenerange)))
  plotend = end(selgenerange) + (0 * (end(selgenerange)-start(selgenerange)))
  legendwidth=1/10
  plottitle <- paste0('ORFquant: ',selgene)
  #write to pdf
  #Gviz needs 3 points per row of a stacked track: a track with more than 21 transcripts gets more height, and the page with it (1 inch per unit of sizes)
  plotsizes <- pmax(1,c(1,length(disctxs),length(seltxs),1,1,length(disctxsint),length(selorfs))/21)
  pdf(plotfile,width=14+2,height=sum(plotsizes))
  pdfdev <- grDevices::dev.cur()
  on.exit(dev.off(pdfdev),add=TRUE)
  #code for arranging legend next to the locus plot
  grid.newpage()
  vp1 <- viewport(x = 0, y = 0, width = 1-legendwidth*1.5, height = 1,
                  just = c("left", "bottom"), name = "vp1")
  vp2 <- viewport(x = 1-legendwidth*1.5, y = 0, width = legendwidth*1.5, height = 1,
                  just = c("left", "bottom"))
  vp3 <- viewport(x = 1-legendwidth*1.5, y = 0, width = legendwidth*1.5, height = 1/7,
                  just = c("left", "bottom"))
  pushViewport(vp1)
  #finally plot the locus
  Gviz::plotTracks(main=plottitle,cex.main=2,legend=TRUE,add=TRUE,
                   from=plotstart,to=plotend,#zoomed in on the orf in question
                   sizes=plotsizes,rot.title=0,cex.title=1,title.width=2.5,
                   c(
                     Gviz::GenomeAxisTrack(range=selgenerange),
                     # Gviz::rnaseqtrack, # plot the riboseq signal
                     # Gviz::txs_discarded_Track,
                     # Gviz::txs_selected_track,
                     # Gviz::DataTrack(riboseqcoutput$P_sites_all%>%subsetByOverlaps(selgenerange),type='hist'),
                     Gviz::GeneRegionTrack(name='discarded\ntranscripts',anno$exons_txs[disctxs]%>%unlist%>%{.$transcript=names(.);.$feature=rep('exon',length(.));.},fill='#F7CAC9',
                                           transcriptAnnotation='transcript'),
                     Gviz::GeneRegionTrack(exon='forestgreen',name='selected\ntranscripts',anno$exons_txs[seltxs]%>%unlist%>%{.$transcript=names(.);.$feature=rep('exon',length(.));.},fill='#F7CAC9',
                                           transcriptAnnotation='transcript'),
                     Gviz::DataTrack(legend=TRUE,name='\t\t P-Sites',col.histogram='forestgreen',riboseqcoutput$P_sites_all%>%subsetByOverlaps(selgenerange),type='hist'),
                     Gviz::AlignmentsTrack(name='\nJunction Reads\n\n\n',col.sashimi='forestgreen',fakejreads[,'cigar'],type='sashimi',sashimiNumbers=TRUE),
                     # Gviz::GeneRegionTrack(discarded_orfs_gen),
                     Gviz::GeneRegionTrack(name='Discarded\nORFs',disc_orfquantgrfix,
                                           transcriptAnnotation='symbol',collapse=FALSE,thinBoxFeature='utr',CDS='blue',utr='white'),
                     Gviz::GeneRegionTrack(name='Selected\nORFs',
                                           range=orfquantgr_sorted,collapse=FALSE,
                                           thinBoxFeature='utr',CDS='red',utr='white',
                                           transcriptAnnotation='symbol'
                     )%>%
                       # identity
                       {Gviz::displayPars(.)[names(orfcols)]<-orfcols;.}
                   ),
                   col.labels='black',
                   chromosome=as.character(seqnames(selgenerange))
  )
  #create barchart of intensities
  popViewport(1)
  pushViewport(vp2)
  cols = I(c(orfcols[names(which.min(orfscores))],orfcols[names(which.max(orfscores))]))
  grid.draw(lemon::g_legend(qplot(x=1:2,y=1:2,color=range(orfscores,na.rm=T))+
                              scale_color_gradient(name='Normalized ORF Expr\n(ORFs_pM)',
                                                   breaks = setNames(sort(na.omit(orfscores)),floor(na.omit(sort(orfscores)))%>% format(big.mark=",",scientific=FALSE) ),
                                                   low=cols[1],high=cols[2])+theme(text=element_text(size=14),legend.key.size=unit(.5,'inches'))
  ))
  popViewport(1)
  pushViewport(vp3)
  normalizePath(plotfile)
  #return file name
  return(plotfile)
} 
