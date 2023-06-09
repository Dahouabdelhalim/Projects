enrichment.function <- function (infile){
  IPR_all_table <- as.data.frame(table(IPR_gene_SNP$iprId))
  colnames(IPR_all_table)<- c("iprId", "freq_all")
  genes_of_interest_annotated<-unique(merge(infile, IPR_gene_SNP)[,c(1,ncol(merge(infile, IPR_gene_SNP)))])
  IPR_interest <- unique(genes_of_interest_annotated$iprId)
  genes_of_interest_table <- as.data.frame(table(genes_of_interest_annotated$iprId))
  colnames(genes_of_interest_table)<- c("iprId", "freq_cl1")
  genes_of_interest_all <- merge(genes_of_interest_table, IPR_all_table)
  for (x in IPR_interest){
    genes_of_interest_in_group=as.numeric(genes_of_interest_all[which(genes_of_interest_all$iprId==x),][,2])
    genes_of_interest_not_in_group=as.numeric(length(unique(genes_of_interest_annotated$gene)))-genes_of_interest_in_group
    genes_in_group_remainder = as.numeric(genes_of_interest_all[which(genes_of_interest_all$iprId==x),][,3]) - genes_of_interest_in_group
    genes_remainder = as.numeric(length(unique(IPR_gene_SNP$gene))) - as.numeric(genes_of_interest_all[which(genes_of_interest_all$iprId==x),][,3]) - genes_of_interest_not_in_group
    a = matrix(c(genes_of_interest_in_group, genes_in_group_remainder, genes_of_interest_not_in_group, genes_remainder), nrow=2, ncol=2)
    p = fisher.test(a, alternative='g')$p.value
    out_frame<- data.frame("iprId" = x, "Pvalue"= p)
    if(x==IPR_interest[1]){out_frame2<-out_frame}else{out_frame2<-rbind(out_frame2,out_frame)}
  }
  out_frame2 <- merge(out_frame2,  genes_of_interest_all)
  out_frame2<-out_frame2[order(out_frame2$Pvalue),]
  out_frame2$P_adjusted<-p.adjust(p =out_frame2$Pvalue, method = "BH")
  return(out_frame2)
}
