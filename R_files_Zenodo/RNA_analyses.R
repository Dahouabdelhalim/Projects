# Matt Teater
# 12/16/2022
# Analytical functions for Mlynarczyk et al. Science


require(igraph)
require(edgeR)
require(EDASeq)
require(RUVSeq)

midcolor <- function(color1,color2) {colorRampPalette(c(color1,color2))(3)[2]}

#Figure 3A
GSEA_network <- function(rnk) {
      originalDir <- getwd()
      #run CP gsea
      system(paste0("~/Downloads/GSEA_4.1.0/gsea-cli.sh GSEAPreranked -gmx ~/Downloads/c2.cp.v7.4.symbols.gmt -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnk ",rnk," -scoring_scheme weighted -rpt_label ",unlist(strsplit(basename(rnk),"\\\\."))[1]," -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ~/gsea_home/output/cp_vs_",unlist(strsplit(basename(rnk),"\\\\."))[1]))
      cpDir <- system(paste0("ls ~/gsea_home/output/cp_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1], "/"),intern=T)
      setwd(paste0("~/gsea_home/output/cp_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1], "/", cpDir))
      gseaResults <- table(system(paste0("ls ~/gsea_home/output/cp_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1], "/*/edb/results.edb"),intern=T), sep="\\t", stringsAsFactors=F, skip=3)
      cpResults <- t(sapply(gseaResults[,1], function(x) {c(unlist(strsplit(unlist(strsplit(x,"NES="))[2]," "))[1], unlist(strsplit(unlist(strsplit(x,"FDR="))[2]," "))[1])})); rownames(cpResults) <- sapply(gseaResults[,1], function(x) {unlist(strsplit(unlist(strsplit(x,"GENESET=gene_sets.gmt#"))[2]," "))[1]}); colnames(cpResults) <- c("NES", "FDR"); cpResults <- cpResults[complete.cases(cpResults),]
      if (length(rownames(cpResults)[as.numeric(cpResults[,2])<0.05])>1) {cpResults <- cpResults[as.numeric(cpResults[,2])<0.05,]; cpResults <- cpResults[order(as.numeric(cpResults[,1]),decreasing=T),]} else if (length(rownames(cpResults)[as.numeric(cpResults[,2])<0.05])==0) {cpResults = NULL} else {myPathway <- rownames(cpResults)[as.numeric(cpResults[,2])<0.05]; cpResults <- t(as.matrix(cpResults[as.numeric(cpResults[,2])<0.05,])); rownames(cpResults) <- myPathway}
      #perform CP leading edge analysis
      if (is.null(cpResults)) {cp.le=list()} else {
          system(paste0("~/Downloads/GSEA_4.1.0/gsea-cli.sh LeadingEdgeTool -dir ~/gsea_home/output/cp_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1], "/",cpDir,"/edb -gsets ", paste(rownames(cpResults),collapse=","), " -imgFormat png -zip_report false -enrichment_zip false -out cp_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1],"leadingEdge"))
          cp_conv <- table(system(paste0("ls ~/gsea_home/output/cp_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1],"/",cpDir,"/cp_vs_",unlist(strsplit(basename(rnk),"\\\\."))[1],"leadingEdge/*/conv.gct"),intern=T), sep="\\t", row.names=1, header=T, stringsAsFactors=F, skip=2); cp_conv <- cp_conv[,-1]
          cp.le = list()
          for (i in 1:nrow(cp_conv)) {
             cp.le <- c(cp.le, list(colnames(cp_conv[i,][,cp_conv[i,]>0])))
              names(cp.le) <- c(names(cp.le)[-length(names(cp.le))], rownames(cp_conv)[i])
          }
      }
      
      #run h gsea
      system(paste0("~/Downloads/GSEA_4.1.0/gsea-cli.sh GSEAPreranked -gmx ~/Downloads/h.all.v7.4.symbols.gmt -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnk ",rnk," -scoring_scheme weighted -rpt_label ",unlist(strsplit(basename(rnk),"\\\\."))[1]," -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ~/gsea_home/output/hallmark_vs_",unlist(strsplit(basename(rnk),"\\\\."))[1]))
      hDir <- system(paste0("ls ~/gsea_home/output/hallmark_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1], "/"),intern=T)
      setwd(paste0("~/gsea_home/output/hallmark_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1], "/", hDir))
      gseaResults <- table(system(paste0("ls ~/gsea_home/output/hallmark_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1], "/", hDir,"/edb/results.edb"),intern=T), sep="\\t", stringsAsFactors=F, skip=3)
      hResults <- t(sapply(gseaResults[,1], function(x) {c(unlist(strsplit(unlist(strsplit(x,"NES="))[2]," "))[1], unlist(strsplit(unlist(strsplit(x,"FDR="))[2]," "))[1])})); rownames(hResults) <- sapply(gseaResults[,1], function(x) {unlist(strsplit(unlist(strsplit(x,"GENESET=gene_sets.gmt#"))[2]," "))[1]}); colnames(hResults) <- c("NES", "FDR"); hResults <- hResults[complete.cases(hResults),]
      if (length(rownames(hResults)[as.numeric(hResults[,2])<0.05])>1) {hResults <- hResults[as.numeric(hResults[,2])<0.05,]; hResults <- hResults[order(as.numeric(hResults[,1]),decreasing=T),]} else if (length(rownames(hResults)[as.numeric(hResults[,2])<0.05])==0) {hResults = NULL} else {myPathway <- rownames(hResults)[as.numeric(hResults[,2])<0.05]; hResults <- t(as.matrix(hResults[as.numeric(hResults[,2])<0.05,])); rownames(hResults) <- myPathway}
      #perform h leading edge analysis
      if (nrow(hResults)<1) {h.le=list()} else {
          system(paste0("~/Downloads/GSEA_4.1.0/gsea-cli.sh LeadingEdgeTool -dir ~/gsea_home/output/hallmark_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1], "/",hDir,"/edb -gsets ", paste(rownames(hResults),collapse=","), " -imgFormat png -zip_report false -enrichment_zip false -out h_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1],"leadingEdge"))
          h_conv <- table(system(paste0("ls ~/gsea_home/output/hallmark_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1],"/",hDir,"/h_vs_",unlist(strsplit(basename(rnk),"\\\\."))[1],"leadingEdge/*/conv.gct"),intern=T), sep="\\t", row.names=1, header=T, stringsAsFactors=F, skip=2); h_conv <- h_conv[,-1]
          h.le = list()
          for (i in 1:nrow(h_conv)) {
             h.le <- c(h.le, list(colnames(h_conv[i,][,h_conv[i,]>0])))
              names(h.le) <- c(names(h.le)[-length(names(h.le))], rownames(h_conv)[i])
          }
      }
      
      #run melnicklab gsea
      system(paste0("~/Downloads/GSEA_4.1.0/gsea-cli.sh GSEAPreranked -gmx ~/Downloads/melnick_human_signature_revised717.gmt -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnk ",rnk," -scoring_scheme weighted -rpt_label ",unlist(strsplit(basename(rnk),"\\\\."))[1]," -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ~/gsea_home/output/melnicklab_vs_",unlist(strsplit(basename(rnk),"\\\\."))[1]))
      melnicklabDir <- system(paste0("ls ~/gsea_home/output/melnicklab_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1], "/"),intern=T)
      setwd(paste0("~/gsea_home/output/melnicklab_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1], "/", melnicklabDir))
      gseaResults <- table(system(paste0("ls ~/gsea_home/output/melnicklab_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1], "/*/edb/results.edb"),intern=T), sep="\\t", stringsAsFactors=F, skip=3)
      melnicklabResults <- t(sapply(gseaResults[,1], function(x) {c(unlist(strsplit(unlist(strsplit(x,"NES="))[2]," "))[1], unlist(strsplit(unlist(strsplit(x,"FDR="))[2]," "))[1])})); rownames(melnicklabResults) <- sapply(gseaResults[,1], function(x) {unlist(strsplit(unlist(strsplit(x,"GENESET=gene_sets.gmt#"))[2]," "))[1]}); colnames(melnicklabResults) <- c("NES", "FDR"); melnicklabResults <- melnicklabResults[complete.cases(melnicklabResults),]
      melnicklabResults <- melnicklabResults[melnicklabResults[,2]<0.05,]
      melnicklabResults <- melnicklabResults[order(as.numeric(melnicklabResults[,1]),decreasing=T),]
      #perform melnicklab leading edge analysis
      if (nrow(melnicklabResults)<1) {melnicklab.le=list()} else {
          system(paste0("~/Downloads/GSEA_4.1.0/gsea-cli.sh LeadingEdgeTool -dir ~/gsea_home/output/melnicklab_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1], "/",melnicklabDir,"/edb -gsets ", paste(rownames(melnicklabResults),collapse=","), " -imgFormat png -zip_report false -enrichment_zip false -out melnicklab_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1],"leadingEdge"))
          melnicklab_conv <- table(system(paste0("ls ~/gsea_home/output/melnicklab_vs_", unlist(strsplit(basename(rnk),"\\\\."))[1],"/",melnicklabDir,"/melnicklab_vs_",unlist(strsplit(basename(rnk),"\\\\."))[1],"leadingEdge/*/conv.gct"),intern=T), sep="\\t", row.names=1, header=T, stringsAsFactors=F, skip=2); melnicklab_conv <- melnicklab_conv[,-1]
          melnicklab.le = list()
          for (i in 1:nrow(melnicklab_conv)) {
             melnicklab.le <- c(melnicklab.le, list(colnames(melnicklab_conv[i,][,melnicklab_conv[i,]>0])))
              names(melnicklab.le) <- c(names(melnicklab.le)[-length(names(melnicklab.le))], rownames(melnicklab_conv)[i])
          }
      }
      
      #make network
      all.list <- c(cp.le, h.le, melnicklab.le)
      myResults <- rbind(cpResults,hResults, melnicklabResults)
      #all.list <- all.list[names(all.list)%in%paste0(rownames(myResults),"_signal")]
      jaccard <- matrix(NA, nrow=length(all.list), ncol=length(all.list))
      for (i in 1:length(all.list)) {
          for (j in 1:length(all.list)) {
              jaccard[i,j] <- length(intersect(all.list[[i]], all.list[[j]]))/length(union(all.list[[i]], all.list[[j]]))
          }
      }
      rownames(jaccard) <- names(all.list); colnames(jaccard) <- names(all.list)
      jaccard[jaccard==1] <- 0

      net <- graph.adjacency(jaccard, weighted=TRUE)
      V(net)$color <- sapply(as.numeric(myResults[sapply(rownames(jaccard),function(x) {unlist(strsplit(x,"_signal"))[1]}),1]), function(x) {if (x<(-2)) {"royalblue1"} else if (x<(-1.75) & x>=(-2)) {"royalblue2"} else if (x<(-1.5) & x>=(-1.75)) {"royalblue3"} else if (x<(-1.25) & x>=(-1.5)) {"royalblue4"} else if (x<1.25 & x>=(-1.25)) {"black"} else if (x<1.5 & x>=1.25) {"firebrick4"} else if (x<1.75 & x>=1.5) {"firebrick3"} else if (x<2 & x>=1.75) {"firebrick2"} else if (x>=2) {"firebrick1"}})
      V(net)$size <- (igraph::degree(net, mode="all")+10)/10
      V(net)$label <- rownames(jaccard)
      V(net)$label.cex <- 0.4
      V(net)$label.dist <- 1
      E(net)$width <- E(net)$weight/6
      E(net)$arrow.size <- 0
      E(net)$edge.color <- "gray80"
      plot(net)
      title(unlist(strsplit(basename(rnk),"\\\\."))[1])
      print(rownames(jaccard))
      
      myResults <- rbind(cpResults, hResults, melnicklabResults)
      myResults <- cbind(myResults, leadingEdge=sapply(rownames(myResults), function(x) {paste(unlist(all.list[names(all.list)%in%paste0(x,"_signal")]),collapse=",")}))
      return(myResults)
   }

#Figure 3B, S8G, S8K, S8M, S8R, S9E, S15J
log2TMM <- function(controlIDs, experimentalIDs, readCountMatrix) {
    myData <- readCountMatrix[,c(controlIDs, experimentalIDs)]
    experimentCategory <- factor(c(rep("Control",length(controlIDs)), rep("Experimental",length(experimentalIDs))))
    design <- model.matrix(~experimentCategory)
    y <- DGEList(counts=myData, genes=rownames(myData))
    y <- calcNormFactors(y, method='TMM') # Normalize library sizes using TMM
    y <- estimateGLMCommonDisp(y,design)
    y <- estimateGLMTrendedDisp(y,design)
    y <- estimateGLMTagwiseDisp(y,design)
    fit <- glmFit(y,design)
    
    lrt_celltype <- glmLRT(fit,coef=2)
    myLog2 <- lrt_celltype$table[,1]; names(myLog2) <- rownames(lrt_celltype$table)
    return(as.matrix(myLog2))
}
log2TPM <- function(controlIDs, experimentalIDs, tpmMatrix) {
    myLog2 <- log2((rowMeans(tpmMatrix[,experimentalIDs])+0.1)/(rowMeans(tpmMatrix[,controlIDs])+0.1)
    return(as.matrix(myLog2))
}
preRankedGSEA_plot <- function(gmx, rnk) {
    originalDir <- getwd()
    setwd("~/GSEA_input/")
    if (!(paste(deparse(substitute(gmx)),".gmx",sep="")%in%list.files())) {write.table(gmx, file=paste("~/GSEA_input/",deparse(substitute(gmx)),".gmx",sep=""), row.names=F, col.names=F, quote=F, sep="\\t")} else {print(paste(deparse(substitute(gmx)),".gmx previously written",sep=""))}
    if (!(paste(deparse(substitute(rnk)),".rnk",sep="")%in%list.files())) {write.table(rnk, file=paste("~/GSEA_input/",deparse(substitute(rnk)),".rnk",sep=""), row.names=T, col.names=F, quote=F, sep="\\t")} else {print(paste(deparse(substitute(rnk)),".rnk previously written",sep=""))}
    if (paste(deparse(substitute(gmx)), "_VS_", deparse(substitute(rnk)),sep="")%in%list.files("~/gsea_home/output/")) {system(paste("rm -rf ~/gsea_home/output/", deparse(substitute(gmx)), "_VS_", deparse(substitute(rnk)),sep=""))}
    system(paste("java -cp ~/PROGRAMS/gsea2-2.0.13.jar -Xmx8g xtools.gsea.GseaPreranked -gmx ~/GSEA_input/", deparse(substitute(gmx)), ".gmx -collapse false -mode Max_probe -norm meandiv -nperm 10000 -rnk ~/GSEA_input/", deparse(substitute(rnk)), ".rnk -scoring_scheme weighted -rpt_label ", deparse(substitute(gmx)), "_vs_", deparse(substitute(rnk))," -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500000 -set_min 0 -zip_report false -out ~/gsea_home/output/", deparse(substitute(gmx)), "_VS_", deparse(substitute(rnk))," -gui false", sep=""))
    setwd(paste("~/gsea_home/output/", deparse(substitute(gmx)), "_VS_", deparse(substitute(rnk)),sep=""))
    setwd(list.files())
    myStats <- read.table("edb/results.edb",sep=" ",fill=T, stringsAsFactors=F)
    myGENESETS <- sapply(myStats[sapply(myStats$V6, function(x) {substr(x,1,7)})%in%"GENESET",6], function(x) {unlist(strsplit(x,"GENESET=gene_sets.gmt#"))[2]})
    pdf(paste(deparse(substitute(gmx)), "_VS_", deparse(substitute(rnk)),"_GSEAplots.pdf",sep=""),width=5,height=6)
    for (i in 1:length(myGENESETS)) {
        myES <- read.table(paste(myGENESETS[i],".xls",sep=""),header=T, sep="\\t", stringsAsFactors=F)
        layout(matrix(c(1,1,2,2,3,3), ncol=2, byrow = TRUE), heights=c(10,1,9))
        par(mar=c(0, 4, 4, 2) + 0.1)
        NES=myStats[myStats$V6%in%names(myGENESETS[i]),8]; NES.sig <- signif(as.numeric(unlist(strsplit(NES,"="))[2]),3)
        print(paste("NES=",NES.sig,sep=""))
        FDR=myStats[myStats$V6%in%names(myGENESETS[i]),10]; FDR.sig <- signif(as.numeric(unlist(strsplit(FDR,"="))[2]),3)
        print(paste("FDR=", FDR.sig,sep=""))
        plot(myES[,5], myES[,7], xlab=NA, xaxt='n', ylab="Enrichment score (ES)", type="l", bty="n", col="springgreen3",lwd=3, main=paste(myGENESETS[i],"\\n n=",length(setdiff(gmx[,toupper(gmx[1,])%in%myGENESETS[i]][complete.cases(gmx[,toupper(gmx[1,])%in%myGENESETS[i]])],""))-1,sep=""))
        abline(h=0,col="black")
        if (FDR.sig==0) {legend("topright", legend=c(paste("NES=",NES.sig,sep=""),"FDR<0.001"), bty="n", fill=NA, border=NA,cex=1.4)} else {legend("topright", legend=c(paste("NES=",NES.sig,sep=""),paste("FDR=", FDR.sig,sep="")), bty="n", fill=NA, border=NA,cex=1.4)}
        myRNK <- cbind(1:length(rnk),rnk[order(rnk,decreasing=T)])
        myGMX <- names(rnk[order(rnk,decreasing=T),]); myGMX[!(myGMX%in%myES$PROBE)] <- NA; myGMX <- cbind(1:length(rnk),myGMX)
        par(mar=c(0, 4, 0, 2) + 0.1)
        plot(as.numeric(myGMX[complete.cases(myGMX),1]),rep(1,length(myGMX[complete.cases(myGMX),1])),xlim=c(1,nrow(myRNK)), type="h", xaxt='n', xlab=NA, yaxt='n', ylab=NA, bty='n', ylim=c(0,1), col="#00000033")
        Pos95 <- max(myRNK[myRNK[,2]>=quantile(myRNK[,2],probs=0.95),1])
        topQuartilePos <- max(myRNK[myRNK[,2]>=quantile(myRNK[,2],probs=0.75),1])
        bottomQuarilePos <- max(myRNK[myRNK[,2]>=quantile(myRNK[,2],probs=0.25),1])
        Pos5 <- max(myRNK[myRNK[,2]>=quantile(myRNK[,2],probs=0.05),1])
        polygon(c(rep(1,2),rep(Pos95,2)), c(0,0.25,0.25,0), col=paste(col2hex("red4"),"AA",sep=""), border=paste(col2hex("red4"),"AA",sep=""))
        polygon(c(rep(Pos95,2),rep(topQuartilePos,2)), c(0,0.25,0.25,0), col=paste(col2hex("red"),"AA",sep=""), border=paste(col2hex("red"),"AA",sep=""))
        polygon(c(rep(topQuartilePos,2),rep(nrow(myRNK[myRNK[,2]>0,]),2)), c(0,0.25,0.25,0), col=paste(col2hex("pink"),"AA",sep=""), border=paste(col2hex("pink"),"AA",sep=""))
        polygon(c(rep(nrow(myRNK[myRNK[,2]>0,]),2),rep(bottomQuarilePos,2)), c(0,0.25,0.25,0), col=paste(col2hex("cornflowerblue"),"AA",sep=""), border=paste(col2hex("cornflowerblue"),"AA",sep=""))
        polygon(c(rep(bottomQuarilePos,2),rep(Pos5,2)), c(0,0.25,0.25,0), col=paste(col2hex("blue"),"AA",sep=""), border=paste(col2hex("blue"),"AA",sep=""))
        polygon(c(rep(Pos5,2),rep(nrow(myRNK),2)), c(0,0.25,0.25,0), col=paste(col2hex("blue4"),"AA",sep=""), border=paste(col2hex("blue4"),"AA",sep=""))
        par(mar=c(5, 4, 0, 2) + 0.1)
        plot(myRNK, xlab=paste("<- log2(",deparse(substitute(rnk)),")",sep=""), ylab="Ranked list metric (PreRanked)", type="h", bty="n", col="grey", xaxt='n')
        axis(1,at=c(1,5000,10000,15000,20000,nrow(myRNK)))
        abline(v=nrow(myRNK[myRNK[,2]>0,]), lty=2, col="grey")
        text(nrow(myRNK[myRNK[,2]>0,]),0, labels=paste("Zero cross at ",nrow(myRNK[myRNK[,2]>0,]),sep=""), pos=3, cex=0.8)
    }
    dev.off()
    setwd(originalDir)
}

#Figure S8L
centrimoMotif_plot <- function(outputDIR, background.fasta, test.fasta, motif.meme, evalThreshold, fisherThreshold) {
    numberSequences <- system(paste0("grep '>' ", test.fasta," | wc -l"),intern=T)
    system(paste0("centrimo --oc ", outputDIR, " --verbosity 1 --score 5.0 --ethresh ", numberSequences," --neg ", background.fasta, " ", test.fasta, " ", motif.meme))
    centrimo <- read.table(pastte0(outputDIR,"/centrimo.txt"), sep="\\t", header=T, stringsAsFactors=F, row.names=NULL); colnames(centrimo) <- colnames(centrimo)[c(2:ncol(centrimo),1)]; centrimo <- centrimo[,-ncol(centrimo)]
    HOCOMOCO <- read.table("http://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv", sep="\\t", header=T, row.names=1, stringsAsFactors=F)
    myC <- cbind(-log10(centrimo$E.value), -log10(centrimo$fisher_adj_pvalue)); rownames(myC) <- sapply(centrimo[,2], function(x) {unlist(strsplit(x,"\\\\."))[1]})
    plot(myC, xlab="-log10(E-value)", ylab="-log10(fisher adj_pvalue)", pch=19)
    text(myC[myC[,2]>-log10(evalThreshold),], labels=sapply(rownames(myC[myC[,2]>-log10(fisherThreshold),]), function(x) {unlist(strsplit(x,"_"))[1]}), pos=c(3,4,1,4,3,3,3,3,3,3), cex=0.7)
    text(myC[myC[,1]>-log10(evalThreshold) & myC[,2]>-log10(fisherThreshold) & myC[,2]<(-log10(fisherThreshold)),], labels=sapply(rownames(myC[myC[,1]>-log10(evalThreshold) & myC[,2]>-log10(fisherThreshold) & myC[,2]<(-log10(fisherThreshold)),]), function(x) {unlist(strsplit(x,"_"))[1]}), pos=c(2,4,4,4,1,4), cex=0.7)
    legend("right", legend=paste0("E-value=",evalThreshold), bty="n", cex=0.7, text.col="blue")
    legend("topleft", legend=paste("Fisher",paste0("adj_p=",fisherThreshold),sep="\\n"), bty="n", cex=0.7, text.col="blue")
    abline(v=-log10(evalThreshold), col="blue"); abline(h=-log10(fisherThreshold), col="blue")
}

#Figure 3H, 4B, S9G
hypergeom_pval <-  function(genesetA, genesetB, totalNumberGenes) {phyper(length(genesetA[genesetA%in%genesetB])-1, length(genesetA), totalNumberGenes-length(genesetA), length(genesetB), lower.tail=FALSE)}

#Figure 4A
differentialAbundance <- function(refSeqReadCountMatrix, ercc92ReadCountMatrix, genoFactors, sourceFactors, FC, qval) {
    #filter out non-expressed genes
    genes <- refSeqReadCountMatrix[apply(refSeqReadCountMatrix, 1, function(x) length(x[x>5])>=2),]
    spikes <- ercc92ReadCountMatrix[apply(ercc92ReadCountMatrix, 1, function(x) length(x[x>5])>=2),]
    filtered <- rbind(genes, spikes)
    #genoFactors {EGFP, W, Q36H}
    #sourceFactors  {input, IP}
    set <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(genoFactors, row.names=colnames(filtered)))
    setNorm <- RUVg(set, rownames(spikes), k=1)
    design <- model.matrix(~ genoFactors - sourceFactors + W_1, data=pData(setNorm))
    y <- DGEList(counts=counts(setNorm), group=x)
    y <- calcNormFactors(y, method="upperquartile")
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmQLFit(y, design)
    qlf.ip <- glmQLFTest(fit, coef=2); results.ip <- qlf.ip$table; results.ip$padj <- p.adjust(results.ip[,4], method="BH")
    ip.up <- rownames(results.ip[results.ip[,1]>log2(FC) & results.ip[,5]<qval,])
    ip.down <- rownames(results.ip[results.ip[,1]<log2(1/FC) & results.ip[,5]<qval,])
    return(list(ip.up, ip.down)
}

#Figure S8F, S9G
differentialExpression <- function(controlIDs, experimentalIDs, readCountMatrix, FC, qval) {
    myData <- readCountMatrix[,c(controlIDs, experimentalIDs)]
    experimentCategory <- factor(c(rep("Control",length(controlIDs)), rep("Experimental",length(experimentalIDs))))
    design <- model.matrix(~experimentCategory)
    y <- DGEList(counts=myData, genes=rownames(myData))
    y <- estimateGLMCommonDisp(y,design)
    y <- estimateGLMTrendedDisp(y,design)
    y <- estimateGLMTagwiseDisp(y,design)
    fit <- glmFit(y,design)
    
    lrt_celltype <- glmLRT(fit,coef=2)
    lrt_celltype$table=cbind(lrt_celltype$table, p_adj=p.adjust(lrt_celltype$table[,4], method="BH"))
    return(list(rownames(lrt_celltype$table)[lrt_celltype$table[,5]<qval & lrt_celltype$table[,1]>log2(FC)],rownames(lrt_celltype$table)[lrt_celltype$table[,5]<qval & lrt_celltype$table[,1]<log2(1/FC)], lrt_celltype$table))
}

#Figure S8H
gsvaAnalysis <- function(normExpr, signatureList, controlIDs, testIDs) {
    p <- nrow(normExpr)   ## number of genes
    n <- ncol(normExpr)   ## number of samples
    nGS <- 3 ## number of gene sets
    min.sz <- 10  ## minimum gene set size
    max.sz <- 100 ## maximum gene set size
    gs <- as.list(sample(min.sz:max.sz, size=nGS, replace=TRUE)) ## sample gene set sizes
    gs <- lapply(gs, function(n, p) sample(1:p, size=n, replace=FALSE), p) ## sample gene sets
    es.dif <- gsva(as.matrix(normExpr), signatureList, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
    
    for (i in 1:nrow(es.dif)) {
        d.test <- density(es.dif[i,testIDs], adjust=1)
        d.control <- density(es.dif[i,controlIDs], adjust=1)
        plot(c(-0.75,d.test$x[d.test$x>=-0.75 & d.test$x<=0.75],0.75), c(0,d.test$y[d.test$x>=-0.75 & d.test$x<=0.75],0), ylab="Density", xlab=paste0(rownames(es.dif)[i], " GSVA score"), main=paste0("p=", signif(wilcox.test(es.dif[i,testIDs], es.dif[i,controlIDs])$p.value,3)), col="darkorange", xlim=c(-0.75,0.75), ylim=c(min(c(d.test$y, d.control$y)), max(c(d.test$y, d.control$y))), type="l")
        polygon(c(-0.75,d.test$x[d.test$x>=-0.75 & d.test$x<=0.75],0.75), c(0,d.test$y[d.test$x>=-0.75 & d.test$x<=0.75],0), col=paste0(col2hex("darkorange"),"AA"), border="darkorange")
        lines(c(-0.75,d.control$x[d.control$x>=-0.75 & d.control$x<=0.75],0.75), c(0,d.control$y[d.control$x>=-0.75 & d.control$x<=0.75],0), col="dodgerblue")
        polygon(c(-0.75,d.control$x[d.control$x>=-0.75 & d.control$x<=0.75],0.75), c(0,d.control$y[d.control$x>=-0.75 & d.control$x<=0.75],0), col=paste0(col2hex("dodgerblue"),"AA"), border="dodgerblue")
    }
}
