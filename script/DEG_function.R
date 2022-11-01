my_DEG_analysis <-function(gse_name,phenodata,ballgownfile,smallprotein,logfc,pvalue)
{

	library(ballgown)
	#library(genefilter)
	#library(dplyr)
	#library(devtools)
	library(ggplot2)
	#library(GSEABase)
  

	###指定分组信息
	sampleGroup <- read.csv(phenodata, header = TRUE)
	#DESeq2说明清楚哪个因子level对应的control，避免后面解释数据遇到麻烦（你不知道到底logFC是以谁做参照的，希望是以control作为参照，这样logFC>1就表示实验组表达大于对照组；但是R不知道，R默认按字母顺序排列因子顺序，所以要对这个因子顺序set一下：
	#factor levels，写在前面的level作为参照
	sampleGroup$Treatment <- factor(sampleGroup$Treatment, levels = c("control", "case"))

	###读入数据
	bg_chrX=ballgown(dataDir =ballgownfile,samplePattern = "SRR", meas='all',pData=sampleGroup)
	whole_tx_table=texpr(bg_chrX, 'all')
	transcript_fpkm=texpr(bg_chrX, 'FPKM')
	colnames(transcript_fpkm)<-substring(colnames(transcript_fpkm), 6)
	rownames(transcript_fpkm)<-whole_tx_table$gene_name
	
	###提取fpkm矩阵
	data<-transcript_fpkm
	group_list <- sampleGroup$Treatment
	
	expMatrix <- data
	fpkmToTpm <- function(fpkm)
	{
	  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
	}
	tpms <- apply(expMatrix,2,fpkmToTpm)
	tpms[1:3,]
	colSums(tpms)
	exprSet <- tpms
	
	#limma 对数据进行归一化处理
	library(limma) 
	exprSet2=normalizeBetweenArrays(exprSet)
	#boxplot(exprSet2,outline=FALSE, notch=T,col=group_list, las=2)
	#判断数据是否需要转换
	exprSet3 <- log2(exprSet2+1)
	
	###主成分分析图
	dat <- t(exprSet2)#画PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换
	dat <- as.data.frame(dat)#将matrix转换为data.frame
	dat <- cbind(dat,group_list) #cbind横向追加，即将分组信息追加到最后一列
	#dat[1:5,1:5]
	
	#BiocManager::install("FactoMineR")
	#BiocManager::install("factoextra")
	library("FactoMineR")
	library("factoextra") 
	
	# before PCA analysis
	dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的
	
	fviz_pca_ind(dat.pca,
	             geom.ind = "point", # show points only (nbut not "text")
	             col.ind = dat$group_list, # color by groups
	             # palette = c("#00AFBB", "#E7B800"),
	             addEllipses = TRUE, # Concentration ellipses
	             legend.title = "Groups"
	)
	ggsave(file=paste(gse_name,"all_samples_PCA.png"))
	
	
	###差异分析
	dat <- exprSet3
	#design <- model.matrix(~0+factor(group_list))
	design=model.matrix(~factor( group_list ))
	fit=lmFit(dat,design)
	fit=eBayes(fit)
	options(digits = 4)
	topTable(fit,coef=2,adjust='BH')
	degene=topTable(fit,coef=2,adjust='BH',number = Inf)
	head(degene) 
	write.csv(degene,paste(gse_name,"gene_results.csv") ,row.names=FALSE)
	
	###筛选差异表达小蛋白
	library(ggrepel)
	smallp<-read.csv(smallprotein,header = T)
	degsp<-subset(degene,degene$ID %in% smallp$Gene | degene$ID %in% smallp$names)
	degsp$g=ifelse(degsp$P.Value<pvalue & abs(degsp$logFC) >logfc,
	               ifelse(degsp$P.Value<pvalue & degsp$logFC > logfc,'UP','DOWN'),'STABLE') 
	deg=degsp
	colnames(deg)<- c('ID','Log2FC','AveExpr','T.statistic','P.value','Adj.P.val','Log.odds','Group')
	deg$Log2FC = round(deg$Log2FC,4)
	deg$AveExpr = round(deg$Log2FC,4)
	deg$T.statistic = round(deg$T.statistic,4)
	deg$P.value = round(deg$P.value,4)
	deg$Adj.P.val = round(deg$Adj.P.val,4)
	deg$Log.odds = round(deg$Log.odds,4)
	write.csv(deg, paste(gse_name,"microprotein_results.csv"),row.names=FALSE)
	  
	table(degsp$g)
	data<-degsp
	data$threshold = as.factor(degsp$g)
	data$label <- ifelse(data$P.Value < pvalue & abs(data$logFC) >= logfc,data$ID,"")
	p<-ggplot(data, aes(logFC, -log10(P.Value),col = threshold)) +
	  ggtitle("Differential genes") +
	  geom_point(alpha=0.3, size=2) +
	  scale_color_manual(values=c("blue", "grey","red")) +
	  labs(x="logFC",y="-log10 (p-value)") +
	  theme_bw()+
	  geom_hline(yintercept = -log10(as.numeric(pvalue)), lty=4,col="grey",lwd=0.6) +
	  geom_vline(xintercept = c(-as.numeric(logfc), as.numeric(logfc)), lty=4,col="grey",lwd=0.6) +
	  theme(plot.title = element_text(hjust = 0.5),
	        panel.grid=element_blank(),
	        axis.title = element_text(size = 12),
	        axis.text = element_text(size = 12))
	p
	p+geom_text_repel(data = data, aes(x = logFC, y = -log10(P.Value), label = label),
	                  size = 3,box.padding = unit(0.5, "lines"),
	                  point.padding = unit(0.8, "lines"),
	                  segment.color = "black",
	                  show.legend = FALSE)
	
	ggsave(file=paste(gse_name,"VolcanoMP.png"),width = 7, height = 7)
	
	
	
	###差异基因注释分析
	#library(ggstatsplot)
	#library(cowplot)
	library(clusterProfiler)
	library(enrichplot)
	#library(ReactomePA)
	#library(stringr)
	#library(tidyr)
	library(org.Sc.sgd.db)
	
	degene$ENSEMBL=degene$ID
	df <- bitr(unique(degene$ENSEMBL), fromType = "ENSEMBL", 
	           toType = c("ENTREZID","GENENAME"),
	           OrgDb = org.Sc.sgd.db)
	#head(df)
	DEG=degene
	#head(DEG)
	DEG$g=ifelse(DEG$P.Value<pvalue & abs(DEG$logFC) >logfc,
	             ifelse( DEG$P.Value<pvalue & DEG$logFC > logfc,'UP','DOWN'),'STABLE') 
	
	DEG=merge(DEG,df,by='ENSEMBL')
	head(DEG)
	
	save(DEG,file = 'anno_DEG.Rdata')
	DEG_diff=DEG[DEG$g == 'UP' | DEG$g == 'DOWN',] 
	gene_diff=DEG_diff$ENSEMBL
	
	###通路与基因之间的关系可视化
	###制作genlist三部曲：
	## 1.获取基因logFC
	geneList <- as.numeric(DEG$logFC)
	## 2.命名
	names(geneList) = as.character(DEG$ENSEMBL)
	## 3.排序很重要
	geneList = sort(geneList, decreasing = TRUE)
	geneList = na.omit(geneList)
	head(geneList)


	geneListgo <- as.numeric(DEG$logFC)
	names(geneListgo) = as.character(DEG$ENTREZID)
	geneListgo = sort(geneListgo, decreasing = T)
	geneListgo <- geneListgo[!is.na(geneListgo)]
	head(geneListgo)

	
	###GO分析
	#library(DOSE)
	#library(ggnewscale)
	#library(topGO)
	GO<-enrichGO(DEG_diff$ENSEMBL, OrgDb = "org.Sc.sgd.db", keyType = "ENSEMBL",ont = "all", pvalueCutoff = 0.5, 
	             pAdjustMethod = "BH", qvalueCutoff = 0.5, minGSSize = 10, 
	             maxGSSize = 500, readable = FALSE, pool = FALSE)
	write.csv(GO, paste(gse_name,"GO_results.csv"),row.names=FALSE)
	if (length(GO$ID) > 0){
		dotplot(GO, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
		ggsave(file=paste(gse_name,"GO_dotplot.png"))
		barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")+labs(y = "GeneNumber", x = NULL)
		ggsave(file=paste(gse_name,"GO_barplot.png"))

		enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE)
		ggsave(file=paste(gse_name,"GO_gene_pathway.png"))
	
	
		enrichplot::heatplot(GO,foldChange=geneListgo,showCategory = 50)
		ggsave(file=paste(gse_name,"GO_gene_pathway_heatmap.png"))
	}
	
	
	
	###KEGG 
	enrichKK <- enrichKEGG(gene   =   as.character(gene_diff),
	                       organism  = 'sce',
	                       keyType = "kegg",
	                       #universe     = gene_all,
	                       pvalueCutoff = 0.5,
	                       qvalueCutoff = 0.5)
	write.csv(enrichKK, paste(gse_name,"KEGG_results.csv"),row.names=FALSE)
	print(enrichKK)
	#head(enrichKK)[,1:6] 
	#气泡图
	if (length(enrichKK$ID) > 0){
		dotplot(enrichKK)
		ggsave(file=paste(gse_name,"KEGG_dotplot.png"))
		##最基础的条形图和点图
		#条带图
		barplot(enrichKK,showCategory=20)+labs(y = "GeneNumber", x = NULL)
		ggsave(file=paste(gse_name,"KEGG_barplot.png"))
		
		enrichplot::cnetplot(enrichKK,circular=FALSE,colorEdge = TRUE)#circluar为指定是否环化，基因过多时建议设置为FALSE
		ggsave(file=paste(gse_name,"KEGG_gene_pathway.png"))

		enrichplot::heatplot(enrichKK,foldChange=geneList,showCategory = 50)
		ggsave(file=paste(gse_name,"KEGG_gene_pathway_heatmap.png"))
	}
	
	GSEA_KEGG <- gseKEGG(geneList, organism = 'sce', nPerm = 1000, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1)
	if (length(GSEA_KEGG$ID) > 0){
 		ridgeplot(GSEA_KEGG)
 		ggsave(file=paste(gse_name,"enrichKEGG_ridgeplot.png"))
		if(length(GSEA_KEGG$ID) < 5){
		enrichplot::gseaplot2(GSEA_KEGG,1:as.numeric(dim(GSEA_KEGG)[1]))
		ggsave(file=paste(gse_name,"enrichKEGG_gseaplot.png"))
		}else{
		enrichplot::gseaplot2(GSEA_KEGG,1:5)
		ggsave(file=paste(gse_name,"enrichKEGG_gseaplot.png"))
		}
	}
	
	GSEA_GO <-  gseGO(geneList     = geneListgo ,
	                OrgDb        = org.Sc.sgd.db,
	                 keyType      = "ENTREZID",
	                 ont          = "all",
	                 nPerm        = 1000,   ## 排列数
	                 minGSSize    = 5,
	                 maxGSSize    = 500,
	                 pvalueCutoff = 0.95,
	                 verbose      = TRUE)
	print(GSEA_GO)
	if (length(GSEA_GO$ID) > 0){
		ridgeplot(GSEA_GO) 
		ggsave(file=paste(gse_name,"enrichGO_ridgeplot.png"))
		if(length(GSEA_GO$ID) < 5){
		enrichplot::gseaplot2(GSEA_GO,1:as.numeric(dim(GSEA_GO)[1]))
		ggsave(file=paste(gse_name,"enrichGO_gseaplot.png"))
		}else{
		enrichplot::gseaplot2(GSEA_GO,1:5)
		ggsave(file=paste(gse_name,"enrichGO_gseaplot.png"))
		}
		
	}
	
}
#my_DEG_analysis("GSE63516","phenodata.csv","ballgown","smallprotein.csv",1.5,0.05)
args = commandArgs(trailingOnly=TRUE)
gse_name<-args[1]
phenodata<-args[2]
ballgownfile<-args[3]
smallprotein<-args[4]
logfc<-args[5]
pvalue<-args[6]

my_DEG_analysis(gse_name,phenodata,ballgownfile,smallprotein,logfc,pvalue)

