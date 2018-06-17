source("http://bioconductor.org/biocLite.R")
# source("/Volumes/Peter/R_functions/Peters_biomart.R")
library(easybiomart)
library(data.table)
library(randomForestSRC)
library(rpart)
# library(pheatmap) 
# library(ChIPseeker)
# library(clusterProfiler)
# library(methyAnalysis)
library(foreach)
library(doParallel)
# library(clusterProfiler)
# library(methylumi)
# library(cummeRbund)
# suppressMessages( library( ExpressionAtlas ) )
# library(ArrayExpress)
# library(RISmed)
library(ggthemes)
# library(miRLAB)
# library(seqinr)
# library(GenomicFeatures)
# library(Biostrings)
# library(mixOmics)
# library(ChIPpeakAnno)
# library(miRNAtap.db)
# library(miRNAtap)
library(topGO)
library(org.Hs.eg.db)
# library(Gviz)
# library(splitstackshape)
# library(Sushi)
# library(survival)
library(BiocParallel)
library(parallel)
# library(grid)
# library(gtable)
# library(GSA)
# library(wordcloud)
# library(xtable)
library(beeswarm)
# library(tm)
# library(PerformanceAnalytics)
# library(vioplot)
library(edgeR)
# library(preprocessCore)
library(ggrepel)
# library(gridExtra)
# library(GenomicRanges)
# library(vioplot)
# library(d3heatmap)
# library(mclust)
# library(rgl)
# library(gplots)
# library(rtracklayer)
# library(igraph)
# require(ReactomePA)
library(plotly)
library(ggplot2)
# library(ComplexHeatmap)
# library(DEXSeq)

#Define custom dist and hclust functions for use with heatmaps
#distance functions for clustering samples
euclidian.dist =function(c) {dist(c,method="euclidian")}
corr.dist = function(x) as.dist(1-cor(t(x)))
#distance functions for clustering genes
avg.clust = function(c) {hclust(c,method="average")}
ward.D.clust = function(x) hclust(x,method="ward.D")

write.delim <- write.table
formals(write.delim)$sep <- "\t"
formals(write.delim)$quote <- F
formals(write.delim)$row.names <- F



ggsave.a4 <- function(file) {
    library(tools)
    ggsave(paste(file_path_sans_ext(file), "pdf", sep = "."), device = "pdf", width = 210, height = 297, units = "mm", dpi = 300)
    ggsave(paste(file_path_sans_ext(file), "png", sep = "."), device = "png", width = 210, height = 297, units = "mm", dpi = 300)
}

ggsave.a4r <- function(file) {
    library(tools)
    ggsave(paste(file_path_sans_ext(file), "pdf", sep = "."), device = "pdf", width = 297, height = 210, units = "mm", dpi = 300)
    ggsave(paste(file_path_sans_ext(file), "png", sep = "."), device = "png", width = 297, height = 210, units = "mm", dpi = 300)
}

ggsave.html <- function(ggplotly.object, file) {
    htmlwidgets::saveWidget(ggplotly.object, file, selfcontained = T, libdir = NULL)
}

# set paper size for quartz
A3_landscape <- quartz.options(width = 16.53, height= 11.69)
A3_portrait <- quartz.options(width = 11.69, height = 16.53)
A4_landscape <- quartz.options(width = 11.69, height = 8.27)
A4_portrait <-  quartz.options(width = 8.27, height = 11.69)

#awk '$3 == "gene"' Homo_sapiens.GRCh38.82.gtf | awk -F "\t" '{print $9}' - | awk -F ";" '{print $1, $3, $5}' - | awk -F " " '{print $2, $4, $6}' - | sed -e 's/\"//g' > gene_annotation.txt
# grch38.file <- read.delim("/Volumes/Peter/gene_models/gene_annotation.txt", header=F, sep=" ")


vioplot_mna <- function(HUGO_gene_name) {
    gene_name <- grch38.file[match(HUGO_gene_name, table = grch38.file$V2), 1]
    gene <- subset(MYCN_CPM_log, row.names(MYCN_CPM_log) == gene_name)
    gene_t <- data.frame(t(gene))
    MYCN_samples <- design[1:2]
    MYCN_amp <- subset(MYCN_samples, MNA == 1)
    MYCN_non_amp <- subset(MYCN_samples, MNA == 0)
    gene_MNA <- subset(gene_t, row.names(gene_t) %in% MYCN_amp$Sample)
    gene_nonMNA <- subset(gene_t, row.names(gene_t) %in% MYCN_non_amp$Sample)

    pval <- signif(MYCNdata.ordered[match(HUGO_gene_name, table = MYCNdata.ordered$external_gene_name), 5],3)

    logFC <- signif(MYCNdata.ordered[match(HUGO_gene_name, table = MYCNdata.ordered$external_gene_name), 2],3)
    FC <- if (logFC < 0) {
        -1/(2^logFC)
    } else {
        2^logFC
    }
    FC <- signif(FC, 3)

    logCPM <- signif(MYCNdata.ordered[match(HUGO_gene_name, table = MYCNdata.ordered$external_gene_name), 3],3)
    CPM <- if (logCPM < 0) {
        -1/(2^logFC)
    } else {
        2^logCPM
    }
    CPM <- round(CPM, 0)

    #get R value
    CPM.an2 <- log2((CPM.an[,4:ncol(CPM.an)] + 1))
    row.names(CPM.an2) <- make.names(CPM.an$V2, unique = T)
    MYCN_vs_gene <- data.frame(t(subset(CPM.an2, row.names(CPM.an2) == HUGO_gene_name | row.names(CPM.an2) == "MYCN")))[4:nrow(CPM.an),]
    MYCN_vs_gene <- MYCN_vs_gene[complete.cases(MYCN_vs_gene),]
    MYCN_vs_gene2 <- merge(MYCN_vs_gene, design[1:2], by.x="row.names", by.y="Sample", all.x=T)
    R_val <- signif(cor(MYCN_vs_gene2[,2:3])[1,2],3)

    quartz(paste("Log2 Counts Per Million, Pval=", pval, "R=", R_val, "FC=", FC, "logCPM=", logCPM))
    chart.Correlation(MYCN_vs_gene2[,2:3], pch = 21, bg=as.factor(MYCN_vs_gene2$MNA), method="pearson", histogram=T)
    #title(paste("Violin Plot:", HUGO_gene_name, "Log2 Counts Per Million \n", "Pval=", pval, "     R=", R_val, "\n          FC=", FC, "      logCPM=", logCPM))
    quartz()
    vioplot(gene_MNA[,1], gene_nonMNA[,1], names = c("MNA", "nonMNA"), col="gray")

    title(paste("Violin Plot:", HUGO_gene_name, "Log2 Counts Per Million \n", "Pval=", pval, "     R=", R_val, "\n          FC=", FC, "      logCPM=", logCPM))
}



read.PWM <- function (file, sep = "\t", header = F, skip = 1, letter.orientation = c("A", "C", "G", "T")) {
    PWM.temp <- read.delim(file, sep = "\t", header = header, skip = skip)
    PWM <- (t(PWM.temp))
    row.names(PWM) <- letter.orientation
    return(PWM)
}

glm_to_chartCor <- function(glm_dataframe, number_of_genes) {
    CPM.an2 <- log2((CPM.an[,4:ncol(CPM.an)] + 1))
    row.names(CPM.an2) <- make.names(CPM.an$V2, unique = T)
    gene_list <- glm_dataframe$external_gene_name[1:number_of_genes]
    MYCN_vs_genes <- data.frame(t(subset(CPM.an2, row.names(CPM.an2) %in% gene_list | row.names(CPM.an2) == "MYCN")))[1:ncol(CPM.an2),]
    MYCN_vs_genes2 <- merge(MYCN_vs_genes, design[1:2], by.x="row.names", by.y="Sample", all.x=T)
    chart.Correlation(MYCN_vs_genes, bg=as.factor(design$MNA), pch = 21, method="pearson", histogram=T)
}

plotMD_FCline <- function(object, column = 1, xlab = "Average log-expression",
ylab = "Expression log-ratio (this sample vs others)",
main = colnames(object)[column], status=NULL)
{
    plotMD(object, column, xlab, ylab, main = colnames(object)[column], status=NULL)
    abline(h=0, col="red", lty=2, lwd=2)
}

plotMD_multiple <- function(object, startcol, endcol, dim)
{
    plot.new()
    samples_to_plot <- endcol - startcol
    par(mfrow=dim)

    for (i in startcol:endcol) {
        plotMD(object, column = i, main = colnames(object)[i], startcol, endcol)
        abline(h=0, col="red", lty=2, lwd=2)
        mtext(side = 2, text = "Expression log-ratio", line = 3)
        mtext(side = 2, text = " (this sample vs others)", line = 2)

        mtext(side = 1, text = "Average", line = 2)
        mtext(side = 1, text = "log-expression", line = 3)
        mtext(side = 1, text = colnames(object)[i], line = 4)
    }
    par(mfrow=c(1,1))
}

matchPWM_in_df <- function(PWM, df) {
    t <- NULL
    for (i in 1:nrow(df)) {
        #t <- print(matchPWM(EGR.t, as.character(TP73AS1.up.seq[1,
        t[i] <- toString(matchPWM(PWM, as.character(df[i,1])))
    }
    t2 <- data.frame(t)
    t3 <- cbind(t2, df$X2)
    t4 <- t3[nchar(as.character(t3$t)) > 0,]
    return(t4)
}

order.cols.by.row <- function(df, rowname)
{
    df[, order(-df[which(rownames(df) == rowname), ]) ]
}

# input requires format as data frame where row.names are ENSG identifiers followed by counts
# cols: ENSG - sample1 - sample2
genewise_corr_plot <- function(df, target_gene, gene_sorted, design_input = "", cex.axis = 1, target_gene_number = 1, gene_sorted_number = 1, height_inch = 8.27, width_inch = 11.69, plot.labels = FALSE,
distance_between_axis_and_labels = 0, square_spacing = 1, square_size = 2, square_text_cex = 1, square_type = 15, text_offset = 0, type = 0, file = 0, dpi = 300, biomart = mart) {

    # df <- log2(CPM.an[4:163] + 1)
    # gene_sorted = "MYCN"
    # target_gene = "TP73-AS1"
    # cex.axis = 1; target_gene_number = 1; gene_sorted_number = 1
    # distance_between_axis_and_labels = 1; square_spacing = 1; square_size = 1; square_text_cex = 1; square_type = 15; text_offset = 0

    # df <- MEG3.bedtools.counts.fc_junc.SRA.txt2.cpm

    # df <- TARGET.log.cpm
    df.temp <- df
    # head(df.temp)
    # gene_sorted = "MYCN"
    # gene_sorted = "MEG3_6_gap"
    # target_gene = "MEG3_7_gap2"; gene_sorted_number = 1; target_gene_number = 1; biomart = mart
    # distance_between_axis_and_labels = 1; square_spacing = 1; square_size = 2; square_text_cex = 1; square_type = 15; text_offset = 0; type = 0; file = 0; dpi = 300; biomart = mart
    # rm(gene_sorted); rm(gene_sorted_number); rm(biomart); rm(distance_between_axis_and_labels) ; rm(square_spacing); rm(square_size); rm(square_text_cex); rm(square_type); rm(text_offset); rm(type); rm(file); rm(dpi); rm(biomart)

    if (is.na(ext_name2ensg(gene_sorted, biomart)[gene_sorted_number,2])) {
        ensg.id <- gene_sorted
    } else {
        ensg.id <- ext_name2ensg(gene_sorted, biomart)[gene_sorted_number,2]
    }

    try(if(grep("ENSG", ensg.id) == 1) print(ext_name2ensg(gene_sorted, biomart)))

    df.sorted <- df.temp[, order(-df.temp[which(rownames(df.temp) == ensg.id), ]) ]

    df1 <- subset(df.sorted, row.names(df.sorted) == ensg.id)
    # df1 <- df1[, order(df1)]
    # df1 <- data.frame(df1)

    df1 <- data.frame(t(rev(df1)))
    df1$V2 <- 1:nrow(df1)
    colnames(df1) <- c(gene_sorted, "V2")
    df1 <- df1[, c("V2", gene_sorted)]

    # shouold now have row.names with samples

    #use this for coloring MNA samples
    if (design_input != "") {
        design_input_matched <- design_input[match(row.names(df1),row.names(design_input)), ]
    }

    try(ensg.target <- ext_name2ensg(target_gene, biomart)[target_gene_number,2])
    if(is.na(ensg.target) == T) ensg.target <- target_gene
    try(print(ext_name2ensg(target_gene, biomart)))

    df2 <- subset(df.sorted, row.names(df.sorted) == ensg.target)
    if (nrow(df2) == 0) print("Your target gene is not present in your dataset")
    df2 <- data.frame(t(rev(df2)))
    df2$V2 <- 1:nrow(df2)
    colnames(df2) <- c(target_gene, "V2")
    df2 <- df2[, c("V2", target_gene)]

    df3 <- cbind(df1[,2], df2[,2]); colnames(df3) <- c(gene_sorted, target_gene)
    ttest <- t.test(df3)
    Pval <- signif(ttest$p.value, 3)
    Rval <- signif(cor(df3)[1,2], 3)

    # start quartz with A4 landscape
    if (type == "pdf") {
        quartz.save(file, type, dpi = dpi)
    } else
    quartz(width = width_inch, height = height_inch)

    par(pin=c(height_inch, width_inch))
    # par(mar=c(5, 4, 4, 5) + 0.1)
    anno_ratio = 4 / 10

    par(mai = c(height_inch * anno_ratio, 1, 0.8, 1)) # margin size specified in in inches

    ## PLOT
    if (design_input != "") {
        plot(df1, col="black", pch=square_type, xlab="", xaxt='n', xpd = NA,
        xlim = c(-(max(nchar(colnames(design_input_matched))) / 2), nrow(df1)), xpd = F)
    } else {
        plot(df1, col="black", pch=square_type, xlab="", xaxt='n', xpd = NA)
    }
    # ylim = c(min(df1[,2]) - ((max(nchar(row.names(df1))) + ncol(design_input_matched)) * 0.3), max(df1[,2]))


    if (plot.labels == TRUE) {
        axis(1, at=1:nrow(df1), pos = min(df1[,2]) - 0.5, labels=row.names(df1), cex.axis = cex.axis, las = 2)
    }
    # axis(2, at=round(min(df1[,2])):round(max(df1[,2])), las = 2)

    ## Allow a second plot on the same graph and add extra space to right margin of plot within frame
    par(new=TRUE) ; par(mai = c(height_inch * anno_ratio, 1, 0.8, 1));
    if  (design_input != "") {
        plot(df2, col="red", pch=square_type, xlab="", ylab="", xaxt = 'n', yaxt = 'n', xpd = NA,
        xlim = c(-(max(nchar(colnames(design_input_matched))) / 2), nrow(df2)))

        # add axis, axis texts and titles
        axis(4, round(min(df2[,2])):round(max(df2[,2])), col="red",col.axis="red", las=1)
        mtext(target_gene,side=4,col="red", line=3)

    } else {
        plot(df2, col="red", pch=square_type, xlab="", ylab="", xaxt = 'n', yaxt = 'n')
        axis(4, round(min(df2[,2])):round(max(df2[,2])), col="red",col.axis="red", las=1)
        mtext(target_gene,side=4,col="red", line=3)
    }
    # ylim = c(min(df2[,2]) - ((max(nchar(row.names(df1))) + ncol(design_input_matched)) * 0.3), max(df2[,2]))
    # ylim = c(min(df2[,2]) - ((max(nchar(row.names(df1))) + ncol(design_input_matched)) * 0.3), max(df2[,2]))


    grid(nx = (nrow(df1)/10), NULL, col = "lightgray", lty = "dotted")
    title(paste(gene_sorted, " vs ", target_gene, "\nR=", Rval))

    if  (design_input != "") {
        print("Possible values for design matrix is: 'YES', 'NO', 'AMPLIFIED', 'NON-AMPLIFIED', MALE', 'FEMALE', 'LOW', 'INTERMEDIATE', 'HIGH', '1', '2', '3', '4', '4S', 'UNKOWN', 'OTHER', 'NA' or empty values")
        ## Add color labeling of clinical characteristics

        # figure out how labelling of characteristics should be done
        #(1) lower y coordinate of figure region
        par()$usr[3] - (1 / par()$plt[4])
        #(2) height of figure region
        par()$usr[3] -(par()$usr[3] - (1 / par()$plt[3]))

        # use (1) and (2) to determine where the first label should be plotted, as well as its distance
        # square_pos = par()$usr[3] - (par()$usr[3] / (par()$usr[3] -(par()$usr[3] - (1 / par()$plt[3]))) * 1/ncol(design_input_matched)) * distance_between_axis_and_labels
        # square_pos = (par()$usr[3] * 7 / 10) * 1/distance_between_axis_and_labels
        # print(par()$usr[3])
        square_pos = par()$mar[3] + distance_between_axis_and_labels
        #- strwidth(row.names(design_input_matched[nchar(design_input_matched) == max(nchar(design_input_matched)), ])[1], cex = cex.axis, units = "inches") * 1.5) * 1/distance_between_axis_and_labels
        # square_distance = par()$usr[3] -(par()$usr[3] - (1 / par()$plt[3] * 1/ncol(design_input_matched)) * square_spacing
        # square_pos = (height_inch * anno_ratio) - (strwidth(row.names(design_input_matched[nchar(design_input_matched) == max(nchar(design_input_matched)), ])[1], cex = cex.axis, units = "inches") * 1.5) * 1/distance_between_axis_and_labels
        ## square_pos = (height_inch * anno_ratio) * 0.5 * 1/distance_between_axis_and_labels
        # square_distance = ((par()$usr[3] -(par()$usr[3] - (1 / par()$plt[3])))  * 1/ncol(design_input_matched)) * square_spacing
        square_distance = height_inch * anno_ratio * 1/ncol(design_input_matched) * square_spacing

        for (i in 0:(ncol(design_input_matched)-1)) {
            for (x in 1:nrow(design_input_matched)) {
                if (design_input_matched[x, (i+1)] == "FEMALE") {
                    color = "red"
                } else if (design_input_matched[x, (i+1)] == "MALE") {
                    color = "blue"

                } else if (design_input_matched[x, (i+1)] == "YES") {
                    color = "green"
                } else if (design_input_matched[x, (i+1)] == "NO") {
                    color = "red"

                } else if (design_input_matched[x, (i+1)] == 'AMPLIFIED') {
                    color = "red"
                } else if (design_input_matched[x, (i+1)] == 'NON-AMPLIFIED') {
                    color = "green"


                } else if (design_input_matched[x, (i+1)] == "LOW") {
                    color = "green"
                } else if (design_input_matched[x, (i+1)] == "INTERMEDIATE") {
                    color = "gold"
                } else if (design_input_matched[x, (i+1)] == "HIGH") {
                    color = "red"

                } else if (design_input_matched[x, (i+1)] == "1") {
                    color = "green"
                } else if (design_input_matched[x, (i+1)] == "2") {
                    color = "darkgreen"
                } else if (design_input_matched[x, (i+1)] == "3") {
                    color = "gold"
                } else if (design_input_matched[x, (i+1)] == "4") {
                    color = "red"
                } else if (design_input_matched[x, (i+1)] == "4S") {
                    color = "blue"

                } else if (design_input_matched[x, (i+1)] == "UNKOWN") {
                    color = "grey"
                } else if (design_input_matched[x, (i+1)] == "OTHER") {
                    color = "grey"
                } else if (design_input_matched[x, (i+1)] == "NA") {
                    color = "grey"

                } else
                color = "black"

                points(x, square_pos - (i * square_distance), pch = square_type, col = color, cex = square_size, xpd = NA)
            }
            # text(0 - nchar(colnames(design_input_matched)[i]), square_pos - (i * square_spacing), colnames(design_input_matched)[i], cex = square_text_size)
            text(0 + text_offset, square_pos - (i * square_distance), colnames(design_input_matched)[i+1], cex = square_text_cex, xpd = NA)
        }}
}


genewise_corr_plot_without_design <- function(df, target_gene, gene_sorted, method = c("spearman", "pearson", "kendall")[1], cex.axis = 1, target_gene_number = 1, gene_sorted_number = 1, height_inch = 8.27, width_inch = 11.69,
type = 0, file = 0, dpi = 300, biomart = mart, square_type = 15, show_labels = F, Rval, col = c("#68382C", "#00A4E6"), use.grid = F, title.with.name = F, cex = 2, cex.y.axis=1.5, gene_sorted_name = "",
y.axis.dist.label = 2.2) {
    # df <- log2(CPM.an[4:163] + 1)
    # gene_sorted = "MYCN"
    # target_gene = "TP73-AS1"
    # cex.axis = 1; target_gene_number = 1; gene_sorted_number = 1
    # distance_between_axis_and_labels = 1; square_spacing = 1; square_size = 1; square_text_cex = 1; square_type = 15; text_offset = 0

    df.temp <- df

    if (is.na(ext_name2ensg(gene_sorted, biomart)[gene_sorted_number,2])) {
        ensg.id <- gene_sorted
    } else {
        ensg.id <- ext_name2ensg(gene_sorted, biomart)[gene_sorted_number,2]
    }

    try(if(grep("ENSG", ensg.id) == 1) print(ext_name2ensg(gene_sorted, biomart)))

    df.sorted <- df.temp[, order(-df.temp[which(rownames(df.temp) == ensg.id), ]) ]

    df1 <- subset(df.sorted, row.names(df.sorted) == ensg.id)
    df1 <- data.frame(t(rev(df1)))
    df1$V2 <- 1:nrow(df1)
    colnames(df1) <- c(gene_sorted, "V2")
    df1 <- df1[, c("V2", gene_sorted)]

    if (is.na(ext_name2ensg(target_gene, biomart)[target_gene_number,2])) {
        ensg.target <- target_gene
    } else {
        ensg.target <- ext_name2ensg(target_gene, biomart)[target_gene_number,2]
        try(print(ext_name2ensg(target_gene, biomart)))
    }
    print(ensg.target)
    # try(ensg.target <- ext_name2ensg(target_gene, biomart)[target_gene_number,2])
    # if(is.na(ensg.target) == T) ensg.target <- target_gene
    # try(print(ext_name2ensg(target_gene, biomart)))

    df2 <- subset(df.sorted, row.names(df.sorted) == ensg.target)
    if (nrow(df2) == 0) print("Your target gene is not present in your dataset")
    df2 <- data.frame(t(rev(df2)))
    df2$V2 <- 1:nrow(df2)
    colnames(df2) <- c(target_gene, "V2")
    df2 <- df2[, c("V2", target_gene)]

    df3 <- cbind(df1[,2], df2[,2]); colnames(df3) <- c(gene_sorted, target_gene)
    ttest <- t.test(df3)
    Pval <- signif(ttest$p.value, 3)

    if (Rval == "") {
    Rval <- signif(cor(df3, method = method)[1,2], 3)
    }

    # start quartz with A4 landscape
    if (type == "pdf") {
        quartz.save(file, type, dpi = dpi)
    } else if  (type == "quartz") {
        quartz(width = width_inch, height = height_inch)
    }

    #par(pin=c(height_inch, width_inch))
    # par(mar=c(5, 4, 4, 5) + 0.1)

    #par(mai = c(height_inch, 1, 0.8, 1)) # margin size specified in in inches

    ## PLOT
    par(mar=c(5, 4, 4, 5) + 0.1)
    # quartz()
    plot(df1, col=col[1], pch=square_type, xlab="", ylab="", xaxt='n', yaxt='n', xpd = NA)

    axis(2, col="black",col.axis=col[1], las=1, cex.axis=cex.y.axis)

    if (gene_sorted_name == "") {
    mtext(gene_sorted,side=2,col=col[1], line=y.axis.dist.label, cex = cex)
    } else {
    mtext(gene_sorted_name,side=2,col=col[1], line=y.axis.dist.label, cex = cex)
    }

    if  (show_labels == T) axis(1, at=1:nrow(df1), pos = min(df1[,2]) - 0.5, labels=row.names(df1), cex.axis = cex.axis, las = 2)
    # axis(2, at=round(min(df1[,2])):round(max(df1[,2])), las = 2)

    ## Allow a second plot on the same graph and add extra space to right margin of plot within frame
    #par(mai = c(height_inch, 1, 0.8, 1));
    par(mar=c(5, 4, 4, 5) + 0.1)

    par(new=TRUE) ;  plot(
    df2, col=col[2], pch=square_type, xlab="", ylab="", xaxt = 'n', yaxt = 'n', xpd = NA)

    # add axis, axis texts and titles
    axis(4, col="black",col.axis=col[2], las=1, cex.axis=cex.y.axis)

    mtext(target_gene,side=4,col=col[2], line=y.axis.dist.label, cex = cex)

    if (use.grid == T) grid(nx = (nrow(df1)/10), NULL, col = "lightgray", lty = "dotted")
    if (title.with.name == T) {
        title(paste(gene_sorted, " vs ", target_gene, "\nR=", Rval))
        } else {
        title(paste("R =", Rval), line = 0.5, cex = cex)
        }
}

ggenewise_corr_plot <- function(df, target_gene, gene_sorted, method = c("spearman", "pearson", "kendall")[1], cex.axis = 1, target_gene_number = 1, gene_sorted_number = 1, height_inch = 8.27, width_inch = 11.69,
                                              type = 0, file = 0, dpi = 300, biomart = mart, square_type = 15, show_labels = F) {
    # df <- log2(CPM.an[4:163] + 1)
    # gene_sorted = "MYCN"
    # target_gene = "LINC00657"
    # cex.axis = 1; target_gene_number = 1; gene_sorted_number = 1; biomart = mart; method = "pearson"
    # distance_between_axis_and_labels = 1; square_spacing = 1; square_size = 1; square_text_cex = 1; square_type = 15; text_offset = 0

    df.temp <- df

    if (is.na(ext_name2ensg(gene_sorted, biomart)[gene_sorted_number,2])) {
        ensg.id <- gene_sorted
    } else {
        ensg.id <- ext_name2ensg(gene_sorted, biomart)[gene_sorted_numberr,2]
    }

    try(if(grep("ENSG", ensg.id) == 1) print(ext_name2ensg(gene_sorted, biomart)))

    df.sorted <- df.temp[, order(-df.temp[which(rownames(df.temp) == ensg.id), ]) ]

    df1 <- subset(df.sorted, row.names(df.sorted) == ensg.id)
    df1 <- data.frame(t(rev(df1)))
    df1$V2 <- 1:nrow(df1)
    colnames(df1) <- c(gene_sorted, "V2")
    df1 <- df1[, c("V2", gene_sorted)]

    if (is.na(ext_name2ensg(target_gene, biomart)[target_gene_number,2])) {
        ensg.target <- target_gene
    } else {
        ensg.target <- ext_name2ensg(target_gene, biomart)[target_gene_number,2]
        try(print(ext_name2ensg(target_gene, biomart)))
    }
    print(ensg.target)
    # try(ensg.target <- ext_name2ensg(target_gene, biomart)[target_gene_number,2])
    # if(is.na(ensg.target) == T) ensg.target <- target_gene
    # try(print(ext_name2ensg(target_gene, biomart)))

    df2 <- subset(df.sorted, row.names(df.sorted) == ensg.target)
    if (nrow(df2) == 0) print("Your target gene is not present in your dataset")
    df2 <- data.frame(t(rev(df2)))
    df2$V2 <- 1:nrow(df2)
    colnames(df2) <- c(target_gene, "V2")
    df2 <- df2[, c("V2", target_gene)]

    df3 <- cbind(df1[,2], df2[,2]); colnames(df3) <- c(gene_sorted, target_gene)
    ttest <- t.test(df3)
    Pval <- signif(ttest$p.value, 3)
    Rval <- signif(cor(df3, method = method)[1,2], 3)

    # start quartz with A4 landscape
    if (type == "pdf") {
        quartz.save(file, type, dpi = dpi)
    } else if  (type == "quartz") {
        quartz(width = width_inch, height = height_inch)
    }

    #par(pin=c(height_inch, width_inch))
    # par(mar=c(5, 4, 4, 5) + 0.1)

    #par(mai = c(height_inch, 1, 0.8, 1)) # margin size specified in in inches

    ## PLOT
    colnames(df1) <- c("pos", "expression")
    colnames(df2) <- c("pos", "expression")

    grid.newpage()

    p1 <- ggplot(df1, aes(pos, expression)) + geom_point()
    # paste("R =", round(Rval,2)
    p2 <- ggplot(df2, aes(pos, expression)) + geom_point(color = "red") + labs(x = "", y = gene_sorted) + guides(colour = guide_legend(NULL), shape = guide_legend(NULL)) + theme_base() %+replace%
        theme(panel.background = element_rect(fill = NA))

    {
    # extract gtable
    g1 <- ggplot_gtable(ggplot_build(p1))
    g2 <- ggplot_gtable(ggplot_build(p2))

    # overlap the panel of 2nd plot on that of 1st plot
    pp <- c(subset(g1$layout, name == "panel", se = t:r))
    g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t,
                         pp$l, pp$b, pp$l)

    # axis tweaks
    ia <- which(g2$layout$name == "axis-l")
    ga <- g2$grobs[[ia]]
    ax <- ga$children[[2]]
    ax$widths <- rev(ax$widths)
    ax$grobs <- rev(ax$grobs)
    ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
    g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
    g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)

    # draw it
    grid.draw(g)
        }


    if  (show_labels == T) axis(1, at=1:nrow(df1), pos = min(df1[,2]) - 0.5, labels=row.names(df1), cex.axis = cex.axis, las = 2)
    # axis(2, at=round(min(df1[,2])):round(max(df1[,2])), las = 2)

    ## Allow a second plot on the same graph and add extra space to right margin of plot within frame
    #par(mai = c(height_inch, 1, 0.8, 1));
    par(mar=c(5, 4, 4, 5) + 0.1)

    par(new=TRUE) ;  plot(
        df2, col="red", pch=square_type, xlab="", ylab="", xaxt = 'n', yaxt = 'n', xpd = NA)

    round(min(df2[,2])):round(max(df2[,2]))

    rep(round(min(df2[,2])) + 1, 2)

    # add axis, axis texts and titles
    axis(4, col="red",col.axis="red", las=1)

    mtext(target_gene,side=4,col="red", line=3)
    grid(nx = (nrow(df1)/10), NULL, col = "lightgray", lty = "dotted")
    title(paste(gene_sorted, " vs ", target_gene, "\nR=", Rval))
}

genewise_corr_Rval <- function(df, gene_sorted, gene_list, target_gene_number = 1, gene_sorted_number = 1, biomart = mart) {
    # df = TARGET.log.cpm; gene_sorted = "MYCN"; gene_list = list(t3); target_gene_number = 1; gene_sorted_number = 1; biomart = mart

    df.temp <- df

    ensg.id <- ext_name2ensg(gene_sorted, biomart)[gene_sorted_number,2]
    print(ext_name2ensg(gene_sorted, biomart))

    df.sorted <- df.temp[, order(-df.temp[which(rownames(df.temp) == ensg.id), ]) ]

    df1 <- subset(df.sorted, row.names(df.sorted) == ensg.id)
    df1 <- data.frame(t(rev(df1)))
    df1$V2 <- 1:nrow(df1)
    colnames(df1) <- c(gene_sorted, "V2")
    df1 <- df1[, c("V2", gene_sorted)]

    corr.res <- data.frame(ensembl_gene_id = character(0), Rval = numeric(0), stringsAsFactors = FALSE)

    for (i in 1:length(gene_list)) {

        target_gene = gene_list[i]

        df2 <- subset(df.sorted, row.names(df.sorted) == target_gene)
        if (nrow(df2) == 0) print("Your target gene is not present in your dataset")
        df2 <- data.frame(t(rev(df2)))
        df2$V2 <- 1:nrow(df2)
        colnames(df2) <- c(target_gene, "V2")
        df2 <- df2[, c("V2", target_gene)]

        df3 <- cbind(df1[,2], df2[,2]); colnames(df3) <- c(gene_sorted, target_gene)
        ttest <- t.test(df3)
        # Pval <- signif(ttest$p.value, 3)
        Rval <- signif(cor(df3)[1,2], 3)

        corr.res[i, ] <- list(target_gene, Rval)
    }

    corr.res.an <- ensg2ext_name_biotype(corr.res$ensembl_gene_id, biomart)
    corr.res.an.R <- merge(corr.res.an, corr.res, by = "ensembl_gene_id", all = T)

    corr.res.an.R[order(corr.res.an.R$Rval) , ]

}

ddct <- function(file, reference_genes, controls, treated, external_data = "") {

    comparisons <- data.frame(
        Sample1 = controls,
        Sample2 = treated
    )

    cq_values <- read.csv(file, skip = grep("Well", readLines(file)) - 1)
    cq_values <- cq_values[c("Sample.Name", "Detector", "Ct")]
    cq_values$Ct <- as.numeric(gsub("Undetermined", 40, cq_values$Ct))

    print(cq_values$Sample.Name[!duplicated(cq_values$Sample.Name)])
    print(cq_values$Detector[!duplicated(cq_values$Detector)])

    ## Calculate average CT for reference genes
    samples <- cq_values[!duplicated(cq_values$Sample.Name), ]
    samples <- samples[grep("NTC", samples$Detector, invert = T), "Sample.Name"]

    t <- cq_values[cq_values$Detector %in% reference_genes, ]

    Ref_CT <- data.frame(Sample.Name = samples)
    for (i in 1:length(samples) ) {
        Ref_CT$Ref_CT[i] <- 2^mean(log2(t[t$Sample.Name == samples[i], "Ct"]))
        Ref_CT$sd[i] <- sd(t[t$Sample.Name == samples[i], "Ct"])
    }

    cq_values <- cq_values[grep("NTC", cq_values$Detector, invert = T), ]
    cq_values_ref <- merge(cq_values, Ref_CT, by = "Sample.Name", all = T)

    cq_values_ref$dCT <- round(2^(cq_values_ref$Ref_CT - cq_values_ref$Ct), 6)

    genes <- levels(cq_values_ref$Detector)
    genes <- genes[!genes %in% reference_genes]
    genes <- genes[grep("NTC", genes, invert = T)]

    replicates <- nrow(cq_values_ref[cq_values_ref$Detector == genes[i] & cq_values_ref$Sample.Name == samples[i], ])
    ddct_res <- data.frame(Detector = genes)

    for (k in 1:length(genes)) {
        for (j in 1:nrow(comparisons)) {
            for (z in 1:replicates) {
                ddct_res[k, paste(comparisons[j, 2], "_vs_", comparisons[j, 1], "_rep", z, sep = "")] <- cq_values_ref[cq_values_ref$Detector == genes[k] & cq_values_ref$Sample.Name == as.character(comparisons[j, 2]), "dCT"][z] / cq_values_ref[cq_values_ref$Detector == genes[k] & cq_values_ref$Sample.Name == as.character(comparisons[j, 1]), "dCT"][z]
            }
        }
    }

    ddct_res[2:ncol(ddct_res)] <- apply(ddct_res[2:ncol(ddct_res)], c(1,2), function(x) ifelse(x < 1, -1/x, x))
    row.names(ddct_res) <- ddct_res$Detector
    ddct_res$Detector <- NULL

    t <- t(ddct_res)

    row.names(t) <- gsub("_rep1", "", row.names(t))
    row.names(t) <- gsub("_rep2", "", row.names(t))
    row.names(t) <- gsub("_rep3", "", row.names(t))
    row.names(t) <- gsub("_rep4", "", row.names(t))

    library(data.table)

    t <- data.table::melt(t)

    ggplot(t, aes(X2, value)) + geom_boxplot() + geom_point(aes(colour = factor(t$X1), shape = factor(t$X1)), size = 4) +
        labs(x = "", y = "Fold change") + guides(colour = guide_legend(NULL), shape = guide_legend(NULL)) + theme_base()

    ggsave("ddCT_plot.pdf", device = "pdf", path = dirname(file), dpi = 300)

    if(external_data != "") {
        print("External data must be a data frame with two columns: external_gene_name and FoldChange")
        external_data <- seq16_P_cpm.df2[c("external_gene_name", "TP73AS1_log2FC")]
        colnames(external_data) <- c("external_gene_name", "foldchange")
        external_data <- external_data[external_data$external_gene_name %in% t$X2, ]

        t2 <- data.frame(X1 = as.factor(c(as.character(t$X1), rep("external_data", nrow(external_data)))),
                         X2 = c(as.character(t$X2), external_data$external_gene_name),
                         value = c(t$value, external_data$foldchange))

        z <- merge(external_data, t, by.x = "external_gene_name", by.y = "X2")

        ggplot(t2, aes(X2, value)) + geom_boxplot(data = t, aes(X2, value)) + geom_point(aes(colour = factor(t2$X1), shape = factor(t2$X1)), size = 4) +
            labs(x = paste("Rs = ", round(min(cor(z[c(2,4)], method = "spearman")), 2)), y = "Fold change") +
            guides(colour = guide_legend(NULL), shape = guide_legend(NULL)) + theme_base()
        ggsave("ddCT_plot_external_data.pdf", device = "pdf", path = dirname(file), dpi = 300)
    }

    # Produce a table of means and sd
    ddct_means <- cast(t, X1~X2, mean)
    ddct_sd <- cast(t, X1~X2, sd)
    ddct_means_melt <- melt(ddct_means)
    row.names(ddct_means_melt) <- paste(ddct_means_melt$X1, ddct_means_melt$X2, sep = "_")
    colnames(ddct_means_melt) <- c("comparison", "mean", "Detector")
    ddct_sd_melt <- melt(ddct_sd)
    row.names(ddct_sd_melt) <- paste(ddct_sd_melt$X1, ddct_sd_melt$X2, sep = "_")
    colnames(ddct_sd_melt) <- gsub("value", "sd", colnames(ddct_sd_melt))
    ddct_means_melt <- merge(ddct_means_melt, ddct_sd_melt, by = "row.names")
    row.names(ddct_means_melt) <- ddct_means_melt$Row.names
    ddct_means_melt <- ddct_means_melt[c("comparison", "Detector", "mean", "sd")]
    ddct_means_melt <- ddct_means_melt[order(ddct_means_melt$Detector), ]

    write.table(ddct_means_melt, file = paste(dirname(file), "ddct_means.txt", sep ="/"), sep = "\t", quote = F)
    write.table(ddct_res, file = paste(dirname(file), "ddct_res.txt", sep ="/"), sep = "\t", quote = F)
    write.table(cq_values_ref, file = paste(dirname(file), "cq_values_dct.txt", sep ="/"), sep = "\t", quote = F)
}

# input format is a DGEList (row.names corresponding to ENSG identifiers) and a gene name with HUGO identifier
group_by_expression <- function(df, gene_name, number_of_groups = 2, cex.axis = 0.5, square_type = 15, biomart = mart) {

    # df = SEQC
    # gene_name = "MYCN"
    # number_of_groups = 2; label.cex = 0.3; cex.axis = 0.5; distance_between_axis_and_labels = 0; square_spacing = 0.2; square_size = 2; square_type = 15; biomart = mart

    gene <- ifelse (is.na(hugo2ensg(gene_name, biomart)[1,2]) == FALSE, hugo2ensg(gene_name, biomart)[1,2], gene_name)

    # DGEList.temp <- cpm(DGEList, log=T)
    df.temp <- subset(df, row.names(df) == gene)
    df.temp.sorted <- t(data.frame(df.temp[, order(df.temp[which(rownames(df.temp) == gene), ]) ])); colnames(df.temp.sorted) <- "CPM"

    # initalize variables for the for loop
    # initial_group = round(0.05 * nrow(DGEList.gene.sorted))

    # p.values <- data.frame(c(rep(1, initial_group - 1))); colnames(p.values) <- "pval"

    # for (i in initial_group:nrow(DGEList.gene.sorted)) {
    #    wilcox.rank <- wilcox.test(DGEList.gene.sorted[1:i, ], DGEList.gene.sorted[(i+1):nrow(DGEList.gene.sorted), ])
    #   p.values[nrow(p.values)+1, ] <- wilcox.rank$p.value
    # }

    cluster.df <- data.frame(df.temp.sorted)
    cluster.df$samplenumber <- 1:nrow(cluster.df)

    # assign group by Gaussian distribution
    library(mclust)
    fit <- Mclust(cluster.df$CPM, G=number_of_groups)
    cluster.fit <- cluster.df
    cluster.fit$group <- fit$classification

    # row.temp <- row.names(p.values[p.values$pval == min(p.values[,1]), ])
    # p.values$group <- c(rep(1, row.temp[length(row.temp)]),
    #  rep(0 , nrow(p.values) - as.numeric(row.temp[length(row.temp)])))
    # plot(p.values[4:5], col="black", pch=ifelse(p.values$group==1, 19, 21))

    design_input_matched <- design_input[match(cluster.fit$samples,row.names(design_input)), ]

    plot(cluster.fit$samplenumber, cluster.fit$CPM, col="black", xlab = "", ylab = "CPM (log2)", yaxt = 'n', xaxt='n',
    pch = ifelse(cluster.fit$group > 2, 1, ifelse(cluster.fit$group == 1, 4, 19)))
    title(gene_name)

    # mtext(paste(round(countMatches(1, cluster.fit$group) / nrow(cluster.fit) * 100, 1), "%", round(countMatches(2, cluster.fit$group) / nrow(cluster.fit) * 100, 1), "%", round(countMatches(3, cluster.fit$group) / nrow(cluster.fit) * 100, 1), "%"), side = 3, cex = 0.4)
    # par(las=2)
    # axis(1, at=1:nrow(cluster.fit), pos = min(cluster.fit$CPM) - 0.5, labels=cluster.fit$samples, cex.axis = cex.axis)
    # axis(2, at=round(min(cluster.fit$CPM)):round(max(cluster.fit$CPM)), cex.axis = 0.5)

    print(paste("Number of clusters:", max(fit$classification)))
    return(cluster.fit)
}


group_by_expression_with_design_matrix <- function(df, gene_name, number_of_groups = 2, design_input, label.cex = 0.3, cex.axis = 0.5, distance_between_axis_and_labels = 0, square_spacing = 0.2, square_size = 2, square_type = 15, biomart = mart){

    df = SEQC
    gene_name = "MYCN"
    number_of_groups = 2; label.cex = 0.3; cex.axis = 0.5; distance_between_axis_and_labels = 0; square_spacing = 0.2; square_size = 2; square_type = 15; biomart = mart

    gene <- ifelse (is.na(hugo2ensg(gene_name, biomart)[1,2]) == FALSE, hugo2ensg(gene_name, biomart)[1,2], gene_name)

    # DGEList.temp <- cpm(DGEList, log=T)
    df.temp <- subset(df, row.names(df) == gene)
    df.temp.sorted <- t(data.frame(df.temp[, order(df.temp[which(rownames(df.temp) == gene), ]) ])); colnames(df.temp.sorted) <- "CPM"

    # initalize variables for the for loop
    # initial_group = round(0.05 * nrow(DGEList.gene.sorted))

    # p.values <- data.frame(c(rep(1, initial_group - 1))); colnames(p.values) <- "pval"

    # for (i in initial_group:nrow(DGEList.gene.sorted)) {
    #    wilcox.rank <- wilcox.test(DGEList.gene.sorted[1:i, ], DGEList.gene.sorted[(i+1):nrow(DGEList.gene.sorted), ])
    #   p.values[nrow(p.values)+1, ] <- wilcox.rank$p.value
    # }

    cluster.df <- data.frame(df.temp.sorted)
    cluster.df$samplenumber <- 1:nrow(cluster.df)


    # assign group by Gaussian distribution
    library(mclust)
    fit <- Mclust(cluster.df$CPM, G=number_of_groups)
    cluster.fit <- cluster.df
    cluster.fit$group <- fit$classification

    # row.temp <- row.names(p.values[p.values$pval == min(p.values[,1]), ])
    # p.values$group <- c(rep(1, row.temp[length(row.temp)]),
    #  rep(0 , nrow(p.values) - as.numeric(row.temp[length(row.temp)])))
    # plot(p.values[4:5], col="black", pch=ifelse(p.values$group==1, 19, 21))

    design_input_matched <- design_input[match(cluster.fit$samples,row.names(design_input)), ]

    plot(cluster.fit$samplenumber, cluster.fit$CPM, xlim = c(-(max(nchar(colnames(design_input_matched))) / 2), nrow(cluster.fit)), ylim = c(min(cluster.fit$CPM) - ((max(nchar(cluster.fit$samples)) + ncol(design_input_matched)) * 0.3), max(cluster.fit$CPM)),  col="black", xlab = "", ylab = "CPM (log2)", yaxt = 'n', xaxt='n',
    pch = ifelse(cluster.fit$group > 2, 1, ifelse(cluster.fit$group == 1, 4, 19)))
    title(gene_name)

    # mtext(paste(round(countMatches(1, cluster.fit$group) / nrow(cluster.fit) * 100, 1), "%", round(countMatches(2, cluster.fit$group) / nrow(cluster.fit) * 100, 1), "%", round(countMatches(3, cluster.fit$group) / nrow(cluster.fit) * 100, 1), "%"), side = 3, cex = 0.4)
    par(las=2)
    axis(1, at=1:nrow(cluster.fit), pos = min(cluster.fit$CPM) - 0.5, labels=cluster.fit$samples, cex.axis = cex.axis)
    axis(2, at=round(min(cluster.fit$CPM)):round(max(cluster.fit$CPM)), cex.axis = 0.5)

    print("Possible values for design matrix is: 'YES', 'NO', 'AMPLIFIED', 'NON-AMPLIFIED', MALE', 'FEMALE', 'LOW', 'INTERMEDIATE', 'HIGH', '1', '2', '3', '4', '4S', 'UNKOWN', 'OTHER', 'NA' or empty values")
    ## Add color labeling of clinical characteristics

    # initiliaze coloring
    square_pos = min(cluster.fit$CPM)-(max(nchar(cluster.fit$samples))/2) - distance_between_axis_and_labels

    for (i in 0:(ncol(design_input_matched)-1)) {
        for (x in 1:nrow(design_input_matched)) {
            if (design_input_matched[x, (i+1)] == "FEMALE") {
                color = "red"
            } else if (design_input_matched[x, (i+1)] == "MALE") {
                color = "blue"

            } else if (design_input_matched[x, (i+1)] == "YES") {
                color = "green"
            } else if (design_input_matched[x, (i+1)] == "NO") {
                color = "red"

            } else if (design_input_matched[x, (i+1)] == 'AMPLIFIED') {
                color = "red"
            } else if (design_input_matched[x, (i+1)] == 'NON-AMPLIFIED') {
                color = "green"


            } else if (design_input_matched[x, (i+1)] == "LOW") {
                color = "green"
            } else if (design_input_matched[x, (i+1)] == "INTERMEDIATE") {
                color = "gold"
            } else if (design_input_matched[x, (i+1)] == "HIGH") {
                color = "red"

            } else if (design_input_matched[x, (i+1)] == "1") {
                color = "green"
            } else if (design_input_matched[x, (i+1)] == "2") {
                color = "darkgreen"
            } else if (design_input_matched[x, (i+1)] == "3") {
                color = "gold"
            } else if (design_input_matched[x, (i+1)] == "4") {
                color = "red"
            } else if (design_input_matched[x, (i+1)] == "4S") {
                color = "blue"

            } else if (design_input_matched[x, (i+1)] == "UNKOWN") {
                color = "grey"
            } else if (design_input_matched[x, (i+1)] == "OTHER") {
                color = "grey"
            } else if (design_input_matched[x, (i+1)] == "NA") {
                color = "grey"

            } else
            color = "white"

            # print(color)
            points(x, square_pos - (i * square_spacing), pch = square_type, col = color)
        }
        # text(0 - nchar(colnames(design_input_matched)[i]), square_pos - (i * square_spacing), colnames(design_input_matched)[i], cex = square_text_size)
        text(0 - max(nchar(colnames(design_input_matched))) / 2, square_pos - (i * square_spacing), colnames(design_input_matched)[i+1], cex = label.cex)
    }


    print(paste("Number of clusters:", max(fit$classification)))
    hugo2ensg(gene, biomart)
    res <- merge(design_input_matched, p.values2, by.x="row.names", by.y="samples")
    res$pval <-NULL
    res$samplenumber <- NULL
    return(res)
    # return(p.values2[c(2,5,4)])
}

correlation.matrix <- function(df, gene_one, gene_list, biomart = mart, gene_sorted_number = 1, target_sorted_number = 1) {

    # gene_list <- MEG3$external_gene_name[3:5]
    # gene_one = "MEG3_6_gap"
    # df = TARGET.log.cpm

    if (is.na(ext_name2ensg(gene_one, biomart)[gene_sorted_number,2]) == FALSE) gene_one_tf <- ext_name2ensg(gene_one, biomart)[gene_sorted_number,2]

    # cor.round <- function(df) {
    #    t[i,1] <- data.frame(round(cor(t(df[row.names(df) == gene_one_tf, ]), t(df[row.names(df) == gene_two, ])), 2))
    #    row.names(t)[i] <- paste(gene_two, i, sep = "_")
    #    t$name[i] <- gene_two
    # }

    # t <- NULL

    for (i in 1:length(gene_list)) {
        print(gene_list[i])
        if (is.na(ext_name2ensg(gene_list[i], biomart)[target_sorted_number,2]) == FALSE) {
            gene_two <- ext_name2ensg(gene_list[i], biomart)[gene_sorted_number,2]
        } else {
            gene_two <- gene_list[i]
        }

        if (i == 1) {
            t <- data.frame(round(cor(t(df[row.names(df) == gene_one_tf, ]), t(df[row.names(df) == gene_two, ])), 2))
            row.names(t) <- paste(gene_two, i, sep = "_")
            colnames(t) <- gene_one
            t$name <- gene_two

        } else {
            t[i,1] <- data.frame(round(cor(t(df[row.names(df) == gene_one_tf, ]), t(df[row.names(df) == gene_two, ])), 2))
            row.names(t)[i] <- paste(gene_two, i, sep = "_")
            t$name[i] <- gene_two
            # tryCatch(cor.round(df)) #, error = function() i = i + 1)
            # rbind(cor.round(df), t)
        }
    }

    # colnames(t) <- gene_one

    t <- t[order(t),]
    t <- t[!is.na(t[1]), ]
    return(t)
}

correlation.matrix.ensg <- function(df, gene_one_supplied_as_ensg_id, gene_list_with_ensg_ids, biomart = mart) {

    # function to get the name of an object
    myname <- function(object) {out<-deparse(substitute(object)); out}

    for (i in 1:length(gene_list_with_ensg_ids)) {

        if (i == 1) {
            t <- data.frame(round(cor(t(df[row.names(df) == gene_one_supplied_as_ensg_id, ]), t(df[row.names(df) == gene_list_with_ensg_ids[i], ])), 2))
            row.names(t) <- paste(gene_list_with_ensg_ids[i], i, sep = "_")
            colnames(t) <- paste(myname(df), gene_one_supplied_as_ensg_id, sep = "_")
            t$name <- gene_one_supplied_as_ensg_id[i]

        } else {
            t[i,1] <- data.frame(round(cor(t(df[row.names(df) == gene_one_supplied_as_ensg_id, ]), t(df[row.names(df) == gene_list_with_ensg_ids[i], ])), 2))
            row.names(t)[i] <- paste(gene_list_with_ensg_ids[i], i, sep = "_")
            t$name[i] <- gene_list_with_ensg_ids[i]
            # tryCatch(cor.round(df)) #, error = function() i = i + 1)
            # rbind(cor.round(df), t)
        }
    }

    # colnames(t) <- gene_one

    t <- t[order(t),]
    t <- t[!is.na(t[1]), ]
    return(t)
}

correlation.matrix.ensg.to.file <- function(df, gene_one_supplied_as_ensg_id, gene_list_with_ensg_ids, biomart = mart, rowStart =1, file = "test.txt") {

    # df <- head(TARGET.log.cpm, 10)
    # gene_one_supplied_as_ensg_id <- "ENSG00000260032"
    # gene_one_supplied_as_ensg_id <- "ENSG00000223972"
    # gene_list_with_ensg_ids <- row.names(df)
    # rowStart = 1
    # biomart = mart
    # file = "/Volumes/Peter/MEG3/test.txt"
    # rm(df); rm(gene_one_supplied_as_ensg_id); rm(gene_list_with_ensg_ids); rm(biomart); rm(rowStart)

    gene_one <- t(df[row.names(df) == gene_one_supplied_as_ensg_id, ])

    for (i in rowStart:nrow(df)) {
        if (i == 1) {
            cat(round(cor(gene_one, t(df[row.names(df) == gene_list_with_ensg_ids[i], ])), 2), gene_list_with_ensg_ids[i], sep = "\t", file = file, append = F)
            cat("", sep = "\n", file = file, append = T)
        } else {
            cat(round(cor(gene_one, t(df[row.names(df) == gene_list_with_ensg_ids[i], ])), 2), gene_list_with_ensg_ids[i], sep = "\t", file = file, append = T)
            cat("", sep = "\n", file = file, append = T)
        }}
}

cor.parallell.data.frame <- function(df, gene_one_id, gene_list_ids, file = "test.txt", correlation_type = "pearson", read.file = F, annotate = F, no_cores = "") {

    ## This function iterates one row at a time and then writes out each correlation to a file. Its input requires a data frame (can not be a matrix)
    ## It is done in a manner that requires low memory due to the fact that only one row is is read and calculated at a time, and each single correlation is then appended
    ## to a file (instead of heaping up in memory)
    ## All your cores, except for one will be used -- unless specified manually.

    library(foreach)
    library(doParallel)

    print("Correlation methods may be one of pearson, kendall, spearman")
    print("Input (df) must be a data frame with row.names = gene names, colnames = samples")
    print("No cores used may be specified manually or it will be automatically designated using all available cores - 1")

    gene_one <- t(df[row.names(df) == gene_one_id, ])

    if (no_cores == "") {
        no_cores <- detectCores() - 1
        registerDoParallel(no_cores)
    }

    # Creates and overwrites file if it already exists
    system(paste("mkdir -p", dirname(file)), ignore.stdout = T, ignore.stderr = T)
    file.create(file)

    nr_times = floor(nrow(df) / no_cores) - 1
    remainder = nrow(df) %% no_cores

    ## Clean-up any existing temp files
    for (i in 0:nr_times) {
        system(paste("trash", paste(file, z, "temp.txt", sep = "_")), ignore.stdout = T, ignore.stderr = T)
    }

    for (i in 0:nr_times) {
        foreach (z = 1:no_cores,
                 .combine = rbind)  %dopar%
            cat(
                round(cor(gene_one, # correlate gene of interest to ...
                          t(df[row.names(df) == gene_list_ids[((i * no_cores) + z)], ]), method = correlation_type), 2), #  gene x in the gene list, and round
                "\t", gene_list_ids[((i * no_cores) + z)], "\n",
                sep = "", file = paste(file, z, "temp.txt", sep = "_"), append = T) # add tab and newline, then append to file
    }

    z = length(grep(paste(file, ".", "temp.txt", sep = "_"), list.files(dirname(file)))) + 1

    file.rem = paste(file, z, "temp.txt", sep = "_")

    # add in the rest
    for (i in 1:remainder) {
        cat(
            round(cor(gene_one, # Pearson correlation with gene one ...
                      t(df[row.names(df) == gene_list_ids[((nr_times + 1) * no_cores) + i], ]), method = correlation_type), 2), # compared against gene x in the gene list, and round
            "\t", gene_list_ids[((nr_times + 1) * no_cores) + i], "\n",
            sep = "", file = file.rem, append = T)
    }

    stopImplicitCluster()

    file.cat <- gsub(paste(z, "temp.txt", sep = "_"), "*_temp.txt", file.rem)

    print("Number of correlations in temp files:")
    system(paste("wc -l", file.cat))

    ## combine all files
    system(paste("cat", file.cat, ">", file))
    system(paste("trash", file.cat), ignore.stdout = T, ignore.stderr = T)

    if (read.file == T || annotate == T) {
        assign.name <- gsub(paste(dirname(file), "/", sep = ""), "", file)
        assign(assign.name, read.delim(file, header = F))

        if (annotate == T) {
            t <- get(assign.name)
            colnames(t) <- c("correlation_type", "gene_name")

            t <- ensg2ext_name_biotype(t$gene_name, combine = T)

            t <- t[order(-abs(t$correlation_type)), ]
            colnames(t) <- c("ensembl_gene_id", "external_gene_name", "gene_biotype", correlation_type)
            assign(assign.name, t, envir = .GlobalEnv)
            rm(t)
        } else {
            t <- get(assign.name)
            colnames(t) <- c("correlation_type", "gene_name")
            t <- t[order(-abs(t$correlation_type)), ]
            colnames(t) <- c(correlation_type, "gene_name")
            assign(assign.name, t, envir = .GlobalEnv)
            rm(t)
        }
    }

    print("Number of correlations in final file:")
    system(paste("wc -l", file))
}



waterfall.plot <- function(DGEList, gene_name, number_of_groups = 2, design_input = 1, label.cex = 0.3, cex.axis = 0.5, distance_between_axis_and_labels = 0, square_spacing = 0.2, square_size = 2, square_type = 15){


    # DGEList = Qlnc_MNA; gene_name = "ENSG00000165197"; number_of_groups = 2; design_input = MYCN_high_expressing_samples
    # label.cex = 0.3; cex.axis = 0.5; distance_between_axis_and_labels = 0; square_spacing = 0.2; square_size = 2; square_type = 15
    # rm(DGEList); rm(gene_name); rm(number_of_groups); rm(design_input); rm(label.cex); rm(cex.axis); rm(distance_between_axis_and_labels); rm(square_spacing); rm(square_size); rm(square_pos); rm(square_type)

    gene <- ifelse (is.na(hugo2ensg(gene_name)[1,2]) == FALSE, hugo2ensg(gene_name)[1,2], gene_name)

    DGEList.temp <- cpm(DGEList, log=T)
    DGEList.gene <- subset(DGEList.temp, row.names(DGEList.temp) == gene)
    DGEList.gene.sorted <- data.frame(DGEList.gene[, order(DGEList.gene[which(rownames(DGEList.gene) == gene), ]) ]); colnames(DGEList.gene.sorted) <- "CPM"

    DGEsorted.counts <- data.frame(row.names = row.names(DGEList.gene.sorted))
    DGEsorted.counts$samplenumber <- 1:nrow(DGEsorted.counts)
    DGEsorted.counts$CPM <- DGEList.gene.sorted$CPM

    # assign group by Gaussian distribution
    library(mclust)
    fit <- Mclust(DGEsorted.counts$CPM, G=number_of_groups)
    DGEsorted.counts$group <- fit$classification
    # row.temp <- row.names(p.values[p.values$pval == min(p.values[,1]), ])
    # p.values$group <- c(rep(1, row.temp[length(row.temp)]),
    #  rep(0 , nrow(p.values) - as.numeric(row.temp[length(row.temp)])))
    # plot(p.values[4:5], col="black", pch=ifelse(p.values$group==1, 19, 21))

    if (design_input != 1) {
        design_input_matched <- design_input[match(row.names(DGEsorted.counts),row.names(design_input)), ]

        plot(DGEsorted.counts$samplenumber, DGEsorted.counts$CPM, xlim = c(-(max(nchar(colnames(design_input_matched))) / 2), nrow(DGEsorted.counts)), ylim = c(min(DGEsorted.counts$CPM) - ((max(nchar(row.names(DGEsorted.counts))) + ncol(design_input_matched)) * 0.3), max(DGEsorted.counts$CPM)),  col="black", xlab = "", ylab = "CPM (log2)", yaxt = 'n', xaxt='n',
        pch = ifelse(DGEsorted.counts$group > 2, 1, ifelse(DGEsorted.counts$group == 1, 4, 19)))
    } else {
        plot(DGEsorted.counts$samplenumber, DGEsorted.counts$CPM,  col="black", xlab = "", ylab = "CPM (log2)", axes = F)
    }

    if ((length(grep("ENSG", gene_name)) > 0) == T) {
        title(ensg2hugo(gene_name))
    } else {
        title(gene_name)
    }

    # mtext(paste(round(countMatches(1, DGEsorted.counts$group) / nrow(DGEsorted.counts) * 100, 1), "%", round(countMatches(2, DGEsorted.counts$group) / nrow(DGEsorted.counts) * 100, 1), "%", round(countMatches(3, DGEsorted.counts$group) / nrow(DGEsorted.counts) * 100, 1), "%"), side = 3, cex = 0.4)
    par(las=2)

    axis(1, at=1:nrow(DGEsorted.counts), pos = min(DGEsorted.counts$CPM) - 0.5, labels=row.names(DGEsorted.counts), cex.axis = cex.axis)
    axis(2, at=round(min(DGEsorted.counts$CPM)):round(max(DGEsorted.counts$CPM)), cex.axis = 0.5)

    if (design_input != 1) {
        design_input_matched <- design_input[match(row.names(DGEsorted.counts),row.names(design_input)), ]

        print("Possible values for design matrix is: 'YES', 'NO', 'AMPLIFIED', 'NON-AMPLIFIED', MALE', 'FEMALE', 'LOW', 'INTERMEDIATE', 'HIGH', '1', '2', '3', '4', '4S', 'UNKOWN', 'OTHER', 'NA' or empty values")
        ## Add color labeling of clinical characteristics

        # initiliaze coloring
        square_pos = min(DGEsorted.counts$CPM)-(max(nchar(row.names(DGEsorted.counts)))/2) - distance_between_axis_and_labels

        for (i in 0:(ncol(design_input_matched)-1)) {
            for (x in 1:nrow(design_input_matched)) {
                if (design_input_matched[x, (i+1)] == "FEMALE") {
                    color = "red"
                } else if (design_input_matched[x, (i+1)] == "MALE") {
                    color = "blue"

                } else if (design_input_matched[x, (i+1)] == "YES") {
                    color = "green"
                } else if (design_input_matched[x, (i+1)] == "NO") {
                    color = "red"

                } else if (design_input_matched[x, (i+1)] == 'AMPLIFIED') {
                    color = "red"
                } else if (design_input_matched[x, (i+1)] == 'NON-AMPLIFIED') {
                    color = "green"


                } else if (design_input_matched[x, (i+1)] == "LOW") {
                    color = "green"
                } else if (design_input_matched[x, (i+1)] == "INTERMEDIATE") {
                    color = "gold"
                } else if (design_input_matched[x, (i+1)] == "HIGH") {
                    color = "red"

                } else if (design_input_matched[x, (i+1)] == "1") {
                    color = "green"
                } else if (design_input_matched[x, (i+1)] == "2") {
                    color = "darkgreen"
                } else if (design_input_matched[x, (i+1)] == "3") {
                    color = "gold"
                } else if (design_input_matched[x, (i+1)] == "4") {
                    color = "red"
                } else if (design_input_matched[x, (i+1)] == "4S") {
                    color = "blue"

                } else if (design_input_matched[x, (i+1)] == "UNKOWN") {
                    color = "grey"
                } else if (design_input_matched[x, (i+1)] == "OTHER") {
                    color = "grey"
                } else if (design_input_matched[x, (i+1)] == "NA") {
                    color = "grey"

                } else {
                    color = "white"
                }
                # print(color)
                points(x, square_pos - (i * square_spacing), cex = square_size, pch = square_type, col = color)
            }
            # text(0 - nchar(colnames(design_input_matched)[i]), square_pos - (i * square_spacing), colnames(design_input_matched)[i], cex = square_text_size)
            text(0 - max(nchar(colnames(design_input_matched))) / 2, square_pos - (i * square_spacing), colnames(design_input_matched)[i+1], cex = label.cex)
        }

        design_input_matched$group <- DGEsorted.counts$group
        return(design_input_matched)
    }

    print(paste("Number of clusters:", max(fit$classification)))
    hugo2ensg(gene)
}

genewise_corr_plot_m <- function(df, gene_list, gene_sorted, cex.axis = 1, labels = c(0, "target", "ordered"), cex.labels = 1, target_gene_number = 1, gene_sorted_number = 1, identify_points_sorted = FALSE, identify_points_target = FALSE)
{
    df_genes <- data.frame(gene_list)
    par(mfrow=c((nrow(df_genes)/2),2))

    for (i in 1:nrow(df_genes)) {
        target_gene <- gene_list[i]
        genewise_corr_plot(df, target_gene, gene_sorted, cex.axis, labels)
    }
    par(mfrow=c(1,1))
}

pubmed.genelist.query <- function(genelist, additional.search.query = "") {
    library(RISmed)
    for (i in 1:length(genelist)) {
        res.temp <- EUtilsSummary(paste(genelist[i], additional.search.query), type="esearch", db="pubmed", datetype='pdat')
        if (i == 1) {
            f <- ArticleTitle(EUtilsGet(res.temp))
        } else {
            f <- rbind(f, data.frame(ArticleTitle(EUtilsGet(res.temp))))
        }
    }
    return(f)
}

pubmed.genelist.query <- function(genelist, additional.search.query = "") {
    res.temp[i] <- EUtilsSummary(paste(genelist[i], additional.search.query), type="esearch", db="pubmed", datetype='pdat')
    for (i in 1:length(genelist)) {
        res.temp <- EUtilsSummary(paste(genelist[i], additional.search.query), type="esearch", db="pubmed", datetype='pdat')
        if (i == 1) {
            f <- ArticleTitle(EUtilsGet(res.temp))
        } else {
            f <- rbind(f, data.frame(ArticleTitle(EUtilsGet(res.temp))))
        }
    }
    return(f)
}

overview.table.of.proteins.and.lncRNAs <- function(proteins.dir.or.files, lncRNA.dir.or.files, row.names_remove.text = "", row.names_prepend.to.numeric = "", p.value.cutoff = 0.05, max.genes = 20, file.html = "") {
    library(gtools)
    library(tableHTML)

    if (grepl(".txt", toString(proteins.dir.or.files)) == TRUE ) {
        protein.files <- proteins.dir.or.files
        lnc.files <- lncRNA.dir.or.files
    } else {
        protein.files <- paste(proteins.dir.or.files, list.files(proteins.dir.or.files), sep = "")
        lnc.files <- paste(lncRNA.dir.or.files, list.files(lncRNA.dir.or.files), sep = "")
    }

        for (i in 1:length(protein.files)) {

            proteins.edgeR.results <- read.delim(protein.files[i])
            proteins.edgeR.results.up <- proteins.edgeR.results[proteins.edgeR.results$logFC > 0, ]
            if (nrow(proteins.edgeR.results.up) > max.genes) proteins.edgeR.results.up <- proteins.edgeR.results.up[1:max.genes, ]
            proteins.edgeR.results.down <- proteins.edgeR.results[proteins.edgeR.results$logFC < 0, ]
            if (nrow(proteins.edgeR.results.down) > max.genes) proteins.edgeR.results.down <- proteins.edgeR.results.down[1:max.genes, ]
            proteins.edgeR.results <- rbind(proteins.edgeR.results.up, proteins.edgeR.results.down)
            # if (nrow(proteins.edgeR.results) > max.genes) proteins.edgeR.results <- proteins.edgeR.results[1:max.genes, ]

            lncRNAs.edgeR.results <- read.delim(lnc.files[i])
            lncRNAs.edgeR.results.up <- lncRNAs.edgeR.results[lncRNAs.edgeR.results$logFC > 0, ]
            if (nrow(lncRNAs.edgeR.results.up) > max.genes) lncRNAs.edgeR.results.up <- lncRNAs.edgeR.results.up[1:max.genes, ]
            lncRNAs.edgeR.results.down <- lncRNAs.edgeR.results[lncRNAs.edgeR.results$logFC < 0, ]
            if (nrow(lncRNAs.edgeR.results.down) > max.genes) lncRNAs.edgeR.results.down <- lncRNAs.edgeR.results.down[1:max.genes, ]
            lncRNAs.edgeR.results <- rbind(lncRNAs.edgeR.results.up, lncRNAs.edgeR.results.down)
            # if (nrow(lncRNAs.edgeR.results) > max.genes) lncRNAs.edgeR.results <- lncRNAs.edgeR.results[1:max.genes, ]

            row.name <- gsub(dirname(protein.files[i]), "",protein.files[i])
            row.name <- gsub("/", "", row.name)
            row.name <- gsub(row.names_remove.text, "", row.name)
            row.name <- gsub(".txt", "", row.name)
            row.name <- paste(row.names_prepend.to.numeric, row.name, sep = "")

            #lncRNAs.edgeR.results$col <- ifelse(lncRNAs.edgeR.results$logFC > 0,
            #        rgb(0, abs(lncRNAs.edgeR.results$logFC) / max(abs(lncRNAs.edgeR.results$logFC)),0),
            #        rgb(abs(lncRNAs.edgeR.results$logFC) / -min(lncRNAs.edgeR.results$logFC),0,0))

            proteins.edgeR.results$max <- ifelse(proteins.edgeR.results$logFC >= 0,
                                                 max(proteins.edgeR.results$logFC),
                                                 min(proteins.edgeR.results$logFC))
            proteins.edgeR.results$logFC.frac <- proteins.edgeR.results$logFC / proteins.edgeR.results$max
            proteins.edgeR.results$col <- ifelse(proteins.edgeR.results$logFC >= 0,
                                                 rgb(0, proteins.edgeR.results$logFC.frac,0),
                                                 rgb(proteins.edgeR.results$logFC.frac,0,0))

            lncRNAs.edgeR.results$max <- ifelse(lncRNAs.edgeR.results$logFC >= 0,
                                                max(lncRNAs.edgeR.results$logFC),
                                                min(lncRNAs.edgeR.results$logFC))
            lncRNAs.edgeR.results$logFC.frac <- lncRNAs.edgeR.results$logFC / lncRNAs.edgeR.results$max
            lncRNAs.edgeR.results$col <- ifelse(lncRNAs.edgeR.results$logFC >= 0,
                                                rgb(0, lncRNAs.edgeR.results$logFC.frac,0),
                                                rgb(lncRNAs.edgeR.results$logFC.frac,0,0))

            lncRNAs.edgeR.results$external_gene_name <- ifelse(lncRNAs.edgeR.results$FDR <= 0.001, paste(lncRNAs.edgeR.results$external_gene_name, "***", sep = ""),  as.character(lncRNAs.edgeR.results$external_gene_name))
            lncRNAs.edgeR.results$external_gene_name <- ifelse(lncRNAs.edgeR.results$FDR <= 0.01 & lncRNAs.edgeR.results$FDR > 0.001, paste(lncRNAs.edgeR.results$external_gene_name, "**", sep = ""),  as.character(lncRNAs.edgeR.results$external_gene_name))
            lncRNAs.edgeR.results$external_gene_name <- ifelse(lncRNAs.edgeR.results$FDR <= 0.05 & lncRNAs.edgeR.results$FDR > 0.01, paste(lncRNAs.edgeR.results$external_gene_name, "*", sep = ""),  as.character(lncRNAs.edgeR.results$external_gene_name))
            proteins.edgeR.results$external_gene_name <- ifelse(proteins.edgeR.results$FDR <= 0.001, paste(proteins.edgeR.results$external_gene_name, "***", sep = ""),  as.character(proteins.edgeR.results$external_gene_name))
            proteins.edgeR.results$external_gene_name <- ifelse(proteins.edgeR.results$FDR <= 0.01 & proteins.edgeR.results$FDR > 0.001, paste(proteins.edgeR.results$external_gene_name, "**", sep = ""),  as.character(proteins.edgeR.results$external_gene_name))
            proteins.edgeR.results$external_gene_name <- ifelse(proteins.edgeR.results$FDR <= 0.05 & proteins.edgeR.results$FDR > 0.01, paste(proteins.edgeR.results$external_gene_name, "*", sep = ""),  as.character(proteins.edgeR.results$external_gene_name))

            proteins.edgeR.results$external_gene_name <- paste('<font color="', proteins.edgeR.results$col, '">', proteins.edgeR.results$external_gene_name, '</font>', sep = "")
            lncRNAs.edgeR.results$external_gene_name <- paste('<font color="', lncRNAs.edgeR.results$col, '">', lncRNAs.edgeR.results$external_gene_name, '</font>', sep = "")

            p.lnc <- lncRNAs.edgeR.results
            p <- proteins.edgeR.results

            p.lnc <- p.lnc[p.lnc$FDR < p.value.cutoff, ]
            p <- p[p$FDR < p.value.cutoff, ]

            if ( i == 1 ) {
                t2 <- data.frame( lncRNAs.up = toString(head(p.lnc[p.lnc$logFC > 0 , "external_gene_name"], 20)),
                                  lncRNAs.down = toString(head(p.lnc[p.lnc$logFC < 0 , "external_gene_name"], 20)),
                                  proteins.up = toString(head(p[p$logFC > 0 , "external_gene_name"], 20)),
                                  proteins.down = toString(head(p[p$logFC < 0 , "external_gene_name"],20)), row.names = row.name)
                colnames(t2)
            } else {
                t <- data.frame( lncRNAs.up = toString(head(p.lnc[p.lnc$logFC > 0 , "external_gene_name"], 20)),
                                 lncRNAs.down = toString(head(p.lnc[p.lnc$logFC < 0 , "external_gene_name"], 20)),
                                 proteins.up = toString(head(p[p$logFC > 0 , "external_gene_name"],20)),
                                 proteins.down = toString(head(p[p$logFC < 0 , "external_gene_name"],20)), row.names = row.name)
                t2 <- rbind(t2, t)
            }
        }

    # &nbsp; designates space in HTML output
    colnames(t2) <- c(paste("lncRNAs", '&nbsp;', "up", sep = ""),
                      paste("lncRNAs", '&nbsp;', "down", sep = ""),
                      paste("proteins", '&nbsp;', "up", sep = ""),
                      paste("proteins", '&nbsp;', "down", sep = ""))

    t2 <- t2[mixedorder(row.names(t2)), ]

    tableHTML::tableHTML(t2, theme = "scientific", footer = "*** FDR < 0.005, ** FDR < 0.001, * FDR < 0.05")

        if (file.html != "") {
            tryCatch(system(paste("trash", file.html), ignore.stdout = T, ignore.stderr = T))
            tryCatch(tableHTML::write_tableHTML(tableHTML::tableHTML(t2, theme = "scientific", footer = "*** FDR < 0.005, ** FDR < 0.001, * FDR < 0.05"), file.html))
            system(paste("open", file.html))
        }

    return(t2)
}

rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

get_args <- function () {
    as.list( match.call(
        def = sys.function( -1 ),
        call = sys.call(-1)) )[-1]

}

boxplot_pp_row <- function(df, row){
    boxplot(df[row,1:5], df[row,6:10], df[row,11:15], df[row,16:20],
    names = c("cell_pre", "cell_post", "TARGET_pri", "TARGET_rec"))
    gene_name <- ensg2hugo(row.names(top_pp_TARGET)[row])[1,2]
    gene_biotype <- ensg2gene_biotype(row.names(top_pp_TARGET)[row])[1,2]
    title(main = paste(gene_name, " ", gene_biotype))
}

boxplot_pp_gene <- function(df, gene){
    ensg <- hugo2ensg(gene)[1,2]
    df.temp <- data.frame(df[ensg,])
    df.temp[1:5,]
    boxplot(df.temp[1:5,], df.temp[6:10, 1], df.temp[11:15, 1], df.temp[16:20, 1],
    names = c("cell_pre", "cell_post", "TARGET_pri", "TARGET_rec"))
    gene_name <- gene
    gene_biotype <- ensg2gene_biotype(ensg)[1,2]
    title(main = paste(gene_name, " ", gene_biotype))
}

plot.corr <- function(df, gene1, gene2, cex.axis = 1) {
    # df = mir_mrna_exp; gene1 = "ENSG00000220635"; gene2 = "mir.490"
    # rm(df); rm(gene1); rm(gene2)

    gene1 <- subset(mir_mrna_exp, row.names(mir_mrna_exp) == gene1)
    gene2 <- subset(mir_mrna_exp, row.names(mir_mrna_exp) == gene2)
    gene1 <- log2(gene1 + 1)
    gene2 <- log2(gene2 + 1)

    Rval <- cor(t(gene1), t(gene2), method = c("pearson"))

    if ((length(grep("ENSG", row.names(Rval))) > 0) == T) {
        name_one <- ensg2hugo(row.names(Rval))[1,2]
        ensg2hugo(row.names(Rval))
        biotype_one <- ensg2gene_biotype(row.names(Rval))
    } else {
        name_one <- row.names(Rval)
    }

    if ((length(grep("ENSG", colnames(Rval))) > 0) == T) {
        name_two <- ensg2hugo(colnames(Rval))[1,2]
        ensg2hugo(colnames(Rval))
        biotype_two <- ensg2gene_biotype(colnames(Rval))
    } else {
        name_two <- colnames(Rval)
    }

    par(mar=c(5, 4, 4, 5) + 0.1)

    plot(t(gene1), col = "blue", pch = 2, cex = 1, xaxt='n', ylab=name_one, xlab = "")
    par(new=TRUE)
    plot(t(gene2), xaxt='n', yaxt='n', ylab = "", xlab = "", col = "red", pch = 6)
    mtext(name_two,side=4,col="red", line=3, las = 3)
    axis(4, at=round(as.numeric(gene2)), col = "black", col.axis="red")
    # mtext(side = 1, at=colnames(gene1))
    axis(1, at=1:ncol(gene1), labels=colnames(gene1), cex.axis = cex.axis, las = 2)

    if(exists("biotype_one") == T) {
        print(paste(name_one, " " ,biotype_one))
    }

    if(exists("biotype_two") == T) {
        print(paste(name_two, " " ,biotype_two))
    }

    title(paste("R = ", round(Rval[1], 2)))
}

plot.clusters <- function(df){
    wss <<- (nrow(df)-1)*sum(apply(df,2,var))
    for (i in 2:15) wss[i] <- sum(kmeans(df, centers=i)$withinss)
    plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
}

K_means_clustering <- function(df, number_of_clusters){
    fit <<- kmeans(df, number_of_clusters) # 5 cluster solution
    # get cluster means
    agg <<- aggregate(df,by=list(fit$cluster),FUN=mean)
    # append cluster assignment
    df2 <- data.frame(df, fit$cluster)
    return(df2[order(df2$fit.cluster),])
}

read.gmt2 <- function(file, Title = "", Footnote = "") {
    file.in <- scan(file, what="", sep="\n")
    # # Separate elements by one or more tabs
    file.temp <- strsplit(file.in, "\t")
    names(file.temp) <- sapply(file.temp, `[[`, 1)
    file.temp <- lapply(file.temp, `[`, -1)
    lapply(file.temp, `[`, -1)
}

#' @title table2display
#' @param input table1
#' @export
#' @import gtable
#' @examples table2display()
table2display <- function(table1, Title = NULL, Footnote = NULL, display.rownames = T)
{
    if (display.rownames == T) table <- tableGrob(table1, theme= ttheme_minimal())
    if (display.rownames == F) table <- tableGrob(table1, theme= ttheme_minimal(), rows = NULL)

    title <- textGrob(Title,gp=gpar(fontsize=14))
    footnote <- textGrob(Footnote, x=0, hjust=0,
    gp=gpar( fontface="italic", fontsize=10))

    padding <- unit(0.5,"line")
    table <- gtable_add_rows(table,
    heights = grobHeight(title) + padding,
    pos = 0)
    table <- gtable_add_rows(table,
    heights = grobHeight(footnote)+ padding)
    table <- gtable_add_grob(table, list(title, footnote),
    t=c(1, nrow(table)), l=c(1,2),
    r=ncol(table))
    grid.newpage()
    grid.draw(table)
}

d3heatmap.romer <- function(Qlnc_name, cluster1, cluster2, plot.cluster = "both")
{
    test <- entrez2ensg(cluster1)
    test2 <- entrez2ensg(cluster2)
    genes <- subset(cpm(Qlnc_name, log = T), row.names(Qlnc_name) %in% test$ensembl_gene_id)
    genes2 <- subset(cpm(Qlnc_name), row.names(Qlnc_name) %in% test2$ensembl_gene_id)
    genes3 <- rbind(genes, genes2)
    genes2.new.names <- ensg2ext_name(ensg = row.names(genes3))
    genes4 <<- merge(genes3, genes2.new.names, by.x="row.names", by.y="ensembl_gene_id")
    row.names(genes4) <- make.names(genes4$external_gene_name, unique = T)
    genes4$external_gene_name <- NULL
    genes4$Row.names <- NULL
    test3 <- ensg2ext_name(test)
    test4 <- ensg2ext_name(test2)

    if (plot.cluster == "both"){
        d3heatmap(genes4)
    } else if (plot.cluster == "first") {
        d3heatmap(genes4, row.names(genes4) %in% test3$external_gene_name)
    } else if (plot.cluster == "second") {
        d3heatmap(genes4, row.names(genes4) %in% test4$external_gene_name)
    }

}

heatmap.3.1 <- function(x,
Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
distfun = dist,
hclustfun = hclust,
dendrogram = c("both","row", "column", "none"),
symm = FALSE,
scale = c("none","row", "column"),
na.rm = TRUE,
revC = identical(Colv,"Rowv"),
add.expr,
breaks,
symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
col = "heat.colors",
colsep,
rowsep,
sepcolor = "white",
sepwidth = c(0.05, 0.05),
cellnote,
notecex = 1,
notecol = "cyan",
na.color = par("bg"),
trace = c("none", "column","row", "both"),
tracecol = "cyan",
hline = median(breaks),
vline = median(breaks),
linecol = tracecol,
margins = c(5,5),
ColSideColors,
RowSideColors,
side.height.fraction=0.3,
cexRow = 0.2 + 1/log10(nr),
cexCol = 0.2 + 1/log10(nc),
labRow = NULL,
labCol = NULL,
key = TRUE,
keysize = 1.5,
density.info = c("none", "histogram", "density"),
denscol = tracecol,
symkey = max(x < 0, na.rm = TRUE) || symbreaks,
densadj = 0.25,
main = NULL,
xlab = NULL,
ylab = NULL,
lmat = NULL,
lhei = NULL,
lwid = NULL,
ColSideColorsSize = 1,
RowSideColorsSize = 1,
KeyValueName="Value",...){

    invalid <- function (x) {
        if (missing(x) || is.null(x) || length(x) == 0)
        return(TRUE)
        if (is.list(x))
        return(all(sapply(x, invalid)))
        else if (is.vector(x))
        return(all(is.na(x)))
        else return(FALSE)
    }

    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
    "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
    "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE

    marttest =  useMart("ENSEMBL_MART_ENSEMBL", dataset='hsapiens_gene_ensembl', host="jul2016.archive.ensembl.org")
    rm()
    hugo2ensg("MYCN", mart)

    if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
        c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
            dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
            dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
        c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
            dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
            dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
    (1:nr)[rowInd]
    else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
    (1:nc)[colInd]
    else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
        breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
        length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
    col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)

        if (!missing(ColSideColors)) {
            #if (!is.matrix(ColSideColors))
            #stop("'ColSideColors' must be a matrix")
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
            stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            #lhei <- c(lhei[1], 0.2, lhei[2])
            lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
        }

        if (!missing(RowSideColors)) {
            #if (!is.matrix(RowSideColors))
            #stop("'RowSideColors' must be a matrix")
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
            stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
            #lwid <- c(lwid[1], 0.2, lwid[2])
            lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }

    if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
            par(mar = c(margins[1], 0, 0, 0.5))
            image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[,rowInd, drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            if (length(rownames(RowSideColors)) > 0) {
                axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
            }
        }
    }

    if (!missing(ColSideColors)) {

        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[2]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc = ColSideColors[colInd, , drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            image(csc, col = as.vector(csc.colors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
            }
        }
    }

    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
        ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
    retval$rowDendrogram <- ddr
    if (exists("ddc"))
    retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
        col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
    cex.axis = cexCol)
    if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
    cex.axis = cexRow)
    if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
    eval(substitute(add.expr))
    if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
    col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }

        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
        xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
        mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
        mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
    high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}

countMatchesBetweenTwoCols <- function(match_col, lookup_col) {
    t3 <- data.frame(nr_matches = as.numeric(), value = as.numeric())
    t3 <- as.data.table(t3)
    
    for (i in 1:length(match_col)) {
        # t3[i, 1] <- countMatches(match_col[i], lookup_col)
        # t3[i, 2] <- match_col[i]
        t3 <- rbindlist(list(t3, as.list(c(countMatches(match_col[i], lookup_col), match_col[i]))))
    
    }
    t3 <- t3[!duplicated(t3$value), ]
    t3 <- t3[order(-t3$nr_matches), ]
    return(t3)
}

plot.distribution.of.two.samples <- function(input.sample1, input.sample2, subtitle = "", type = "n", name.sample1 = "", name.sample2 = "", background.col.1 = "transparent", line.col.1 = "red", background.col.2 = "transparent", line.col.2 = "blue", paired.t.test = T) {
    quartz()

    d <- density(input.sample1)
    f <- density(input.sample2)
    min <- ifelse((min(f$x) > min(d$x)) == T, min(d$x), min(f$x))
    max <- ifelse((max(f$x) < max(d$x)) == T, max(d$x), max(f$x))

    name.sample1 <- ifelse(name.sample1 == "", deparse(substitute(input.sample1)), "sample1")
    name.sample2 <- ifelse(name.sample2 == "", deparse(substitute(input.sample2)), "sample2")

    subtitle <- ifelse(subtitle == "",
    paste("N (", name.sample1, ") = ", d$n, "\t Bandwidth (", name.sample1, ") = ", round(d$bw,3),
    "\n N (", name.sample2, ") = ", f$n, "\t Bandwidth (", name.sample2, ") = ", signif(f$bw,3), sep = ""),
    subtitle)
    plot(d, type=type, xlim = c(min, max), main = "", xlab = subtitle)
    polygon(d, col=background.col.1, border=line.col.1)
    polygon(f, col=background.col.2, border=line.col.2)

    t.test <- t.test(input.sample1, input.sample2, paired = paired.t.test)
    wilcox.test <- wilcox.test(input.sample1, input.sample2)

    title(
    paste(
    name.sample1, " vs ", name.sample2,
    "\n T.test pval = ", format.pval(t.test$p.value, scientific = T, digits = 3),
    "\t R (pearson) = ", format.pval(cor(input.sample1, input.sample2, method = "pearson"), scientific = T, digits = 3),
    "\n Wilcox pval = ", format.pval(wilcox.test$p.value, scientific = T, digits = 3),
    "\t R (spearman) = ",
    format.pval(cor(input.sample1, input.sample2, method = "spearman"), scientific = T, digits = 3),
    sep = ""
    )
    )
}

venn.jaccard <- function(input) {
    venn(input)
    title(paste(
    "N (", names(input)[[1]], ") = ", length(input[[1]]),
    "\nN (", names(input)[[2]], ") = ", length(input[[2]]),
    "\nJaccard index = ", round(table(input[[1]] %in% input[[2]])[2] / (length(input[[1]]) + length(input[[2]])),2), sep = ""))
}



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots==1) {
        print(plots[[1]])

    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

tot_count <- function(x, labs, digits, varlen)
{
  paste(labs, "\n\nn =", x$frame$n)
}

only_count <- function(x, labs, digits, varlen)
{
  paste(x$frame$n)
}

## pull the standardized error from a forest object (randomForestSRC)
    get.error <- function(obj) {
      100 * c(sapply(obj$yvar.names, function(nn) {
        o.coerce <- randomForestSRC:::coerce.multivariate(obj, nn)
        if (o.coerce$family == "class") {
          tail(o.coerce$err.rate[, 1], 1)
        }
        else {
          tail(o.coerce$err.rate, 1) / var(o.coerce$yvar, na.rm = TRUE)
} }))
}
    ## pull the standardized VIMP from a forest object (randomForestSRC)
    get.vimp <- function(obj) {
      vimp <- 100 * do.call(cbind, lapply(obj$yvar.names, function(nn) {
        o.coerce <- randomForestSRC:::coerce.multivariate(obj, nn)
        if (o.coerce$family == "class") {
          o.coerce$importance[, 1]
        }
        else {
          o.coerce$importance / var(o.coerce$yvar, na.rm = TRUE)
} }))
      colnames(vimp) <- obj$yvar.names
vimp }

surv.plot <- function(df, time.col.name, event.col.name, variable.col.name, covariate.col.name = "", subset.variables = "", rho = 0, censor_level = 2, legend.text.cex = 10, names = "", plot.pval = F, plot.only.logrank.pval = T, validate.coxph = F, pval.text.size = 10) {
    library(survival)

    # df <- k; time.col.name = "Overall.Survival.Time.in.Days"; event.col.name = "Vital.Status"; variable.col.name = "group"; covariate.col.name = ""; subset.variables = ""; rho = 0; censor_level = 2; legend.text.cex = 10; names = ""; plot.pval = T; plot.only.logrank.pval = F; validate.coxph = F

    z <- cbind(df[time.col.name], df[event.col.name], df[variable.col.name])
    colnames(z) <- c("time", "event", "variable")
    if (covariate.col.name != "") z <- cbind(df[time.col.name], df[event.col.name], df[variable.col.name], df[covariate.col.name])
    if (covariate.col.name != "") colnames(z) <- colnames(z) <- c("time", "event", "variable", "covariable")
    tryCatch(if(subset.variables != "") z <- z[z$variable %in% subset.variables, ])

    type = levels(factor(z$variable))
    type.covariable = levels(factor(z$covariable))

    if(class(z$time) %in% c("factor", "character")) {
        z$time <- factor(z$time)
        } else {
        z$time <- as.numeric(z$time)
        }

    if(class(z$event) %in% c("factor", "character")) {
        z$event <- factor(z$event)
        } else {
        z$event <- as.numeric(z$event)
        }

    if(class(z$variable) %in% c("factor", "character")) {
        z$variable <- factor(z$variable)
        } else {
        z$variable <- as.numeric(z$variable)
        }

    if(covariate.col.name != "" & class(z$variable) %in% c("factor", "character")) {
        z$covariable <- factor(z$covariable)
        } else if (covariate.col.name != "" & !class(z$variable) %in% c("factor", "character")) {
        z$covariable <- as.numeric(z$covariable)
        }


    cat(paste(levels(z$event)[censor_level], "will be used as event"))

    # while(length(z$time) - length(unique(z$time)) != 0) z$time <- z$time + (as.numeric(as.factor(duplicated(z$time)))-1)
    if (covariate.col.name == "") {

            surv_survfit = survfit(Surv(time, event == levels(factor(z$event))[censor_level]) ~ variable, data = z)
            surv_survdiff <- survdiff(Surv(time, event == levels(factor(z$event))[censor_level]) ~ variable, data = z, rho = rho)
            surv_coxphfit <- coxph(Surv(time, event == levels(factor(z$event))[censor_level]) ~ variable, data = z, robust = T)

                if (names == "") {

                    names(surv_survfit$strata)  <- gsub("variable=", "", names(surv_survfit$strata))
                    names(surv_survdiff$n) <- gsub("variable=", "", names(surv_survdiff$n))
                    names(surv_coxphfit$coefficients) <- gsub("variable", "" , names(surv_coxphfit$coefficients))

                    } else {
                    names(surv_survfit$strata) <- names
                    names(surv_survdiff$n) <- names

                    }

            cat("surv_survdiff\n")
            print(surv_survdiff)
            cat("summary(surv_coxphfit)\n")
            print(summary(surv_coxphfit))

        } else {
            surv_survfit = survfit(Surv(time, event == levels(factor(z$event))[censor_level]) ~ variable + covariable, data = z)
            surv_survdiff <- survdiff(Surv(time, event == levels(factor(z$event))[censor_level]) ~ variable + covariable, data = z, rho = rho)
            surv_coxphfit <- coxph(Surv(time, event == levels(factor(z$event))[censor_level]) ~ variable + covariable, data = z, robust = T)

            if (names == "") {

                    names(surv_survfit$strata)  <- gsub("covariable=", "", names(surv_survfit$strata))
                    names(surv_survfit$strata)  <- gsub("variable=", "", names(surv_survfit$strata))
                    names(surv_survdiff$n) <- gsub("covariable=", "", names(surv_survdiff$n))
                    names(surv_survdiff$n) <- gsub("variable=", "", names(surv_survdiff$n))
                    names(surv_coxphfit$coefficients) <- gsub("covariable", "" , names(surv_coxphfit$coefficients))
                    names(surv_coxphfit$coefficients) <- gsub("variable", "" , names(surv_coxphfit$coefficients))

                    } else {
                    names(surv_survfit$strata) <- names
                    names(surv_survdiff$n) <- names

                    }

            cat("surv_survdiff\n")
            print(surv_survdiff)
            cat("summary(surv_coxphfit)\n")
            print(summary(surv_coxphfit))

        # pval <- summary(surv_coxphfit)
        # pval <- data.frame(pval$coefficients)
        # pval <- round(pval[5],3)
        # colnames(pval) <- "pval"
        # names(surv_coxphfit$coefficients) <- c(paste(type[2], type[1], sep = "-"),
        # paste(type.covariable[2], type.covariable[1], sep = "-"))
        # tryCatch(names(surv_survfit$strata) <- c(paste(type[1], type.covariable[1], sep = ", "),
        # paste(type[1], type.covariable[2], sep = ", "),
        # paste(type[2], type.covariable[1], sep = ", "),
        # paste(type[2], type.covariable[2], sep = ", ")))
        # tryCatch(names(surv_survdiff$n) <- c(paste(type[1], type.covariable[1], sep = ", "),
        # paste(type[1], type.covariable[2], sep = ", "),
        # paste(type[2], type.covariable[1], sep = ", "),
        # paste(type[2], type.covariable[2], sep = ", ")))
        }

    library(broom)
    surv_tidy = tidy(surv_survfit)
    mx = max(surv_tidy$n.censor)

    # plot(survfit(surv_coxphfit))
    # surv_survfit
    # surv_survdiff
    # nrow(z[z$variable ==1 & z$event == 2, ])
    # nrow(z[z$variable ==2 & z$event == 2, ])
    # median(z[z$variable ==1, "time"])
    # median(z[z$variable ==2, "time"])
    # median(z[z$variable ==1 & z$event == 2, "time"])
    # median(z[z$variable ==2 & z$event == 2, "time"])

    # surv_tidy$strata <- gsub(levels(factor(surv_tidy$strata))[1], paste(type[1], " (", table(surv_tidy$strata)[1], ")", sep = ""), surv_tidy$strata)
    # surv_tidy$strata <- gsub(levels(factor(surv_tidy$strata))[2], paste(type[2], " (", table(surv_tidy$strata)[2], ")", sep = ""), surv_tidy$strata)

    cat("tidy(surv_coxphfit)\n")
    tidy(surv_coxphfit)


    library(ggthemes)
    p <- ggplot(surv_tidy, aes(time, estimate, fill = strata), geom="step") +
        geom_step(aes(linetype =  factor(strata), color = factor(strata)), size = 0.5) +
        geom_point(shape = "+", size = 5, show.legend = F) +
      geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.1, show.legend = F) +
      scale_linetype_discrete(name = "") +
      xlab("OS (days)") +
      ylab("Proportion Survival") +
      theme_classic() +
        theme(legend.position = "bottom", legend.text = element_text(size = legend.text.cex)) +
      # theme(legend.justification = c(1,1), legend.position = c(1,1), legend.text = element_text(size = 10), legend.key.size = unit(1, "cm"), legend.background = ) +
        guides(linetype=guide_legend(override.aes=list(size=0.25))) +
        scale_color_discrete("")

    pval <- data.frame(glance(surv_coxphfit))

        if (plot.only.logrank.pval == T) {
            pval <- paste("p=", signif(pval$p.value.sc, 3), sep = "") # Logrank test
        } else {
            pval <- paste(paste("p (robust) =", signif(pval$p.value.robust, 3)), # Robust
            paste("p (logrank) =", signif(pval$p.value.sc, 3)), # Logrank test
            paste("p (LR) =", signif(pval$p.value.log, 3)), sep = "\n") # Likelihood ratiotest
        }

            if (plot.pval == T & covariate.col.name == "") {
                p + geom_text(aes(x=Inf,y=-Inf, hjust = 1, vjust = -1 , label = pval), size = pval.text.size)
            } else {
                p
            }

    # if(covariate.col.name == "" && validate.coxph == T) {
        #validate_coxph = cox.zph(surv_coxphfit, transform = "rank")
        #print(validate_coxph)
        #plot(validate_coxph)
    # }
}

plot.venn <- function(input, text = "", file = "", elements.on.each.line = 10, display.common.names = T) {
# setRepositories(); 1 8 9 6 2 5 # http://matticklab.com/index.php?title=Weighted_Venn_diagrams_in_R
# chooseBioCmirror(); 3
# install.packages("Vennerable", repos="http://R-Forge.R-project.org")
library(Vennerable)

input.venn <- Vennerable::Venn(input)
C3 <- compute.Venn(input.venn, doWeights = TRUE, doEuler = T)
gp <- VennThemes(C3)

tryCatch(gp$Set$Set1$lwd <- 0); tryCatch(gp$Set$Set2$lwd <- 0); tryCatch(gp$Set$Set3$lwd <- 0); tryCatch(gp$Set$Set4$lwd <- 0); tryCatch(gp$Set$Set5$lwd <- 0); tryCatch(gp$Set$Set6$lwd <- 0); tryCatch(gp$Set$Set7$lwd <- 0); tryCatch(gp$Set$Set8$lwd <- 0); tryCatch(gp$Set$Set9$lwd <- 0)
tryCatch(gp$SetText$Set1$col <- "black"); tryCatch(gp$SetText$Set2$col <- "black"); tryCatch(gp$SetText$Set3$col <- "black"); tryCatch(gp$SetText$Set4$col <- "black"); tryCatch(gp$SetText$Set5$col <- "black"); tryCatch(gp$SetText$Set6$col <- "black"); tryCatch(gp$SetText$Set7$col <- "black"); tryCatch(gp$SetText$Set8$col <- "black"); tryCatch(gp$SetText$Set9$col <- "black")

    if(file == "") {
        plot(C3, gpList = gp)
    } else {
        pdf(file, onefile = F)
        plot(C3, gpList = gp)
        dev.off()
    }

    if (display.common.names == T) {

        if (text == "") text <- Reduce(intersect, input)

        if (length(text) < elements.on.each.line) elements.on.each.line = length(text)

        for (i in 1:floor(length(text) / elements.on.each.line)) {
        assign(paste("text", i, sep  = "."), text[(elements.on.each.line*(i-1)+1):(elements.on.each.line*i)])
        }

        if (length(text) %% elements.on.each.line != 0) {
            i = i + 1
            assign(paste("text", i, sep  = "."), text[(length(text)-length(text) %% elements.on.each.line+1):length(text)])
        }

        for (k in 1:i) {
                    text <- get(paste("text", k, sep  = "."))
                    grid.text(toString(text), x = 0.04, y = 0.15 + ((k-1)*0.02), gp=gpar(fontsize=8), just = c("left", "bottom"))
        }

    }
}

write.list <- function(x, file.txt, sep.names = "\t", sep.val = "\t") {
    z <- deparse(substitute(x))
    cat(z, "\n", file=file.txt)

    nams=names(x)

        if (sep.names == sep.val) {
                for (i in seq_along(x) ) {
                    cat(nams[i], sep.names,  x[[i]], "\n",
                    file=file.txt, append=TRUE)
                }

            } else {

            }
       # awk '{gsub(" \t ","\t",$0); print;}' /Volumes/Peter/gene_models/ensembl.gene.ontology.golist.txt | sed 's/ /, /g' > /Volumes/Peter/gene_models/ensembl.gene.ontology.golist2.txt
}

read.myGO2genes <- function(file) {
myGO2genes <- read.delim(file)
row.names(myGO2genes) <- myGO2genes$ensembl_gene_id; myGO2genes$ensembl_gene_id <- NULL
myGO2genes.list <- list()
    for (x in 1:nrow(myGO2genes)) {
    t <- myGO2genes[x,]
    list <- t[t != 0]
    myGO2genes.list[[row.names(t)]] <- list
    return(myGO2genes.list)
    }
}

plot.diagnostic.pvalue.plot <- function(pvalues, pvalue.column.name = "", overlay.padj = F, padj.column.name = "", binwidth = 0.01, xlim.max = 1, lookup.intervall = F, add.y.line = F, title = "") {
    if ( pvalue.column.name == "") {

            if(grepl("$", strsplit(as.character(match.call(plot.diagnostic.pvalue.plot, envir = parent.frame())), "=")[[2]], fixed = T) == F) {
                if (title == "") title = strsplit(as.character(match.call(plot.diagnostic.pvalue.plot, envir = parent.frame())), "=")[[2]]
                df <- data.frame(pvalues)
                colnames(df) <- "pvalue"
                pvalue.column.name <- "pvalue"
            } else {
            pvalue.column.name <- gsub(".*\\$","", strsplit(as.character(match.call(plot.diagnostic.pvalue.plot, envir = parent.frame())), "=")[[2]])

            df <- get(gsub("\\$.*","", strsplit(as.character(match.call(plot.diagnostic.pvalue.plot, envir = parent.frame())), "=")[[2]]  ))
            df <- as.data.frame(df)
            if (title == "") title = gsub("\\$.*","", strsplit(as.character(match.call(plot.diagnostic.pvalue.plot, envir = parent.frame())), "=")[[2]]  )
            }
    } else {
        df <- as.data.frame(pvalues)
    }


            if (overlay.padj == T & padj.column.name == "") padj.column.name <- colnames(df)[grep(c("fdr|bonf|holm|hoch|hommel|padj|p.adj|p_adj|fwer"), tolower(colnames(df)))]

            if (overlay.padj == F) {
                df <- df[c(pvalue.column.name)]
                colnames(df) <- "pvalue"
            } else {
                df <- df[c(pvalue.column.name, padj.column.name)]
                colnames(df) <- c("pvalue", "padj")
            }


    xlim <- c(0, xlim.max)

     print(paste("Number of tests found in interval:", length(df$pvalue < xlim.max)))

            if (overlay.padj == F) {
            p <- ggplot(df, aes(x = pvalue)) + geom_histogram(binwidth = binwidth, boundary = 0) + xlim(xlim) + theme_classic() + ggtitle(label = title) + theme(text = element_text(size = 15),  plot.title = element_text(hjust = 0.5))
            } else {
            p <- ggplot(df, aes(x = pvalue)) + geom_histogram(binwidth = binwidth, boundary = 0) + xlim(xlim) + theme_classic() +
            geom_histogram(aes(x = padj), data = df, binwidth = binwidth, boundary = 0, fill = alpha("red", 0.5)) + ggtitle(label = title) + theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5))
            }

            if (add.y.line == T) p <- p + geom_line(y = length(df$pvalue < xlim.max) / (xlim.max / binwidth), col = "red")

  return(p)

    if (lookup.intervall == T) {
  cat("Number of tests in each interval:\n")
  tryCatch(e <- t(sapply(by(as.list(df[1]), as.list(df[1]), cut, breaks = seq(0, xlim.max, by = binwidth)),table)))
  tryCatch(print(colSums(e)))
    }
}

dist.fc.range <- function(df, FCcolName = "FC", RangeColName = "range", convert.FC.to.logFC = T, prior.count = 1, convert.Range.to.logRange = T, as.absolute.distance = F, combine = F) {
    df <- df[c(FCcolName, RangeColName)]

    if (convert.FC.to.logFC == T) {
        if (min(df$FC) > 0) {
        df$log2FC <- log2(df$FC)
        } else {
        df$log2FC <- ifelse(df$FC >= 0, log2(df$FC+prior.count), -log2(abs(df$FC)+prior.count))
        }
    }

    if (convert.Range.to.logRange == T) {
        df$log2range <- ifelse(df$range >= 0, log2(df$range+prior.count), -log2(abs(df$range)+prior.count))
        }

    for (i in 1:nrow(df)) {
    df$euc.dist[i] <- dist(rbind(c(df$range[i], 0), c(0,df$FC[i])), method = "euclidian")
    }
    for (i in 1:nrow(df)) {
    df$man.dist[i] <- dist(rbind(c(df$range[i], 0), c(0,df$FC[i])), method = "manhattan")[1]
    }

    for (i in 1:nrow(df)) {
    df$euc.dist.log[i] <- dist(rbind(c(df$log2range[i], 0), c(0,df$log2FC[i])), method = "euclidian")
    }
    for (i in 1:nrow(df)) {
    df$man.dist.log[i] <- dist(rbind(c(df$log2range[i], 0), c(0,df$log2FC[i])), method = "manhattan")[1]
    }

    df$man.dist.rank <- rank(df$man.dist) / nrow(df)
    df$euc.dist.rank <- rank(df$euc.dist) / nrow(df)

    df$euc.dist.in.unit.circle <- df$euc.dist / max(df$euc.dist)
    df$man.dist.in.unit.circle <- df$man.dist / max(df$man.dist)

    if (as.absolute.distance == F) {
        for (r in 1:nrow(df)) {
            for (c in 5:ncol(df)) {
            df[r,c] <- ifelse(df$range[r] >= 0, df[r,c], -1*df[r,c])
            }
        }
    }

    if (combine == T) {
    df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
    df <- merge(df,df2, by = "row.names")
    row.names(df) <- df$Row.names; df$Row.names <- NULL
    }

    tryCatch(df <- df[order(df$man.dist.in.unit.circle), ])

    return(df)
}

fread.url.gz <- function(url, header=F, remove_empty_rows=F, ...) {

    
if(grepl(".gz", url)) test <- fread(paste(paste('URL=',url, sep = ""), ' ; curl "$URL" | gunzip -c', sep = ""), ...)
if(!grepl(".gz", url)) test <- fread(paste(paste('URL=',url, sep = ""), ' ; curl "$URL"', sep = ""), ...)

if(remove_empty_rows == T) {
    counter.frame = data.frame(counter=1:50)
    for (i in 1:50) {
        counter=1
            for (k in 1:ncol(test)) {
                if(test[i,k, with=F] != "") counter <- counter+1
            }
        counter.frame$counter[i] = counter
    }
}

    if(header == T) {
        test <- test[c(min(which(counter.frame$counter == median(counter.frame$counter))):nrow(test)), ]
        colnames(test) <- as.character(test[1, ])
        test <- test[2:nrow(test), ]
    }

    return(test)
}

#' @title Cast multiple columns
#' @param df which data.frame to cast
#' @param cast.formula formula, e.g. gene ~ sample, to be used in casting
#' @param cols the names of the columns to cast, e.g. cols = c("colname1", "colname2")
#' @export
#' @import 
dcast.multiple <- function(df, cast.formula, cols) {
    cast.formula = formula(cast.formula)
    
    for (i in 1:length(cols)) {
        df.temp <- dcast(df, cast.formula, value.var = cols[i])
        colnames(df.temp) <- c(colnames(df.temp)[1], paste(colnames(df.temp)[2:ncol(df.temp)], cols[i], sep = "."))
        
        nCols = ncol(df.temp) - 1
        k <- list(NULL)
        
        ## remove columns that have all values identical
        for (j in 2:nCols) {
            for (jj in 3:nCols) {
                if( "TRUE" %in% all.equal(df.temp[j], df.temp[jj]) & j != jj )  k[length(k)+1] <- j
            }
        }
        
        k <- unlist(k)
        df.temp <- df.temp[colnames(df.temp)[!colnames(df.temp) %in% colnames(df.temp)[k] ]]
        
        if(i == 1) df2 <- df.temp
        
        df.temp <- as.data.table(df.temp)
        df2 <- as.data.table(df2)
        
        if (i > 1 ) df2 <- merge(df2, df.temp, by = colnames(df2)[1])
    }
    
    return(df2)
}

read.quant.sf <- function(dir, subset.files=NULL,rename.samples=NULL,reorder.samples=T,as.wide.format=F, pretty.print=F) {
    files <- paste(file.path(dir), list.files(dir), "", sep = "/")
    files = paste(files, "quant.sf", sep = "")
    names(files) <- gsub(paste(dir, "quant.sf", "/", "_quant", sep="|"), "", files)
    
    
    if(is.null(subset.files)) {
        cat("You can subset files by specifying folder names in subset.files:\n")
        cat(names(files))
    } else {
        files <- files[names(files) %in% subset.files]
    }
    
    if(!is.null(rename.samples)) names(files) <- rename.samples
    if(reorder.samples==T) files <- files[order(names(files))]
    
    
    if(as.wide.format==F) {
        
        for(i in 1:length(files)) {
            temp <- fread(files[i])
            temp$Sample <- names(files[i])
            if(i == 1) temp.all <- temp
            if(i != 1) temp.all <- rbind(temp.all, temp)
        }
        
    } else {
    
        for(i in 1:length(files)) {
            
            temp <- fread(files[i])
            if(i==1) initial.colnames=colnames(temp)[3:ncol(temp)]
            
            colnames(temp) <- c("Name", "Length", paste(colnames(temp)[3:ncol(temp)], names(files[i]), sep = "_"))
            
            if(i == 1) {
                temp.all <- temp
                names <- temp[, 1, with = F]
            } else {
                if (all.equal(names, temp[, 1, with = F])) { 
                    temp.all <- cbind(temp.all, temp[, 3:ncol(temp), with = F])
                } else {
                    print("Exiting, Name column not equal in all files")
                    break()
                }
            }
        }
        
        if(pretty.print == T) {
            colnams = colnames(temp.all)[3:ncol(temp.all)]
            
            colum.reorder = colnams[ grep(initial.colnames[1], colnams, fixed = T) ]
            for (ii in 2:length(initial.colnames)) {
                colum.reorder = c(colum.reorder, colnams[ grep(initial.colnames[ii], colnams, fixed = T) ])
            }
            
            temp.all <- temp.all[, c("Name", "Length", colum.reorder[!duplicated(colum.reorder)]), with = F]
        }
    }
    
    return(temp.all)

}

sshfs <- function(remote_path, local_path, servername=NULL, id_rsa_path = NULL, useStallo = T) {
if(is.null(servername)) servername = "put001@stallo.uit.no:"
remote_path=paste(servername, remote_path, sep = "")
suppressWarnings(if(file.exists(local_path)) system(paste("umount -f", local_path)))
system(paste( "mkdir -p", local_path))
if(is.null(id_rsa_path)) id_rsa_path="~/.ssh/id_rsa"
system(paste("sshfs -o Ciphers=arcfour,Compression=no,auto_cache,reconnect,allow_other,defer_permissions,IdentityFile=", id_rsa_path, " ", remote_path, " ", local_path, sep = ""))
}


fread.sshfs <- function(file_with_remote_path, servername=NULL, id_rsa_path = NULL, force = F, ...) {
if(is.null(servername)) servername = "put001@stallo.uit.no:" # put your most used servername here to avoid typing in servername everytime instantiating this process
remote_path=paste(servername, file_with_remote_path, sep = "")
HOME=system("printf $HOME", intern = T)

locally_mounted_file_path = paste(HOME, "/share/", gsub(".*/", "", dirname(file_with_remote_path)), sep = "")
if(file.exists(locally_mounted_file_path) & force == F) cat("\nMount point already exists. Use another mountpath or set FORCE = TRUE")
if(file.exists(locally_mounted_file_path) & force == F) break() 
suppressWarnings(dir.create(locally_mounted_file_path, recursive = T))
if(is.null(id_rsa_path)) id_rsa_path="~/.ssh/id_rsa" # put path to id_rsa file here to avoid typing in path everytime instantiating this process
cat(paste("Mounting", dirname(remote_path), "to", locally_mounted_file_path, "..."))
system(paste("sshfs -o Ciphers=arcfour,Compression=no,auto_cache,reconnect,allow_other,defer_permissions,IdentityFile=", id_rsa_path, " ", dirname(remote_path), " ", locally_mounted_file_path, sep = ""))

file = paste(locally_mounted_file_path, gsub(dirname(file_with_remote_path), "", file_with_remote_path), sep = "")
cat("\nReading", file)
library(data.table)
temp <- fread(file, ...)

cat(paste("Umounting and removing", locally_mounted_file_path), "...")
system(paste("umount -f", locally_mounted_file_path))
unlink(locally_mounted_file_path, recursive = T, force = T)

return(temp)
}

source_https <- function(url, ...) {
  # load package
  require(RCurl)
 
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}

tdt <- function(DT, transpose.col, ...) {
# The transpose function is efficient, but lacks the keeping of row and colnames
new.row.names <- colnames(DT)
new.row.names <- new.row.names[!new.row.names %in% transpose.col]
new.col.names <- DT[, transpose.col, with = F]
DT <- DT[, !colnames(DT) %in% transpose.col, with = F]
DT <- transpose(DT, ...)
colnames(DT) <- unlist(new.col.names)
DT$var <- new.row.names
# change order of DT after transposing 
setcolorder(DT, c("var", setdiff(names(DT), "var")))
colnames(DT)[1] <- transpose.col
return(DT)
}

# input df is a data.table with gene in rows and samples in col. first col has to contain name of gene. var.list is extracted from first column
# time.event.data is data.frame with cols binary data (event) in col1 and time in col2. The row order has to match col order in df
# sub.idx has to be subset containg index on samples to subset on (e.g. 1,3,5). Example which(type == "highrisk")
coxph.parallell <- function(df, time.event.data, var.list = NULL, file = "test.txt", sub.idx = NULL, read.file = F, no_cores = "") {
    start = Sys.time()
    
    if (!is.null(sub.idx)) {
        df <- df[, c(1, sub.idx+1), with = F]
        time.event.data <- time.event.data[sub.idx, ]
    }
    
    colnames(time.event.data) <- c("event", "time")
    colnames(df)[1] <- "var"
    setkey(df, var)
    
    if (is.null(var.list)) {
        var.list <- df$var
    }

    if (no_cores == "") {
        no_cores <- detectCores() - 1
        registerDoParallel(no_cores)
    }
    
    # Creates and overwrites file if it already exists
    dir.create(dirname(file), showWarnings = F)
    file.create(file, overwrite = T)
    
    nr_times = floor(nrow(df) / no_cores) - 1
    remainder = nrow(df) %% no_cores
    
    ## Clean-up any existing temp files
     for (i in 1:(no_cores+1)) {
        tryCatch(unlink(paste(file, i, "temp.txt", sep = "_"), force = T))
    }
    
    for (i in 0:nr_times) {
        foreach (z = 1:no_cores,
                 .combine = rbind)  %dopar%
                 {
            name = var.list[((i * no_cores) + z)] # contains name
            t1 <- df[name, 2:ncol(df)] # contains only binary or continous values
            
            case.table = table(unlist(t1)) # how many yes no. This is removed later if var is continous
            
            survdata <- cbind(time.event.data, t(t1))
            colnames(survdata) <- c("event", "time", "var")
            
            res.cox <- coxph(Surv(time, event) ~ var, data = survdata, na.action = na.omit)
            chi.sq.p <- cox.zph(res.cox) # test the proportional hazards assumption of cox regression.
            
            s <- summary(res.cox)
            
            res <- list(var=name, 
                n=sum(case.table), 
                yes=round(case.table[2]/sum(case.table),2)*100, 
                no=round(case.table[1]/sum(case.table),2)*100, 
                estimate = round(s$coefficients[1,1],2), 
                se.coef=round(s$coefficients[1,3],2), 
                hazard.ratio=round(s$coefficients[1,2],2), 
                zvalue=round(s$coefficients[1,4],2), 
                p.wald=signif(s$coefficients[1,5], 3),  
                p.lrt = signif(s$logtest[3], 3),
                p.logrank = signif(s$sctest[3], 3),
                rsq=round(s$rsq[1], 3),
                concordance = round(s$concordance[1], 2),
                concordance.se = round(s$concordance[2], 2),
                CI.lower = round(s$conf.int[3],2), 
                CI.upper = round(s$conf.int[4],2), 
                chi.sq.p = signif(chi.sq.p$table[3], 3)
                )
            
            # if var is continous we removed the number of yes's and no's
            if (length(table(unlist(t1))) != 2) res$yes <- NULL
            if (length(table(unlist(t1))) != 2) res$no <- NULL
            
            cat(unlist(res), sep ="\t", file = paste(file, z, "temp.txt", sep = "_"), append = T) # add tab and newline, then append to file
            cat("", sep ="\n", file = paste(file, z, "temp.txt", sep = "_"), append = T) # add tab and newline, then append to file
                 }
            }
    
    
        if (remainder != 0) {
                z <- length(grep(paste(file, ".", "temp.txt", sep = "_"), list.files(dirname(file), full.names = T))) + 1
            
                file.rem = paste(file, z, "temp.txt", sep = "_")
            
                # add in the rest
                for (i in 1:remainder) {
                
            name = var.list[((nr_times + 1) * no_cores) + i] # contains name
            t1 <- df[name, 2:ncol(df)] # contains only binary 
            case.table = table(unlist(t1)) # how many yes no 
            survdata <- cbind(time.event.data, t(t1))
            colnames(survdata) <- c("event", "time", "var")
            
            res.cox <- coxph(Surv(time, event) ~ var, data = survdata, na.action = na.omit)
    
            chi.sq.p <- cox.zph(res.cox) # test the proportional hazards assumption of cox regression.
            
            s <- summary(res.cox)
            
            res <- list(var=name, 
                n=sum(case.table), 
                yes=round(case.table[2]/sum(case.table),2)*100, 
                no=round(case.table[1]/sum(case.table),2)*100, 
                estimate = round(s$coefficients[1,1],2), 
                se.coef=round(s$coefficients[1,3],2), 
                hazard.ratio=round(s$coefficients[1,2],2), 
                zvalue=round(s$coefficients[1,4],2), 
                p.wald=signif(s$coefficients[1,5], 3),  
                p.lrt = signif(s$logtest[3], 3),
                p.logrank = signif(s$sctest[3], 3),
                rsq=round(s$rsq[1], 3),
                concordance = round(s$concordance[1], 2),
                concordance.se = round(s$concordance[2], 2),
                CI.lower = round(s$conf.int[3],2), 
                CI.upper = round(s$conf.int[4],2),
                chi.sq.p = signif(chi.sq.p$table[3], 3)
                )
            
            # if var is continous we removed the number of yes's and no's
            if (length(table(unlist(t1))) != 2) res$yes <- NULL
            if (length(table(unlist(t1))) != 2) res$no <- NULL
            
            cat(unlist(res), sep ="\t", file = file.rem, append = T) # add tab and newline, then append to file
            cat("", sep ="\n", file = file.rem, append = T) # add tab and newline, then append to file
                
            }
        } 
    
    stopImplicitCluster()
    
    # file.cat <- gsub(paste(z, "temp.txt", sep = "_"), "*_temp.txt", file.rem)
    file.cat <- paste(file, '*', "temp.txt", sep = "_")
    
    cat("\nNumber of correlations in temp files:\n")
    tryCatch(system(paste("wc -l", file.cat)))
    
    cat(paste("\nNumber of variables provided to this function:", length(var.list)))
    
    ## combine all files 
    system(paste("cat", file.cat, ">", file))
    ## clean up temp files
    unlink(file.cat, force = T)

    cat(paste("\nNumber of correlations in final file:\n", system(paste("wc -l", file), intern = T)))
    
    end = Sys.time()
    time = round((end - start), 2)
    cat(paste("\nFinished in", time, units(time)))

}

plot.surv <- function(genename, clindata, counts, sub.idx = NULL, textsize = 18, count.index = NA, log.transform = T, cutoff.modus = c("mean", "scan", "median", "Kmeans", "user.defined"),scan.nr = 5, user.defined = NULL, max.time = NA, min.time = NA, ...) {
        if(log.transform == T) cat("Log transforming data with pseudocount 1\n")
        if(log.transform == F) cat("Using raw data. Not log transforming\n")
         
        if(length(cutoff.modus) == 5) cutoff.modus = "scan"
        if(!is.null(user.defined)) cutoff.modus = "user.defined"
        cat(paste("Using", cutoff.modus, "as cutoff.modus\n"))
        cat("Cutoff modus may be mean, median, Kmeans, scan or user.defined\n")
        
        cat(paste("Testing", genename, "\n"))
        colnames(counts)[1] <- "ensembl_gene_id"
        
        if (is.na(count.index)) {
            survdata <- cbind(clindata, t(counts[counts$ensembl_gene_id == hugo2ensg(genename)[1,2], 2:ncol(counts)])) 
        } else {
            survdata <- cbind(clindata, t(counts[count.index, 2:ncol(counts)]))
        }
        
        survdata <- survdata[, 1:3]
        colnames(survdata) <- c("sbin", "time", "var")

         if (cutoff.modus == "user.defined") {
                survdata$bin <- as.factor(user.defined)
         } 
        
        
        if (!is.null(sub.idx)) {
            survdata <- survdata[sub.idx, ]    
            survdata$bin <- factor(as.character(survdata$bin)) # redo this to remove any levels that should not be present
        }
        
        if (!is.na(max.time)) survdata <- survdata[survdata$time <= max.time, ]
        if (!is.na(min.time)) survdata <- survdata[survdata$time >= min.time, ]
        
        
        # order survdata in increasing expression
        survdata <- survdata[order(survdata$var), ]
        
        if(log.transform == T) survdata$var <- log2(survdata$var+1)
        
            # categorize data based on of the function scan, kmeans, median or mean
            if(cutoff.modus == "scan") {
                survdata$bin = 0
                survdata$bin[1:scan.nr] <- 1 # set 1 for the first lowly expressed samples
        
                pval <- list(pval = rep(1,length(survdata$bin)))
                raw.pval <- list(pval = rep(1,length(survdata$bin)))
                
                for(i in 1:nrow(survdata)) {
                    res.cox <- suppressWarnings(coxph(Surv(time, sbin) ~ bin, data = survdata, na.action = na.omit))
                    
                    s <- summary(res.cox)
                    raw.pval$pval[i] <- ifelse(s$logtest[3] == 0, 2e-16, s$logtest[3])
                    pval$pval[i] <- p.adjust(raw.pval$pval[i], method = "bonferroni", n = i)
                    survdata$bin[i] <- 1
                }
                
                # we would not like to select among first and last samples 
                choose.from = (scan.nr+1):(length(raw.pval$pval)-scan.nr-1)
                pvalues = which(pval$pval[choose.from] == min(pval$pval[choose.from]))
                # if all pvalues were adjusted to 1, we choose one of from raw pvalues
                use.raw.pvalues=F
                if(length(pvalues) > 1) use.raw.pvalues = T
                if(use.raw.pvalues == T) pvalues = which(raw.pval$pval[choose.from] == min(raw.pval$pval[choose.from]))
                pvalues <- pvalues+(scan.nr)
                
                survdata$bin <- 2
                survdata$bin[1:(pvalues[length(pvalues)]-1)] <- 1
                col =  c(rep("grey", table(survdata$bin)[1]), rep("black", nrow(survdata) - table(survdata$bin)[1]))
                
                if(use.raw.pvalues == F) df.bar <- barplot(-log10(pval$pval), ylab = "-log10 PValue", xlab = "sample", col = col, border = col)
                if(use.raw.pvalues == T) df.bar <- barplot(-log10(raw.pval$pval), ylab = "-log10 PValue", xlab = "sample", col = col, border = col)
                # plot points of raw values (not log transformed)
                if(log.transform == T) point.exps = 2^(survdata$var-1)
                if(log.transform == F) point.exps = survdata$var
                if(use.raw.pvalues == F) point.exps <- point.exps / (max(point.exps) / max(-log10(pval$pval)))
                if(use.raw.pvalues == T) point.exps <- point.exps / (max(point.exps) / max(-log10(raw.pval$pval)))
                points(df.bar, point.exps, col = "orange", pch = 15)
                title(paste("Scan: logrank iterative test\n", "P=", signif(min(raw.pval$pval),2),"P.bonf=", signif(min(pval$pval),2)))
            }
            
            if (cutoff.modus == "Kmeans") {
                k.fit <- kmeans(survdata[,3], 2, iter.max = 100, nstart = 50) 
                # get cluster means
                agg <- aggregate(survdata[,3],by=list(k.fit$cluster),FUN=mean)
                # append cluster assignment
                cl <- data.frame(survdata[,3], k.fit$cluster)
                cl <- cl[order(cl$k.fit.cluster),]
                cl <- cl[order(cl[,1]), ]
                survdata <- cbind(survdata, cl[, 2])
                colnames(survdata)[4] <- "bin"
            }
            
            if (cutoff.modus == "median") {
                if (log.transform == T) median.sample <- which(round(survdata[,3], 1) == median( unlist(round(survdata[,3], 1))))
                if (log.transform == F) median.sample <- which(round(survdata[,3], 0) == median( unlist(round(survdata[,3], 0))))
                
                survdata$bin <- 2
                survdata$bin[1:length(median.sample)] <- 1
            }
            
            if (cutoff.modus == "mean") {
                mean.sample <- which(survdata[,3] <= mean(unlist(survdata[,3])))
                survdata$bin <- 2
                survdata$bin[1:length(mean.sample)] <- 1
            }
        
                    if (cutoff.modus != "user.defined") {
                        # Categorize into LOW or HIGH
                        if(mean(unlist(survdata[survdata$bin == 1, "var"])) > mean(unlist(survdata[survdata$bin == 2, "var"]))) {
                            survdata$bin <- gsub("1", "HIGH", survdata$bin) 
                            survdata$bin <- gsub("2", "LOW", survdata$bin)
                        } else {
                            survdata$bin <- gsub("1", "LOW", survdata$bin) 
                            survdata$bin <- gsub("2", "HIGH", survdata$bin)
                        }
                    
                        survdata$bin <- factor(survdata$bin, levels = c("LOW", "HIGH"))
                    }
        
        print(table(survdata$bin))
        assign("survdata", survdata, envir = .GlobalEnv)
        
        fit <- survfit(Surv(time, sbin) ~ bin, data = survdata, na.action = na.omit)
        assign("survfit", fit, envir = .GlobalEnv)
        print(coxph(Surv(time, sbin) ~ bin, data = survdata, na.action = na.omit))
        
        survdiff <- survdiff(Surv(time, sbin) ~ bin, data = survdata, na.action = na.omit)
        assign("survdiff", survdiff, envir = .GlobalEnv)
        print(survdiff(Surv(time, sbin) ~ bin, data = survdata, na.action = na.omit))
        
        case.table = table(survdata$bin)
        
        res.cox <- coxph(Surv(time, sbin) ~ bin, data = survdata, na.action = na.omit)
        chi.sq.p <- cox.zph(res.cox) # test the proportional hazards assumption of cox regression
                
        s <- summary(res.cox)
                
                res <- list(var=genename, 
                    n=sum(case.table), 
                    yes.percent=round(case.table[2]/sum(case.table),2)*100, 
                    no.percent=round(case.table[1]/sum(case.table),2)*100, 
                    estimate = round(s$coefficients[1,1],2), 
                    se.coef=round(s$coefficients[1,3],2), 
                    hazard.ratio=round(s$coefficients[1,2],2), 
                    zvalue=round(s$coefficients[1,4],2), 
                    p.wald=signif(s$coefficients[1,5], 3),  
                    p.lrt = signif(s$logtest[3], 3),
                    p.logrank = signif(s$sctest[3], 3),
                    rsq=round(s$rsq[1], 3),
                    concordance = round(s$concordance[1], 2),
                    concordance.se = round(s$concordance[2], 2),
                    CI.lower = round(s$conf.int[3],2), 
                    CI.upper = round(s$conf.int[4],2),
                    chi.sq.p = signif(chi.sq.p$table[3], 3))
        
        assign("surv.coxph", res, envir = .GlobalEnv)  
        assign("surv.coxph.summary", s, envir = .GlobalEnv)  
        return(ggsurvplot(fit, ...))
}

getSeriesMatrixCharacteristics <- function(ftpURLof_series_matrix.txt.gz) {
    file = gsub(".*/", "", ftpURLof_series_matrix.txt.gz)
    file = paste(getwd(), file, sep = .Platform$file.sep)
    system(paste("wget ", ftpURLof_series_matrix.txt.gz, " -O ", file, sep = ""))
    system(paste("gzip -f -d ", file, sep = ""))
    file = gsub(".gz", "", file)
    temp.file = readLines(file)
    headers <- read.delim2(file, header = F, skip = grep('!Sample_title', temp.file, fixed = T)-1, nrows = 1)
    headers <- unlist(headers, recursive = F, use.names = F)
    headers <- gsub("\\ .*", "", headers)
    skip=grep('!Sample_characteristics_ch1', temp.file)
    series_matrix <- read.delim2(file, header = F, skip = skip[1]-1, nrows = length(skip), col.names = headers)
    colnames(series_matrix) <- c("characteristic", colnames(SEQC)[2:ncol(SEQC)])
    series_matrix$characteristic <- gsub("\\:.*","", series_matrix$SEQC_NB001)
    n <- colnames(series_matrix)
    series_matrix <- as.data.table(t(series_matrix))
    series_matrix <- cbind(n, series_matrix)
    
    colnames(series_matrix) <- unlist(series_matrix[1])
    series_matrix <- series_matrix[2:nrow(series_matrix), ]
    series_matrix <- as.data.frame(series_matrix)
    row.names(series_matrix) <- series_matrix$characteristic
    series_matrix$characteristic <- NULL
    
    for (i in 1:nrow(series_matrix)) {
        temp <- unlist(series_matrix[i, ], use.names = F)[match(colnames(series_matrix), gsub("\\:.*","", unlist(series_matrix[i, ], use.names = F)))]
        if(i == 1) temp.all <- temp
        if(i != 1) temp.all <- rbind(temp.all, temp)
    }
    
    row.names(temp.all) <- row.names(series_matrix)
    colnames(temp.all) <- colnames(series_matrix)
    series_matrix <- temp.all
    series_matrix <- apply(series_matrix, 2, function(x) gsub(".*\\:","", x))
    series_matrix <- cbind(headers[2:length(headers)], series_matrix)
    colnames(series_matrix)[1] <- "sample"
    series_matrix <- as.data.table(series_matrix)
    colnames(series_matrix) <- gsub(" ", "_", colnames(series_matrix))
    series_matrix <- apply(series_matrix, 2, function(x) gsub(" ", "", x))
    series_matrix <- as.data.frame(series_matrix, stringAsFactors = F)

    return(series_matrix)
}
