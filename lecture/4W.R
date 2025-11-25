setwd('Desktop/Study/Bioinformatics_practice/lecture/data/')

exdata <- read.table('GeneExpExample1000.txt',header = TRUE)
View(exdata)

cal_rpkm <- function(cnts, lens) {
  return(10^9 * (cnts/sum(cnts) * lens))
}

cal_tpm <- function(cnts, lens) {
  tpm <- 10^6*((cnts/lens)/sum(cnts/lens))
  return(tpm)
}

expr <- exdata[,7:ncol(exdata)]
rownames(expr) <- exdata$EnsemblGeneID

gene.len <- as.numeric(exdata[,4] - exdata[,3])

rpkm.data <- apply(expr, 2, function(x) cal_rpkm(x, gene.len))
tpm.data <- apply(expr,2, function(x) cal_tpm(x, gene.len))

apply(rpkm.data, 2, sum)
apply(tpm.data, 2, sum)

rpkm.data.lg2 <- log2(rpkm.data + 1)
tpm.data.lg2 <- log2(tpm.data + 1)

hist(tpm.data[,1])
hist(tpm.data[,1], xlab = 'TPM', main = sprintf("Histogram of %s", colnames(tpm.data[1])))
hist(tpm.data.lg2[,1], xlab = 'log2(TPM+1)', main = sprintf("Histogram of %s", colnames(tpm.data.lg2[1])))


group <- rep('kidney', ncol(tpm.data.lg2))
group[grep("Liver", colnames(tpm.data.lg2))] == 'liver'
View(tpm.data.lg2)

group <- ifelse(grepl("Liver", colnames(tpm.data.lg2)), 'liver', 'kidney')

tres <-  apply(tpm.data.lg2,1,function(x) t.test(x~group))
tp <- sapply(tres, function(x) x$p.value)
td <- sapply(tres, function(x) x$estimate[2]-x$estimate[1])

names(td) <- names(tres)
logtp <- -log10(tp)

up.gene <- names(td[tp < 0.05 & td > 1])
dn.gene <- names(td[tp < 0.05 & td < -1])

plot(x=td, y=logtp, xlab = 'fold change', ylab = '-logP', main = "Volcano plot of kidney vs. liver")
points(x=td[names(td) %in% up.gene], y=logtp[names(tp) %in% up.gene], col='red')
points(x=td[names(td) %in% dn.gene], y=logtp[names(tp) %in% dn.gene], col='blue')
abline(h = -log10(0.05), v = c(1, -1), col = 'red', lty = 2)
legend('topleft', legend = c('Up in liver', 'Up in kidney'), col = c('red', 'blue'), pch=1, cex=0.8)

library(ggplot2)
mat <- data.frame(exdata[,1:6], "logP" = logtp, 'FC'=td)
mat$label <-"Nodiff"
mat$label[mat$logP > -log10(0.05) & mat$FC > 1] <- 'Liver'
mat$label[mat$logP > -log10(0.05) & mat$FC < -1] <- 'Kidney'

ggplot(data=mat, aes(x=FC, y=logP , col=label)) + geom_point()

ggplot(data=mat, aes(x=FC, y=logP , col=label)) + geom_point() +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
  geom_vline(xintercept = c(-1,1), linetype="dashed", color = "red") +
  scale_color_manual(values=c("blue","red","grey")) +
  theme_bw()


order(..., na.last = TRUE, decreasing = FALSE, method = c("auto", "shell", "radix"))
mat_sorted <- mat[order(mat$FC, decreasing = T),]
mat_sorted <- mat_sorted[mat_sorted$Status=="KNOWN",]
top.gene <- mat_sorted$EnsemblGeneID[1]
val <- tpm.data.lg2[top.gene,]
boxplot(val~group)

group <- factor(group, levels = c('liver', 'kidney'))
boxplot(val~group, col=c("red"))
boxplot(val~group, col=c("blue","red"), xlab="organ", ylab="log2(TPM+1)")


expr <- read.table('GSE25297_exprsdata.txt', header = TRUE, row.names = 1)
pheno <- read.table('GSE25297_phenodata.txt', header = TRUE, row.names = 1)
colnames(expr) <- pheno[colnames(expr), 'Sample']

exc.samples <- c("CNSlymphoma_rep4", "CNSlymphoma_rep5", "CNSlymphoma_rep7", "NormalLymphNode_rep7")
expr.ft <- expr[,!(colnames(expr) %in% exc.samples)]
match(colnames(expr.ft),pheno$Sample)
group <- pheno$Type[match(colnames(expr.ft),pheno$Sample)]
tres <- apply(expr.ft, 1, function(x) t.test(x~group))

tp <- sapply(tres, function(x) x$p.value)
td <- sapply(tres, function(x) x$estimate[2] - x$estimate[1])
df <- data.frame('ID' = names(tres), "P" = tp, "FC" = td)
df.sorted <- df[order(df$FC, decreasing = T),]
top.gene <- df.sorted$ID[1]
val <- as.matrix(expr.ft)[top.gene,]

val.avg <- aggregate(val, by=list(group), mean)
val.avg 
val.sd <- aggregate(val, by=list(group), sd)
val.sd

df2 <- data.frame(val.avg, val.sd[,-1])
colnames(df2) <- c("Group","avg","sd")

ggplot(df2, aes(x=Group, y=avg, fill=Group)) + geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.2) +
  scale_fill_manual(values=c("grey","darkred")) +
  theme_bw()

library(dplyr)
df3 <- data.frame('Sample' = colnames(expr.ft), "Group" = group, "Expr" = val)
df4 <- df3 %>% 
  group_by(Group) %>% 
  summarize(avg=mean(Expr), sd =sd(Expr) )

library(ggpubr) 
df3$Group = factor(df3$Group, levels = c("Normal", "PCNSL"))
ggbarplot(df3, x="Group", y="Expr",
          add =c('mean_sd'),
          fill="Group",
          palette = c("grey","darkred"))

top.gene.list <- df.sorted$ID[1:5]
mat <- as.matrix(expr.ft)[top.gene.list,]
mat2 <- t(mat)

library(reshape2)
df5 <- melt(mat2)
df5 <- melt(mat2, varnames = c("Sample","probe"), value.name = "Expr")
df5$Group <- pheno$Type[match(df5$Sample, pheno$Sample)]

df6 <- df5 %>% group_by(probe, Group) %>%
  summarize(avg = mean(Expr), sd = sd(Expr))

ggplot(df6, aes(x=probe, y=avg, fill=Group)) + geom_bar(stat="identity")

ggplot(df6, aes(x=probe, y=avg, fill=Group)) + geom_bar(stat="identity", position=position_dodge())

ggplot(df6, aes(x=probe, y=avg, fill=Group)) + geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=0.2, position = position_dodge(width=0.9)) +
  scale_fill_manual(values=c("grey","darkred"))
                    
ggbarplot(df5, x = "probe", y = "Expr", fill = "Group",
          add = "mean_sd", palette = c("#00AFBB", "#E7B800"),
          position = position_dodge(width=0.9)) +
  stat_compare_means(aes(group = Group), label = "p.signif")

#Quiz
ggplot(df5, aes(x=probe, y=Expr, fill = Group)) + geom_boxplot()
ggboxplot(df5, x = "probe", y = "Expr", fill = "Group")

ggplot(df5, aes(x=Group, y=Expr, fill=Group)) + geom_boxplot() + facet_grid(~probe)

ggboxplot(df5, x = "Group", y = "Expr",
          color = "Group", palette = "npg", #nature publishing group
          add = "jitter",
          facet.by = "probe", nrow=1) +
  stat_compare_means(aes(group = Group), label = "p.format", method = "t.test")

ggboxplot(df5, x = "Group", y = "Expr",
          color = "Group", palette = "npg", #nature publishing group
          add = "jitter",
          facet.by = "probe", nrow=3) +
  stat_compare_means(aes(group = Group), label = "p.format", method = "t.test")
