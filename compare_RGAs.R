# compares RGAs with PAV and SNPs
library(dplyr)
library(ggplot2)
library(patchwork)


df <- read.table('../Final_Genes_Table.csv', header = TRUE, sep='\t')

df <- df %>% filter(X.GeneName != 'NO_SNPS')
# ignore TM-CCs
df <- df %>% filter(Class != 'TM-CC')
df <- df %>% filter(PAV != 'NA')
# First plot, compare low, moderate, high imapct SNPs per bp with RGA class and PAV
ap1 <- ggplot(df) + geom_boxplot(aes(x=Class, y=variants_impact_LOW/Total_exon_length, fill=PAV)) + ylab('Low impact/bp') + theme_bw() 
ap2 <- ggplot(df) + geom_boxplot(aes(x=Class, y=variants_impact_MODERATE/Total_exon_length, fill=PAV)) + ylab('Moderate impact/bp') + theme_bw()
ap3 <- ggplot(df) + geom_boxplot(aes(x=Class, y=variants_impact_HIGH/Total_exon_length, fill=PAV)) + ylab('High impact/bp') + theme_bw()

ap1 + ap2 + ap3 + plot_layout(ncol=1)

ggsave('Variants_exon_length_vs_class.png', plot=ap1 + ap2 + ap3 + plot_layout(ncol=1), dpi=1200)

# Second plot, compare upstream/downstream SNPs with RGA/PAV
bp1 <- ggplot(df) + geom_boxplot(aes(x=Class, y=variants_effect_downstream_gene_variant, fill=PAV)) + ylab('Downstream variants') + theme_bw()
bp2 <- ggplot(df) + geom_boxplot(aes(x=Class, y=variants_effect_upstream_gene_variant, fill=PAV)) + ylab('Upstream variants') + theme_bw()

bp1 + bp2 + plot_layout(ncol=1)
ggsave('Variants_up_down_stream_vs_class.png', plot=bp1 + bp2 + plot_layout(ncol=1), dpi=1200)


# Run MANOVA and two-way ANOVA on the same data
summary(manova(cbind(variants_impact_LOW/Total_exon_length, variants_impact_MODERATE/Total_exon_length, variants_impact_HIGH/Total_exon_length) ~ Class, data=df))
summary(aov(variants_impact_LOW/Total_exon_length ~ Class+PAV, data=df))
summary(aov(variants_impact_MODERATE/Total_exon_length ~ Class+PAV, data=df))
summary(aov(variants_impact_HIGH/Total_exon_length ~ Class+PAV, data=df))

# now lets get rid of outliers
# methods from http://rpubs.com/lcollado/7904 
# removes outliers past the third standard deviation (sets to NA)

findOutlier <- function(data, cutoff = 3) {
  ## Calculate the sd
  sds <- apply(data, 2, sd, na.rm = TRUE)
  ## Identify the cells with value greater than cutoff * sd (column wise)
  result <- mapply(function(d, s) {
    which(d > cutoff * s)
  }, data, sds)
  result
}
removeOutlier <- function(data, outliers) {
  result <- mapply(function(d, o) {
    res <- d
    res[o] <- NA
    return(res)
  }, data, outliers)
  return(as.data.frame(result))
}

df$low <- df$variants_impact_LOW/df$Total_exon_length
df$mod <- df$variants_impact_MODERATE/df$Total_exon_length
df$high <- df$variants_impact_HIGH/df$Total_exon_length
outliers <- findOutlier(df)
dfFilt <- removeOutlier(df, outliers)
dfFilt$Class <- df$Class
dfFilt$PAV <- df$PAV
dfFilt$Cluster <- df$Cluster

# same as above, supplementary figure
ap1 <- ggplot(dfFilt) + geom_boxplot(aes(x=Class, y=low, fill=PAV)) + ylab('Low impact/bp') + theme_bw() 
ap2 <- ggplot(dfFilt) + geom_boxplot(aes(x=Class, y=mod, fill=PAV)) + ylab('Moderate impact/bp') + theme_bw()
ap3 <- ggplot(dfFilt) + geom_boxplot(aes(x=Class, y=high, fill=PAV)) + ylab('High impact/bp') + theme_bw()
ap1 + ap2 + ap3 + plot_layout(ncol=1)
ggsave('Variants_exon_lengths_without_outliers.png', plot=ap1 + ap2 + ap3 + plot_layout(ncol=1), dpi=1200)

summary(manova(cbind(low, mod, high) ~ Class, data=dfFilt))

# same as above, ANOVA
summary(aov(low ~ PAV+Class, data=dfFilt))
summary(aov(mod ~ PAV+Class, data=dfFilt))
summary(aov(high ~ PAV+Class, data=dfFilt))

# now check for physical clusters
cp1 <- ggplot(df) + geom_boxplot(aes(x=Class, y=low, fill=Cluster)) + ylab('Low impact/bp') + theme_bw() 
cp2 <- ggplot(df) + geom_boxplot(aes(x=Class, y=mod, fill=Cluster)) + ylab('Moderate impact/bp') + theme_bw()
cp3 <- ggplot(df) + geom_boxplot(aes(x=Class, y=high, fill=Cluster)) + ylab('High impact/bp') + theme_bw()
cp1 + cp2 + cp3 + plot_layout(ncol=1)

summary(aov(low ~ Cluster, data=df))
summary(aov(mod ~ Cluster, data=df))
summary(aov(high ~ Cluster, data=df))


# time for Fisher's exact test!
# core, variable, low, moderate
core <- filter(df, PAV=='Core')
variable <- filter(df, PAV=='Variable')
var_low <- sum(variable$variants_impact_LOW)
var_med <- sum(variable$variants_impact_MODERATE)
var_high <- sum(variable$variants_impact_HIGH)
core_low <- sum(core$variants_impact_LOW)
core_med <- sum(core$variants_impact_MODERATE)
core_high <- sum(core$variants_impact_HIGH)
variable_total_exon <- sum(variable$Total_exon_length)
core_total_exon <- sum(core$Total_exon_length)

mytable <-  data.frame(low=c(core_low, var_low), high=c(core_high, var_high))
fisher.test(mytable)
mytable <-  data.frame(low=c(core_low, var_low), high=c(core_med, var_med))
fisher.test(mytable)


# genes that have high impact SNPs/have 0 SNPs vs. genes that are core/variable

high_impact <- filter(df, variants_impact_HIGH == 0)
no_high_impact <- filter(df, variants_impact_HIGH != 0)
nrow(filter(high_impact, PAV=='Core'))
core_high <- nrow(filter(high_impact, PAV=='Core'))
variable_high <- nrow(filter(high_impact, PAV=='Variable'))
core_no_high <- nrow(filter(no_high_impact, PAV=='Core'))
variable_no_high <- nrow(filter(no_high_impact, PAV=='Variable'))

mytable <-  data.frame(high_impact=c(core_high, variable_high), no_high=c(core_no_high, variable_no_high))
mytable
fisher.test(mytable)


# genes that have low impact SNPs/have 0 SNPs vs. genes that are core/variable

low_impact <- filter(df, variants_impact_LOW == 0)
no_low_impact <- filter(df, variants_impact_LOW != 0)

core_low <- nrow(filter(low_impact, PAV=='Core'))
variable_low <- nrow(filter(low_impact, PAV=='Variable'))
core_no_low <- nrow(filter(no_low_impact, PAV=='Core'))
variable_no_low <- nrow(filter(no_low_impact, PAV=='Variable'))

mytable <-  data.frame(low_impact=c(core_low, variable_low), no_low=c(core_no_low, variable_no_low))
mytable
fisher.test(mytable)
