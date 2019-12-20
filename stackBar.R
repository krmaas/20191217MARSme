install.packages("reshape2")
install.packages("tidyverse")
install.packages("plyr")

library(reshape2)
library(tidyverse)
library(ggplot2)
library(plyr)


taxa.bars <- read.delim(file="../../mothurProcess/dec17ws.taxa.barplot.txt", 
                        header=T)

taxa.col <- c("Frankiales"="#d94801",
    "ActinobacteriaOther"="#f16913",
    "Betaproteobacteria"="#79C360",
    "Lactobacillales"="#A6CEE3",
    "Lachnospiraceae"="#7DB4D5",
    "Flavobacterium"="#08519c",
    "BacteriodetesOther"="#5C9FC9",
    "Ruminococcaceae"="#3A89BD",
    "Coriobacteriia"="#FB9A99",
    "RF9"="#fdae6b",
    "Bacteria_Other"="#d3d3d3",
    "Bacteria_unclassified"="#000000")

taxa.bars1 <- left_join(alpha.env, taxa.bars, by="group")
taxa.bars <- select(taxa.bars1, group, 
                    Type, Site, Frankiales, ActinobacteriaOther, 
                    Betaproteobacteria, Flavobacterium, BacteriodetesOther,
                    Bacteria_Other)

tax <- melt(taxa.bars, id.vars = c("group", "Type", "Site"))

tax2 <- aggregate(value~variable+Type, data=tax, FUN=sum)




ggplot(tax2,  aes( y=value, x=factor(Type), color=NULL, fill=factor(variable), order=-as.numeric(variable)))+ #order=- makes the stack the same order as the legend
  geom_bar(position="fill", stat="identity")+
    xlab("type")+
    ylab("Percent total community")+
    scale_fill_manual(values=taxa.col, guide=guide_legend(title="Taxonomic level"))+
    theme_bw()
ggsave(file="barplot.tiff", height=5, width=5, dpi=600)


##### anova on each taxa

anova.out <- dlply(tax, .(variable), function(tax) aov(value~Type, data=tax))
juicy.bits <- function(x)
{c(anova(x)[1,4], anova(x)[1,5])} # from http://oardc.osu.edu/culman/wp-content/uploads/2014/05/Reshape-and-Plyr.txt

#V1 is intercept, V2 p-value
ldply(anova.out, juicy.bits)

tukey.out <- dlply(tax, .(variable), function(tax) TukeyHSD(aov(value~Type, data=tax)))
tuk.juice <- function(x)
{(x)$Type[,c(1,4)]}
tuk1 <- ldply(tukey.out, tuk.juice)
write.csv(file="tukey.taxa.csv", tuk1)


