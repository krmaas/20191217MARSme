#### data analysis for 20191219 MARSme workshop

#### load libraries ####

# install.packages("ggplot2")
# install.packages("vegan")
# install.packages("dplyr")
# install.packages("indicspecies")
# install.packages("tidyverse")
# install.packages("tibble")

library(ggplot2)
library(vegan)
library(dplyr)
library(indicspecies)
library(tidyverse)
library(tibble)

parseDistanceDF = function(phylip_file) {
    
    # Read the first line of the phylip file to find out how many sequences/samples it contains
    temp_connection = file(phylip_file, 'r')
    len = readLines(temp_connection, n=1)
    len = as.numeric(len)
    len = len +1
    close(temp_connection)
    
    
    phylip_data = read.table(phylip_file, fill=T, row.names=1, skip=1, col.names=1:len)
    colnames(phylip_data) <- row.names(phylip_data)
    return(phylip_data)
}



#### read in data ####

otu <- read.table(file="../MARSme20191218.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.0.03.subsample.shared",
                  header=T, row.names=2, stringsAsFactors = F)

str(otu)

## remove label columns from otu

otu <- select(otu, -label, -numOtus)

# come back to reading in if we have time
# taxa <- read.table(file="../M")

alph <- read.table(file="../MARSme20191218.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.groups.ave-std.summary",
                   header=T, stringsAsFactors = F)

# can match == or doesn't match !=
# alph1 <- filter(alph, method=="ave")
alph <- filter(alph, method != "std")

expdata <- read.table(file="../may18ws.env.txt", 
                      header = T, stringsAsFactors = T)
samples <- read.table(file="../may18ws.sample.txt",
                      header=T, stringsAsFactors = F)
str(expdata)
expdata$Sample <- as.character(expdata$Sample)
str(expdata)

expdata <- left_join(expdata, samples, by="Sample")

alpha.expdata <- left_join(alph, expdata, by="group")

# read in beta diversity

jc <- parseDistanceDF("../MARSme20191218.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.jest.0.03.lt.ave.dist")
tyc <- parseDistanceDF("../MARSme20191218.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.braycurtis.0.03.lt.ave.dist")
tyc <- parseDistanceDF("../MARSme20191218.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.thetayc.0.03.lt.ave.dist")

# read in taxonomy

taxa <- read.table(textConnection(gsub("\\(.+?\\);", "\t", 
                readLines("../MARSme20191218.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.0.03.cons.taxonomy"))), 
                col.names=c("OTU","Size", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), skip=1)

# filtering taxa based on column headers of otu
taxa <- taxa[taxa$OTU %in% names(otu), ]

# get the subsampled abundance for each otu in the taxa table
sub.size <- data.frame(OTU = names(otu), size.sub = colSums(otu))

taxa <- left_join(taxa, sub.size, by = "OTU")

#### alpha diversity ####

ggplot(data=alpha.expdata, aes(x=sobs))+
    geom_histogram()

ggplot(data=alpha.expdata, aes(x=sobs, y=shannon, color=Site))+
    geom_point(size=3)

ggplot(data=alpha.expdata, aes(x=sobs, y=shannon, color=Type))+
    geom_point(size=3)

ggplot(data=alpha.expdata, aes(x=sobs, y=invsimpson, color=Site))+
    geom_point(size=3)

ggplot(data=alpha.expdata, aes(x=invsimpson, y=shannon, color=Site))+
    geom_point(size=3)

ggplot(data=alpha.expdata, aes(x=Type, y=sobs))+
    geom_boxplot()+
    geom_point()

ggplot(data=alpha.expdata, aes(x=Type, y=shannon))+
    geom_boxplot()+
    geom_point()


ggplot(data=alpha.expdata, aes(x=Type, y=invsimpson))+
    geom_boxplot()+
    geom_point()+
    labs(x="Sample Type", y="Inverse Simpson", 
         title="Alpha diversity Wier Pond")+
    theme_bw()
    


ggplot(data=alpha.expdata, aes(x=Type, y=sobs))+
    geom_violin()+
    geom_point()

site.aov <- aov(alpha.expdata$sobs~alpha.expdata$Site)

site.aov$qr$qr

summary(site.aov)

type.aov <- aov(alpha.expdata$invsimpson~alpha.expdata$Type)
summary(type.aov)

TukeyHSD(type.aov)
TukeyHSD(aov(alpha.expdata$shannon~alpha.expdata$Type))

### how to get help

?aov
?anova.lm

??nonparametric



#### Beta Diversity #####

jc.nms <- metaMDS(as.dist(jc), k=2, trymin=50, trymax=200, wascores=F)
jc.points <- data.frame(jc.nms$points)
jc.plot <- ggplot(data=jc.points, aes(x=MDS1, y= MDS2, 
                                      label=rownames(jc)))

x <- max(jc.points$MDS1)/1.5
y <- min(jc.points$MDS2)*2

jc.plot + geom_point(aes(color = alpha.expdata$Type, 
                         shape=alpha.expdata$Site), size =3)+
    # geom_text()+
    annotate("text", x, y, label = paste( "stress = ", 
                        round(jc.nms$stress, digits = 3)))+
    theme_bw()+
    stat_ellipse(aes(color=alpha.expdata$Type))+
    scale_color_manual(values = c("darkgrey", "brown", "blue"),
                     name = "Sample Type"  )+
    scale_shape_manual(values = c(16, 17), name = "Site" )+
    labs(title = "Jaccard NMS")+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
ggsave(file="jc.nms.jpg", height = 5, width= 9)   



bc.nms <- metaMDS(as.dist(bc), k=2, trymin=50, trymax=200, wascores=F)
bc.points <- data.frame(bc.nms$points)
bc.plot <- ggplot(data=bc.points, aes(x=MDS1, y= MDS2, 
                                      label=rownames(bc)))

x <- max(bc.points$MDS1)/1.5
y <- min(bc.points$MDS2)*2

bc.plot + geom_point(aes(color = alpha.expdata$Type, 
                         shape=alpha.expdata$Site), size =3)+
    # geom_text()+
    annotate("text", x, y, label = paste( "stress = ", 
                                          round(jc.nms$stress, digits = 3)))+
    theme_bw()+
    stat_ellipse(aes(color=alpha.expdata$Type))+
    scale_color_manual(values = c("darkgrey", "brown", "blue"),
                       name = "Sample Type"  )+
    scale_shape_manual(values = c(16, 17), name = "Site" )+
    labs(title = "Bray Curtis NMS")+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
ggsave(file="bc.nms.jpg", height = 5, width= 9)


tyc.nms <- metaMDS(as.dist(tyc), k=2, trymin=50, trymax=200, wascores=F)
tyc.points <- data.frame(tyc.nms$points)
tyc.plot <- ggplot(data=tyc.points, aes(x=MDS1, y= MDS2, 
                                      label=rownames(tyc)))

x <- max(tyc.points$MDS1)/1.5
y <- min(tyc.points$MDS2)*2

tyc.plot + geom_point(aes(color = alpha.expdata$Type, 
                          shape=alpha.expdata$Site), size =3)+
    # geom_text()+
    annotate("text", x, y, label = paste( "stress = ", 
                                          round(jc.nms$stress, digits = 3)))+
    theme_bw()+
    stat_ellipse(aes(color=alpha.expdata$Type))+
    scale_color_manual(values = c("darkgrey", "brown", "blue"),
                       name = "Sample Type"  )+
    scale_shape_manual(values = c(16, 17), name = "Site" )+
    labs(title = "Theta YC NMS")+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
ggsave(file="tyc.nms.jpg", height = 5, width= 9)    


#### hypothesis testing on beta diversity ####

permanova <- adonis(as.dist(jc) ~ alpha.expdata$Type, perm = 99,
                    rm.na=TRUE)
permanova

(permanova <- adonis(as.dist(bc) ~ alpha.expdata$Type, perm = 99,
                    rm.na=TRUE))

(permanova <- adonis(as.dist(tyc) ~ alpha.expdata$Type, perm = 99,
                    rm.na=TRUE))


#### indicator species ####

indic <- multipatt(otu[,-1], alpha.expdata$Type, 
                   control = how(nperm=99))
summary(indic)


write.csv(file="indicator.species.csv", 
          indic$sign %>% 
            rownames_to_column( var = "OTU") %>%
            mutate(p.fdr = round(p.adjust(p.value, "fdr"), 3)) %>%
            right_join(taxa, by= "OTU") %>%
            # filter(p.fdr < 0.6) %>%
            arrange(index))
          
              
              
              
              
            

