# R script for a mesocosm experiment
# on phytoplankton colonisation and variability.
# Necessary R packages
library(ggcyto);library(flowViz);library(flowStats);library(flowDensity)
library(stringr);library(grid);library(gridExtra);library(MASS);library(vegan)
library(plyr);library(dunn.test)
library(FactoMineR);library(factoextra);library(scales)

# The mesocosm system consisted of 12 mesocosm tanks
# treatment indices
hightemp.select <- c(1,3,5,8,9,10) # index number of tanks kept at 3°C above ambient temperature
lowtemp <- 1:12
lowtemp <- lowtemp[-hightemp.select] # index number of tanks kept at ambient temperature

# Flow cytometric data analysis
# (1) Import raw data
# Files of different date are organized into separate folders
setwd("") # add file path
mesoflow19 <- vector(mode = "list", length = length(dir()))
for(i in 1:length(dir())) {
  setwd("") # add file path
  names(mesoflow19)[i] <- paste("2019-",dir()[i], sep = "")
  setwd(paste("/", dir()[i], sep="")) # add file path
  mesoflow19[[i]] <- sapply(dir(), read.FCS, USE.NAMES = TRUE)
  mesoflow19[[i]] <- as(mesoflow19[[i]], "flowSet")
  str_sub(sampleNames(mesoflow19[[i]]), start = -4, end = -1) <- ""
  mesoflow19[[i]]@phenoData@data$name <- sampleNames(mesoflow19[[i]])
}

dates <- c("2019-04-08","2019-04-11","2019-04-14","2019-04-17","2019-04-20",
           "2019-04-23","2019-04-26","2019-04-29","2019-05-02","2019-05-05",
           "2019-05-08","2019-05-11","2019-05-14","2019-05-17","2019-05-20")

# (2) Convert into GatingSet objects and log transform data
gs.meso <-  vector(mode = "list", length = length(mesoflow19))
for(i in 1:length(mesoflow19)) {
  gs.meso[[i]] <- GatingSet(mesoflow19[[i]])
}
names(gs.meso) <- dates

power10 <- function(x) {
  original <- 10^x
  return(original)
}

trans.func <- log10
inv.func <- power10
trans.obj <- trans_new("logtrans", trans.func, inv.func)
trans.list <- transformerList(c("FSC","SSC","FL2","FL3","FL4"), trans.obj)

for(i in 1:length(gs.meso)) {
  gs.meso[[i]] <- transform(gs.meso[[i]], trans.list)
}

# (3) Gating (noise filtering)
rg1 <- rectangleGate("FL3" = c(0.5,Inf), filterId = "chlorophyll")
lapply(gs.meso, FUN = gs_pop_add, gate = rg1)

lapply(gs.meso, recompute)

# (4) Binning functions
# (Step 1) Creating bin matrix
bins <- function(ll,ul,numbin) {
  boundaries <- seq(from = ll, to = ul, length.out = (numbin+1))
  bin.matrix <- matrix(data = c(boundaries[1:numbin],boundaries[-1]),
                       nrow = numbin, ncol = 2, byrow = F)
  return(bin.matrix)
}

# (Step 2) Extract GatingHierarchy data
prepare <- function(gh) {
  gh.matrix <- as.matrix(exprs(gh_pop_get_data(gh, y = "chlorophyll"))[,c(1,5)])
  gh.matrix <- as.matrix(data.frame(gh.matrix,
                                    rep(NA, times = dim(gh.matrix)[1]),
                                    rep(NA, times = dim(gh.matrix)[1])))
  colnames(gh.matrix)[3:4] <- c("FSC.bin","FL3.bin")
  return(gh.matrix)
}

# (Step 3) Assigning bin label to each cell
binning <- function(vec, bins) {
  tot.counter <- 1
  bincounter <- 1
  bin.index <- vector(length = length(vec))
  while(tot.counter <= length(vec)) {
    while(vec[tot.counter] < bins[bincounter,2]) {
      bin.index[tot.counter] <- bincounter
      tot.counter <- tot.counter +1
      if(tot.counter > length(vec)) {
        break
      }
    }
    bincounter <- bincounter +1
  }
  return(bin.index)
}

# (Step 4) Create matrix of bin counts
bin.matrix <- function(x, numbin1, numbin2) {
  temp.matrix <- matrix(data = 0, nrow = numbin1, ncol = numbin2)
  for(i in 1:numbin1) {
    temp.vector <- x[x[,1] == i,2]
    bin.counts <- tapply(temp.vector, as.factor(temp.vector), FUN = length)
    temp.matrix[i,as.numeric(names(bin.counts))] <- bin.counts
  }
  return(temp.matrix)
}

# (Step 5) delete empty bins
del.empty <- function(x) {
  col.total <- colSums(x)
  non.zeros <- which(col.total != 0)
  return(x[,non.zeros])
}

# (5) Binning the data
# create a complex list of bin indices
# within each gating hierarchy (within each sample)
bin.indices <- vector(mode = "list", length = length(gs.meso))
names(bin.indices) <- names(gs.meso)

# create fingerprint matrices and vectors
fp.matrices <- vector(mode = "list", length = length(gs.meso))
names(fp.matrices) <- names(gs.meso)

for(i in 1:length(gs.meso)) {
  bin.indices[[i]] <- vector(mode = "list", length = 12)
  fp.matrices[[i]] <- vector(mode = "list", length = 12)
}

# Setting bin boundaries
fsc.bins <- bins(0,3,64)
fl3.bins <- bins(0.5,3,64)
# assigning bin number to cell-level data
# using the newly defined functions
for(i in 1:length(gs.meso)) {
  for(j in 1:length(bin.indices[[i]])) {
    gh.temp <- prepare(gs.meso[[i]][[j]])
    gh.temp <- gh.temp[order(gh.temp[,1]),]
    gh.temp[,3] <- binning(gh.temp[,1], fsc.bins)
    gh.temp <- gh.temp[order(gh.temp[,2]),]
    gh.temp[,4] <- binning(gh.temp[,2], fl3.bins)
    gh.temp <- gh.temp[order(gh.temp[,3],gh.temp[,4]),]
    bin.indices[[i]][[j]] <- gh.temp[,3:4]
  }
}
rm(gh.temp)

# filling up bin matrices in the FSC-FL3 space
for(i in 1:length(bin.indices)) {
  for(j in 1:length(fp.matrices[[i]])) {
    fp.matrices[[i]][[j]] <- bin.matrix(bin.indices[[i]][[j]], dim(fsc.bins)[1], dim(fl3.bins)[1])
  }
}

# convert bin matrices to bin vectors
fp.vectors <- vector(mode = "list", length(fp.matrices))
names(fp.vectors) <- names(fp.matrices)
for(i in 1:length(fp.matrices)) {
  fp.vectors[[i]] <- matrix(NA, nrow = 12, ncol = dim(fsc.bins)[1]*dim(fl3.bins)[1])
  rownames(fp.vectors[[i]]) <- 1:12
  for(j in 1:length(fp.matrices[[i]])) {
    temp <- NULL
    for(k in 1:dim(fp.matrices[[i]][[j]])[1]) {
      temp <- c(temp, fp.matrices[[i]][[j]][k,])
    }
    fp.vectors[[i]][j,] <- temp
  }
}

# delete empty bins from each date
for(i in 1:length(fp.vectors)) {
  fp.vectors[[i]] <- del.empty(fp.vectors[[i]])
}

# log transform bin counts
fp.log.vectors <- vector(mode = "list", length(fp.matrices))
for(i in 1:length(fp.vectors)) {
  fp.log.vectors[[i]] <- log10(fp.vectors[[i]] + 1)
}
names(fp.log.vectors) <- names(fp.vectors)

# determine distances (dissimilarity)
# between the samples for each date
# using the Bray-Curtis index
# and determine compositional variation
# from the distances using betadisper()
fp.log.dist <- vector(mode = "list", length = length(dates))
fp.log.disp <- vector(mode = "list", length = length(dates))
names(fp.log.dist) <- names(fp.log.vectors)
names(fp.log.disp) <- names(fp.log.vectors)
treatment.vector <- rep("Ambient.Temp", times = 12)
treatment.vector[hightemp.select] <- "Higher.Temp"
for(i in 1:length(fp.log.vectors)) {
  fp.log.dist[[i]] <- vegdist(fp.log.vectors[[i]], method = "bray", binary = FALSE)
  fp.log.disp[[i]] <- betadisper(fp.log.dist[[i]], group = treatment.vector,
                                 type = "centroid", bias.adjust = TRUE)
}

distances.log.df <- as.data.frame(matrix(data = NA, nrow = 12 * length(dates), ncol = 5))
colnames(distances.log.df) <- c("Day","Date","Mesocosm","Treatment","Distance")
distances.log.df$Day <- rep(seq(from = 0, to = 42, by = 3), each = 12)
distances.log.df$Date <- rep(dates, each = 12)
distances.log.df$Mesocosm <- rep(sampleNames(gs.meso[[1]]), times = length(dates))
distances.log.df$Treatment <- rep(treatment.vector, times = length(dates))
distances.log.df$Treatment <- as.factor(distances.log.df$Treatment)

for(i in 1:length(dates)) {
  distances.log.df$Distance[distances.log.df$Date == dates[i]] <- fp.log.disp[[i]]$distances
}

# Fig 1
# Cytometric scatterplots of the communities
# for a specific date
plotlabels <- c(expression(paste("10"^"1")),
                expression(paste("10"^"2")),
                expression(paste("10"^"3")))
flowplot <- vector(mode = "list", length = 12)
for(k in 1:12) {
  if(k %in% hightemp.select == TRUE) {
    flowplot[[k]] <- as.ggplot(ggcyto(gs.meso[[which(names(gs.meso) == "2019-05-05")]][[k]],
                                      aes(x = "FSC", y = "FL3"),
                                      subset = "chlorophyll") +
                                 geom_hex(bins = 200) +
                                 theme_bw() +
                                 theme(panel.grid = element_blank(), aspect.ratio = 1,
                                       panel.background = element_rect(fill = "lightblue2"),
                                       axis.title = element_blank(), plot.title = element_blank(),
                                       plot.margin = unit(c(0.5,0,0.5,0),"mm"),
                                       legend.position = "none") +
                                 ggcyto_par_set(limits = list(x= c(0,3), y = c(0.5,3))) +
                                 scale_x_continuous(labels = plotlabels, breaks = c(1,2,3)) +
                                 scale_y_continuous(labels = plotlabels, breaks = c(1,2,3)) +
                                 annotate(geom = "text", x=0.3, y=2.9,
                                          label = "+3°C", size = 4, fontface = "bold" ))
    } else {
      flowplot[[k]] <- as.ggplot(ggcyto(gs.meso[[which(names(gs.meso) == "2019-05-05")]][[k]],
                                        aes(x = "FSC", y = "FL3"),
                                        subset = "chlorophyll") +
                                   geom_hex(bins = 200) +
                                   theme_bw() +
                                   theme(panel.grid = element_blank(), aspect.ratio = 1,
                                         panel.background = element_rect(fill = "white"),
                                         axis.title = element_blank(), plot.title = element_blank(),
                                         plot.margin = unit(c(0.5,0,0.5,0),"mm"),
                                         legend.position = "none") +
                                   ggcyto_par_set(limits = list(x= c(0,3), y = c(0.5,3))) +
                                   scale_x_continuous(labels = plotlabels, breaks = c(1,2,3)) +
                                   scale_y_continuous(labels = plotlabels, breaks = c(1,2,3)))
      }
}
  
ordered.plots <- arrangeGrob(grobs = flowplot, nrow = 4,
                                     bottom = textGrob("FSC", gp = gpar(fontface = "bold", fontsize = 15)),
                                     left = textGrob("FL3", gp = gpar(fontface = "bold", fontsize = 15), rot = 90))

setwd("") # add file path
tiff(filename = "Fig_1.tiff",
     width = 2100, height = 2500, res = 300)
grid.draw(ordered.plots)
dev.off()

# Fig 2
# Representation of the binning process
plot1 <- as.ggplot(ggcyto(gs.meso[[which(names(gs.meso) == "2019-05-05")]][[5]],
                 aes(x = "FSC", y = "FL3"),
                 subset = "chlorophyll") +
            geom_hex(bins = 200) +
            theme_bw() +
            theme(aspect.ratio = 1, panel.ontop = TRUE,
                  panel.grid.major = element_line(colour = "black", size = 0.1),
                  panel.grid.minor = element_line(colour = "black", size = 0.1),
                  panel.background = element_blank(),
                  axis.title = element_text(size = 10),
                  axis.text = element_text(size = 6),
                  plot.title = element_blank(),
                  plot.margin = unit(c(0.5,0,0.5,0),"mm"),
                  legend.position = "none") +
            xlab("FSC") +
            ylab("FL3") +
            ggcyto_par_set(limits = list(x= c(0,3), y = c(0.5,3))) +
            scale_x_continuous(labels = plotlabels, breaks = c(1,2,3), minor_breaks = seq(0,3,0.2)) +
            scale_y_continuous(labels = plotlabels, breaks = c(1,2,3), minor_breaks = seq(0.5,3,0.125)))

sample.vector <- fp.vectors[[which(names(gs.meso) == "2019-05-05")]][5,]
bin.index <- seq(1,dim(fp.vectors[[which(names(gs.meso) == "2019-05-05")]])[2])
sample.vector <- data.frame(bin.index,sample.vector)

plot2 <- ggplot(data = sample.vector, mapping = aes(x = bin.index, y = sample.vector)) +
  theme_bw() +
  theme(panel.grid = element_blank(), aspect.ratio = 0.6,
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 6),
        plot.margin = unit(c(0.5,0,0.5,20),"mm")) +
  xlab("Bin index") +
  ylab("Cell count") +
  geom_line(colour = "black", size = 0.25)

figure2 <- arrangeGrob(grobs = list(plot1,plot2), nrow = 1,
                       heights = unit(11, units = "cm"), widths = unit(c(5,11), units = "cm"))
setwd("") # add file path
tiff(filename = "Fig_2.tiff",
     width = 1930, height = 720, res = 300)
grid.draw(figure2)
dev.off()

# Mean and standard deviation of centroid distances
dist.mean <- as.data.frame(tapply(distances.log.df$Distance,
                                  INDEX = list(distances.log.df$Date,distances.log.df$Treatment),
                                  FUN = mean))
dist.sd <- as.data.frame(tapply(distances.log.df$Distance,
                                INDEX = list(distances.log.df$Date,distances.log.df$Treatment),
                                FUN = sd))
beta.div <- as.data.frame(matrix(data = NA, nrow = 2 * length(dates), ncol = 5))
colnames(beta.div) <- c("Date","Day","Treatment","Mean","sd")
beta.div$Date <- rep(dates, times = 2)
beta.div$Day <- rep(unique(distances.log.df$Day), times = 2)
beta.div$Treatment <- rep(colnames(dist.mean), each = length(dates))
beta.div$Mean <- c(dist.mean$Ambient.Temp,dist.mean$Higher.Temp)
beta.div$sd <- c(dist.sd$Ambient.Temp,dist.sd$Higher.Temp)
beta.div <- beta.div[beta.div$Day >= 12,]

# Descriptive statistics on FSC and FL3 data
cyto.stat <- as.data.frame(matrix(data = NA, nrow = 12 * length(dates), ncol = 12))
colnames(cyto.stat) <- c("Day","Date","Mesocosm","Treatment",
                         "FSC.mean","FSC.mode","FL3.mean","FL3.mode",
                         "Cyto.Diversity","Cyto.Evenness","Tax.Diversity","Tax.Evenness")
cyto.stat$Day <- rep(seq(from = 0, to = 42, by = 3), each = 12)
cyto.stat$Date <- distances.log.df$Date
cyto.stat$Mesocosm <- distances.log.df$Mesocosm
cyto.stat$Treatment <- as.factor(distances.log.df$Treatment)

extract <- function(gh) {
  gh.matrix <- as.data.frame(exprs(gh_pop_get_data(gh, y = "chlorophyll"))[,c(1,5)])
  colnames(gh.matrix) <- c("FSC","FL3")
  return(gh.matrix)
}

for(i in 1:length(gs.meso)) {
  for(j in 1:length(gs.meso[[i]])) {
    gh.temp <- extract(gs.meso[[i]][[j]])
    fsc.heights <- getPeaks(obj = gh.temp$FSC)$P.h
    fl3.heights <- getPeaks(obj = gh.temp$FL3)$P.h
    cyto.stat$FSC.mean[cyto.stat$Date == dates[i] & cyto.stat$Mesocosm == sampleNames(gs.meso[[i]])[j]] <-
      mean(gh.temp$FSC)
    cyto.stat$FSC.mode[cyto.stat$Date == dates[i] & cyto.stat$Mesocosm == sampleNames(gs.meso[[i]])[j]] <-
      getPeaks(obj = gh.temp$FSC)$Peaks[fsc.heights == max(fsc.heights)]
    cyto.stat$FL3.mean[cyto.stat$Date == dates[i] & cyto.stat$Mesocosm == sampleNames(gs.meso[[i]])[j]] <-
      mean(gh.temp$FL3)
    cyto.stat$FL3.mode[cyto.stat$Date == dates[i] & cyto.stat$Mesocosm == sampleNames(gs.meso[[i]])[j]] <-
      getPeaks(obj = gh.temp$FL3)$Peaks[fl3.heights == max(fl3.heights)]
    rm(list = c("gh.temp","fsc.heights","fl3.heights"))
  }
}

# Diversity and evenness
even <- function(x) {
  div <- diversity(x, index = "shannon")
  species <- specnumber(x)
  evenness <- div/log(species)
  return(evenness)
}

for(i in 1:length(fp.vectors)) {
  cyto.stat$Cyto.Diversity[cyto.stat$Date == dates[i]] <-
    diversity(fp.vectors[[i]], index = "shannon")
  cyto.stat$Cyto.Evenness[cyto.stat$Date == dates[i]] <-
    even(fp.vectors[[i]])
}

# Analyzing taxonomic diversity and variability
# Taxonomic diversity and evenness
dates.taxon <- c("2019-04-20","2019-04-26","2019-05-02","2019-05-08","2019-05-14","2019-05-20")
community <- list()
setwd("") # add file path
for(i in 1:6) {
  community[[i]] <- read.csv(file = paste("Phytoplankton_",dates.taxon[i],".csv", sep = ""),
                             header=TRUE, row.names=1, sep=",")
  community[[i]] <- t(community[[i]])
  rownames(community[[i]]) <- 1:12
}
names(community) <- dates.taxon

for(i in 1:length(community)) {
  cyto.stat$Tax.Diversity[cyto.stat$Date == dates.taxon[i]] <-
    diversity(community[[i]], index = "shannon")
  cyto.stat$Tax.Evenness[cyto.stat$Date == dates.taxon[i]] <-
    even(community[[i]])
}


# Taxonomic beta diversity
tax.dist <- vector(mode = "list", length = length(dates.taxon))
tax.disp <- vector(mode = "list", length = length(dates.taxon))
names(tax.dist) <- dates.taxon
names(tax.disp) <- dates.taxon
for(i in 1:length(dates.taxon)) {
  tax.dist[[i]] <- vegdist(community[[i]], method = "bray", binary = FALSE)
  tax.disp[[i]] <- betadisper(tax.dist[[i]], group = treatment.vector,
                              type = "centroid", bias.adjust = TRUE)
}

distances.tax <- as.data.frame(matrix(data = NA, nrow = 12 * length(dates.taxon), ncol = 5))
colnames(distances.tax) <- c("Day","Date","Mesocosm","Treatment","Distance")
distances.tax$Day <- rep(seq(from = 12, to = 42, by = 6), each = 12)
distances.tax$Date <- rep(dates.taxon, each = 12)
distances.tax$Mesocosm <- rep(sampleNames(gs.meso[[1]]), times = length(dates.taxon))
distances.tax$Treatment <- rep(treatment.vector, times = length(dates.taxon))
distances.tax$Treatment <- as.factor(distances.tax$Treatment)

for(i in 1:length(dates.taxon)) {
  distances.tax$Distance[distances.tax$Date == dates.taxon[i]] <- tax.disp[[i]]$distances
}

cyto.stat$Cyto.Distances <- distances.log.df$Distance
for(i in 1:length(dates.taxon)) {
  cyto.stat$Tax.Distances[cyto.stat$Date == dates.taxon[i]] <- distances.tax$Distance[distances.tax$Date == dates.taxon[i]]
}


# Statistical tests
dunn.p <- as.data.frame(matrix(data = NA, nrow = length(dates), ncol = 12))
colnames(dunn.p) <- colnames(cyto.stat)[c(1,2,5:14)]
dunn.p$Day <- unique(cyto.stat$Day)
dunn.p$Date <- dates

for(i in 3:dim(dunn.p)[2]) {
  for(j in 1:length(dates)) {
    temp <- cyto.stat[cyto.stat$Date == dates[j],c(4,i+2)]
    if(is.na(sum(temp[,2])) == FALSE) {
      dunn.p[j,i] <- dunn.test(temp[,2], g = temp[,1])$P
    } else {
      dunn.p[j,i] <- NA
    }
    rm("temp")
  }
}

dunn.p <- dunn.p[dunn.p$Day >= 12,]
dunn.p

# Plots
significant.FL3mean <- dunn.p$Day[dunn.p$FL3.mean < 0.05]
significant.FL3mode05 <- dunn.p$Day[dunn.p$FL3.mode < 0.05 & dunn.p$FL3.mode >= 0.01]
significant.FL3mode01 <- dunn.p$Day[dunn.p$FL3.mode < 0.01]
significant.cyto.var05 <- dunn.p$Day[dunn.p$Cyto.Distances < 0.05 & dunn.p$Cyto.Distances >= 0.01]
significant.cyto.var01 <- dunn.p$Day[dunn.p$Cyto.Distances < 0.01]
significant.tax.var05 <- dunn.p$Day[dunn.p$Tax.Distances < 0.05 & dunn.p$Tax.Distances >= 0.01]
significant.tax.var01 <- dunn.p$Day[dunn.p$Tax.Distances < 0.01]
significant.tax.var05 <- significant.tax.var05[-which(is.na(significant.tax.var05))]
significant.tax.var01 <- significant.tax.var01[-which(is.na(significant.tax.var01))]

# Means and SD of all variables
# one list / index
cyto.stat <- cyto.stat[cyto.stat$Day >= 12,]
mean.sd <- vector(mode = "list", length = 10)
names(mean.sd) <- colnames(cyto.stat)[5:14]
for(i in 1:length(mean.sd)) {
  mean.sd[[i]] <- as.data.frame(matrix(NA, nrow = length(unique(cyto.stat$Day))*2, ncol = 4))
  colnames(mean.sd[[i]]) <- c("Day","Treatment","Mean","sd")
  mean.sd[[i]]$Day <- rep(unique(cyto.stat$Day), times = 2)
  mean.sd[[i]]$Treatment <- rep(c("Ambient temperature","Elevated temperature"),
                                each = length(unique(cyto.stat$Day)))
  mean.sd[[i]]$Treatment <- as.factor(mean.sd[[i]]$Treatment)
  temp.mean <- tapply(cyto.stat[,i+4], INDEX = list(cyto.stat$Day,cyto.stat$Treatment), FUN = mean)
  temp.sd <- tapply(cyto.stat[,i+4], INDEX = list(cyto.stat$Day,cyto.stat$Treatment), FUN = sd)
  mean.sd[[i]]$Mean <- c(temp.mean[,1],temp.mean[,2])
  mean.sd[[i]]$sd <- c(temp.sd[,1],temp.sd[,2])
  rm(list = c("temp.mean","temp.sd"))
}


# Figure on FSC and FL3 means and modes
treatment.color3 <- c("Ambient temperature" = "white", "Elevated temperature" = "black")
grob.FSC.mean <- ggplot(data = mean.sd[[1]], mapping = aes(x = Day, y = Mean,
                                                           group = Treatment)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "white"),
        axis.title.x = element_text(margin = margin(10,00,0,0), size = 20),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.position = c(0.33,0.88),
        legend.spacing.y = unit(0.1,"mm"),
        legend.background = element_rect(linetype = 1, size = 0.5, color = "black"),
        axis.text.x = element_text(size = 18, hjust = 0.4, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        aspect.ratio = 0.75) +
  geom_errorbar(aes(x = Day,
                    ymin = Mean - sd,
                    ymax = Mean + sd,
                    group = Treatment),
                width = 1, position = position_dodge(width = 1.75)) +
  geom_point(aes(fill = Treatment),
             size = 5, shape = 21, position = position_dodge(width = 1.75)) +
  labs(x = "", y = "log(FSC) mean") +
  scale_fill_manual(name = NULL, values = treatment.color3)

grob.FSC.mode <- ggplot(data = mean.sd[[2]], mapping = aes(x = Day, y = Mean,
                                                           group = Treatment)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "white"),
        axis.title.x = element_text(margin = margin(10,00,0,0), size = 20),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 20),
        axis.text.x = element_text(size = 18, hjust = 0.4, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        aspect.ratio = 0.75) +
  geom_errorbar(aes(x = Day,
                    ymin = Mean - sd,
                    ymax = Mean + sd,
                    group = Treatment),
                width = 1, position = position_dodge(width = 1.75)) +
  geom_point(aes(fill = Treatment),
             size = 5, shape = 21, position = position_dodge(width = 1.75),
             show.legend = FALSE) +
  labs(y = "log(FSC) mode") +
  scale_fill_manual(name = NULL, values = treatment.color3)

grob.FL3.mean <- ggplot(data = mean.sd[[3]], mapping = aes(x = Day, y = Mean,
                                                           group = Treatment)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "white"),
        axis.title.x = element_text(margin = margin(10,00,0,0), size = 20),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 20),
        axis.text.x = element_text(size = 18, hjust = 0.4, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        aspect.ratio = 0.75) +
  geom_errorbar(aes(x = Day,
                    ymin = Mean - sd,
                    ymax = Mean + sd,
                    group = Treatment),
                width = 1, position = position_dodge(width = 1.75)) +
  geom_point(aes(fill = Treatment),
             size = 5, shape = 21, position = position_dodge(width = 1.75),
             show.legend = FALSE) +
  labs(x = "", y = "log(FL3) mean") +
  scale_fill_manual(name = NULL, values = treatment.color3) +
  annotate("text", x = significant.FL3mean, y = 0.95, size = 10, label = "*")

grob.FL3.mode <- ggplot(data = mean.sd[[4]], mapping = aes(x = Day, y = Mean,
                                                           group = Treatment)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "white"),
        axis.title.x = element_text(margin = margin(10,00,0,0), size = 20),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 20),
        axis.text.x = element_text(size = 18, hjust = 0.4, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        aspect.ratio = 0.75) +
  geom_errorbar(aes(x = Day,
                    ymin = Mean - sd,
                    ymax = Mean + sd,
                    group = Treatment),
                width = 1, position = position_dodge(width = 1.75)) +
  geom_point(aes(fill = Treatment),
             size = 5, shape = 21, position = position_dodge(width = 1.75),
             show.legend = FALSE) +
  labs(y = "log(FL3) mode") +
  scale_fill_manual(name = NULL, values = treatment.color3) +
  annotate("text", x = significant.FL3mode05, y = 0.5, size = 10, label = "*") +
  annotate("text", x = significant.FL3mode01, y = 0.5, size = 10, label = "**")

setwd("") # add file path
tiff(filename = "Fig_4.tiff", width = 3500, height = 2800, res = 300)
grid.arrange(grobs = list(grob.FSC.mean,grob.FL3.mean,
                          grob.FSC.mode,grob.FL3.mode), nrow = 2)
dev.off()


# Figure on cytometric and taxonomic diversity
grob.cyto.div<- ggplot(data = mean.sd$Cyto.Diversity, mapping = aes(x = Day, y = Mean,
                                                                    group = Treatment)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "white"),
        axis.title.x = element_text(margin = margin(10,00,0,0), size = 20),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 20),
        axis.text.x = element_text(size = 18, hjust = 0.4, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        aspect.ratio = 0.75) +
  geom_errorbar(aes(x = Day,
                    ymin = Mean - sd,
                    ymax = Mean + sd,
                    group = Treatment),
                width = 1, position = position_dodge(width = 1.75)) +
  geom_point(aes(fill = Treatment),
             size = 5, shape = 21, position = position_dodge(width = 1.75),
             show.legend = FALSE) +
  labs(x = "", y = expression(paste(italic("H"["c"])))) +
  scale_fill_manual(name = NULL, values = treatment.color3)

grob.cyto.even <- ggplot(data = mean.sd$Cyto.Evenness, mapping = aes(x = Day, y = Mean,
                                                                     group = Treatment)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "white"),
        axis.title.x = element_text(margin = margin(10,00,0,0), size = 20),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.position = c(0.67,0.88),
        legend.spacing.y = unit(0.1,"mm"),
        legend.background = element_rect(linetype = 1, size = 0.5, color = "black"),
        axis.text.x = element_text(size = 18, hjust = 0.4, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        aspect.ratio = 0.75) +
  geom_errorbar(aes(x = Day,
                    ymin = Mean - sd,
                    ymax = Mean + sd,
                    group = Treatment),
                width = 1, position = position_dodge(width = 1.75)) +
  geom_point(aes(fill = Treatment),
             size = 5, shape = 21, position = position_dodge(width = 1.75)) +
  labs(y = expression(paste(italic("J"["c"])))) +
  scale_fill_manual(name = NULL, values = treatment.color3)

grob.tax.div <- ggplot(data = mean.sd$Tax.Diversity, mapping = aes(x = Day, y = Mean,
                                                                   group = Treatment)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "white"),
        axis.title.x = element_text(margin = margin(10,00,0,0), size = 20),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 20),
        axis.text.x = element_text(size = 18, hjust = 0.4, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        aspect.ratio = 0.75) +
  geom_errorbar(aes(x = Day,
                    ymin = Mean - sd,
                    ymax = Mean + sd,
                    group = Treatment),
                width = 1, position = position_dodge(width = 1.75)) +
  geom_point(aes(fill = Treatment),
             size = 5, shape = 21, position = position_dodge(width = 1.75),
             show.legend = FALSE) +
  labs(x = "", y = expression(paste(italic("H"["t"])))) +
  scale_fill_manual(name = NULL, values = treatment.color3)

grob.tax.even <- ggplot(data = mean.sd$Tax.Evenness, mapping = aes(x = Day, y = Mean,
                                                                   group = Treatment)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "white"),
        axis.title.x = element_text(margin = margin(10,00,0,0), size = 20),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 20),
        axis.text.x = element_text(size = 18, hjust = 0.4, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        aspect.ratio = 0.75) +
  geom_errorbar(aes(x = Day,
                    ymin = Mean - sd,
                    ymax = Mean + sd,
                    group = Treatment),
                width = 1, position = position_dodge(width = 1.75)) +
  geom_point(aes(fill = Treatment),
             size = 5, shape = 21, position = position_dodge(width = 1.75),
             show.legend = FALSE) +
  labs(y = expression(paste(italic("J"["t"])))) +
  scale_fill_manual(name = NULL, values = treatment.color3)

setwd("") # add file path
tiff(filename = "Fig_5.tiff", width = 3500, height = 2800, res = 300)
grid.arrange(grobs = list(grob.cyto.div,grob.tax.div,
                          grob.cyto.even,grob.tax.even), nrow = 2)
dev.off()

# Figure on compositional variability
grob.cyto.var <- ggplot(data = mean.sd$Cyto.Distances, mapping = aes(x = Day, y = Mean,
                                                                   group = Treatment)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "white"),
        axis.title.x = element_text(margin = margin(10,00,0,0), size = 20),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 20),
        axis.text.x = element_text(size = 18, hjust = 0.4, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        aspect.ratio = 0.75) +
  geom_errorbar(aes(x = Day,
                    ymin = Mean - sd,
                    ymax = Mean + sd,
                    group = Treatment),
                width = 1, position = position_dodge(width = 1.75)) +
  geom_point(aes(fill = Treatment),
             size = 5, shape = 21, position = position_dodge(width = 1.75),
             show.legend = FALSE) +
  labs(y = "Cytometric variation") +
  scale_fill_manual(name = NULL, values = treatment.color3) +
  annotate("text", x = significant.cyto.var05, y = 0.1, size = 10, label = "*") +
  annotate("text", x = significant.cyto.var01, y = 0.1, size = 10, label = "**")

grob.tax.var <- ggplot(data = mean.sd$Tax.Distances, mapping = aes(x = Day, y = Mean,
                                                                   group = Treatment)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "white"),
        axis.title.x = element_text(margin = margin(10,00,0,0), size = 20),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 20),
        axis.text.x = element_text(size = 18, hjust = 0.4, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.position = c(0.69,0.9),
        legend.spacing.y = unit(0.1,"mm"),
        legend.background = element_rect(linetype = 1, size = 0.5, color = "black"),
        aspect.ratio = 0.75) +
  geom_errorbar(aes(x = Day,
                    ymin = Mean - sd,
                    ymax = Mean + sd,
                    group = Treatment),
                width = 1, position = position_dodge(width = 1.75)) +
  geom_point(aes(fill = Treatment),
             size = 5, shape = 21, position = position_dodge(width = 1.75)) +
  labs(y = "Taxonomic variation") +
  lims(y = c(0.2,0.8)) +
  scale_fill_manual(name = NULL, values = treatment.color3) +
  annotate("text", x = significant.tax.var05, y = 0.2, size = 10, label = "*") +
  annotate("text", x = significant.tax.var01, y = 0.2, size = 10, label = "**")

setwd("") # add file path
tiff(filename = "Fig_6b.tiff", width = 3600, height = 1400, res = 300)
grid.arrange(grobs = list(grob.cyto.var,grob.tax.var), nrow = 1)
dev.off()

# Correspondence analysis
# Community data frame
community.df <- as.data.frame(community[[1]])
for(i in 2:length(community)) {
  community.df <- rbind(community.df,community[[i]])
}
treatment <- rep("ambient", times = 12)
treatment[hightemp.select] <- "elevated"
community.df$Treatment <- rep(treatment, times = length(community))
community.df$Date <- rep(names(community), each = length(community))
end <- dim(community.df)[2]
community.df <- community.df[,c(end,end-1,1:(end-2))]

# Omit rare taxa from the analysis
remove <- which(colnames(community.df) == "Chroococcus sp." |
                  colnames(community.df) == "Closterium sp." |
                  colnames(community.df) == "Cosmarium bioculatum" |
                  colnames(community.df) == "Cryptomonas erosa" |
                  colnames(community.df) == "Didymocystis inermis" |
                  colnames(community.df) == "Golenkiniopsis sp." |
                  colnames(community.df) == "Monoraphidium sp." |
                  colnames(community.df) == "Oocystis sp." |
                  colnames(community.df) == "Rhodomonas minuta" |
                  colnames(community.df) == "Schroederia sp." |
                  colnames(community.df) == "Staurastrum sp." |
                  colnames(community.df) == "Trachelomonas volvocina" |
                  colnames(community.df) == "Treubaria schmidlei")

community.df <- community.df[,-remove]

# Simplify names and merge closely related taxa
community.df$"centric diatoms" <- community.df$"Aulacoseira sp." +
  community.df$"Pantocsekiella ocellata"
community.df$"chrysophytes" <- community.df$"Chrysamoeba sp." +
  community.df$"Chrysoccoccus punctiformis"
community.df$"dinophytes" <- community.df$"Gymnodinium sp." +
  community.df$"Peridinium sp."
remove <- which(colnames(community.df) == "Aulacoseira sp." |
                  colnames(community.df) == "Pantocsekiella ocellata" |
                  colnames(community.df) == "Chrysamoeba sp." |
                  colnames(community.df) == "Chrysoccoccus punctiformis" |
                  colnames(community.df) == "Gymnodinium sp." |
                  colnames(community.df) == "Peridinium sp.")
community.df <- community.df[,-remove]
colnames(community.df)[colnames(community.df) == "Chlamydomonas sp."] <- "Chlamydomonas"
colnames(community.df)[colnames(community.df) == "Chrysochromulina parva"] <- "Chrysochromulina"
colnames(community.df)[colnames(community.df) == "Haematococcus pluvialis"] <- "Haematococcus"
colnames(community.df)[colnames(community.df) == "Kirchneriella obesa"] <- "Kirchneriella"
colnames(community.df)[colnames(community.df) == "Koliella longiseta"] <- "Koliella"
colnames(community.df)[colnames(community.df) == "Lagerheimia genevensis"] <- "Lagerheimia"
colnames(community.df)[colnames(community.df) == "Mallomonas sp."] <- "Mallomonas"
colnames(community.df)[colnames(community.df) == "Scenedesmus sp."] <- "Scenedesmus"

ca.community <- CA(community.df[,3:dim(community.df)[2]], graph = F)
summary(ca.community)

col.contrib <- facto_summarize(ca.community, element = "col", result = "contrib", axes = 1:2)
treatment.color <- c("ambient" = "grey75", "elevated" = "grey35")

tiff(filename = "Fig_7.tiff", width = 2500, height = 2500, res = 300)
grid.draw(fviz_ca_biplot(ca.community, axes = c(1,2), label = "none", arrows = c(F,F),
               map = "symmetric", col.col = "black", invisible = c("row","col"),
               labelsize = 5.5, pointsize = 4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 15),
        plot.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.position = c(0.85,0.90),
        legend.key = element_blank(),
        legend.background = element_rect(fill = NULL, colour = "black")) +
  xlim(-0.75,2.4) +
  geom_point(aes(x = ca.community$row$coord[,1],
                 y = ca.community$row$coord[,2],
                 colour = rep(treatment, times = length(unique(community.df$Date)))),
             size = 5, shape = 16) +
  geom_text(aes(label = rownames(ca.community$col$coord),
                x = ca.community$col$coord[,1],
                y = ca.community$col$coord[,2]),
            colour = "black", size = log10(col.contrib$contrib+5)*5,
            fontface = "bold") +
  scale_color_manual(name = "Temperature", values = treatment.color))
dev.off()
