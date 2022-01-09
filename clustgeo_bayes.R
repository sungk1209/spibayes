# *------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: clustgeo_bayes.R
# | DATE: Sep.22.2020
# | CREATED BY:  Kay Sung     
# *----------------------------------------------------------------
# | PURPOSE: Cluster the drought trend using D.Ward clustering method
# |         Input is the Bayesian GAM modeling results
# |
# |
# *------------------------------------------------------------------

require(tidyverse)
require(devtools)
#devtools::install_github('kevinblighe/PCAtools')
require(ClustGeo)
require(factoextra) 
require(cluster)
require(spData)# contains datasets used in this example
require(sf)     # for handling spatial (vector) data
require(tmap)   # map making package: Used to print maps

#require(factoextra)

select <- dplyr::select

setwd("C:/Users/kyungmin/Documents/SPIcmdstan/Spatial_ondemand/final_analysis")
#data_path <- ("C:/Users/kyungmin/Documents/SPIcmdstan/output/betas")
#output_path <- ("C:/Users/kyungmin/Documents/SPIcmdstan/output")

station <- read.csv("GHCNDgauges.csv")
station_list <- station %>%
  filter(variable == "PRCP") %>% 
  filter(X < 1920 & X.1 > 2018) %>%
  mutate(id = as.character(ID)) %>% 
  select(id, Latitude,Longitude)

#################load file
#####I usually create the beta values in the server
load("station_clust.rda")

res.pca <- prcomp(temp[,2:ncol],scale = TRUE)
fviz_eig(res.pca, addlabels = TRUE)
fviz_pca_var(res.pca, col.var = "black")

# d0 <- dist(station_list[,4:159]) 
# tree <- hclustgeo(d0)
# station_list <- station_list %>%
#   mutate(clust = (cutree(tree,10)))
# 
# #Plot clustered trees
# plot(tree,hang=-1,label=FALSE, xlab="",sub="",
#      main="Ward dendrogram with D0 only",cex.main=0.8,cex=0.8,cex.axis=0.8,cex.lab=0.8)
# rect.hclust(tree,k=10,border=seq(1,10,by = 1))
# legend("topright", legend= paste("cluster",1:10), fill=1:10, cex=0.8,bty="n",border="white")
# 
# coordi <- as.data.frame(station_list[,2:3])
# city_label <- as.data.frame(station_list[,1])
# 
# names(P10) <- city_label
# sp::plot(coordi,col=P10,main="10 clusters partition obtained with D0 only",cex.main=0.8)
# legend("topright", legend=paste("cluster",1:10), fill=1:10, cex=0.8,bty="n",border="white")


###########################################################################
#Plotting with various tensor_ cyclic product
##################################################################################
#color brewers
colpal <- c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#dd1c77','#bc80bd')
#Generate us states map
us_states4326 = st_transform(us_states, 4326)

##################################
### s(tensor) +s(year) + s(tensor) with Ward.D clustering
############################################

#Normalize parameters
station_list <- station_clust %>%
  filter(Longitude > -130) %>%
  drop_na()

dim(station_list)

b_jdate <- station_list %>%
  select(contains("mean_jdate"))
b_year <- station_list %>%
  select(contains("mean_year"))
b_tensor <- station_list %>%
  select(contains("mean_tensor"))

b_disp_jdate <- station_list %>%
  select(contains("disp_jdate"))
b_disp_year <- station_list %>%
  select(contains("disp_year"))
b_disp_tensor <- station_list %>%
  select(contains("disp_tensor"))

b_jdate_norm <- apply(b_jdate, 1, scale)
b_year_norm <- apply(b_year, 1, scale)
b_tensor_norm <- apply(b_tensor, 1, scale)

b_djdate_norm <- apply(b_disp_jdate, 1, scale)
b_dyear_norm <- apply(b_disp_year, 1, scale)
b_dtensor_norm <- apply(b_disp_tensor, 1, scale)

### You need to transpose them

b_jdate_norm <- t(b_jdate_norm)
b_year_norm <- t(b_year_norm)
b_tensor_norm <- t(b_tensor_norm)

b_djdate_norm <- t(b_djdate_norm)
b_dyear_norm <- t(b_dyear_norm)
b_dtensor_norm <- t(b_dtensor_norm)

temp_norm <- cbind(b_year_norm)

#####clustering using Ward.D hierarchical method
nclust <- 8

d0 <- dist(temp_norm) 
tree <- hclustgeo(d0)

station_list <- station_list %>%
  mutate(cluster = (cutree(tree,nclust)))
##########################
k3 <- pam(temp_norm, k = 4, stand = TRUE)
k4 <- pam(temp_norm, k = 5, stand = TRUE)
k5 <- pam(temp_norm, k = 6, stand = TRUE)
k6 <- pam(temp_norm, k = 15, stand = TRUE)

# plots to compare
#p1 <- fviz_cluster(k2, geom = "point", data = df) + ggtitle("k = 2")
p3 <- fviz_cluster(k3, geom = "point",  data = temp_norm) + ggtitle("k = 4")
p4 <- fviz_cluster(k4, geom = "point",  data = temp_norm) + ggtitle("k = 5")
p5 <- fviz_cluster(k5, geom = "point",  data = temp_norm) + ggtitle("k = 6")
p6 <- fviz_cluster(k6, geom = "point",  data = temp_norm) + ggtitle("k = 7")

library(gridExtra)
grid.arrange( p3, p4,p5,p6, nrow = 2)
gap_stat <- clusGap(temp_norm, FUN = pam,K.max = 20, B = 20)
fviz_gap_stat(gap_stat)
#################################################################
#########Mapping
#############################################################

coordinates <- as_tibble(station_list)
coordi_sf <- st_as_sf(coordinates,coords = c("Longitude","Latitude"), crs= 4326)

us_states4326 = st_transform(us_states, 4326)

tmap_mode("plot")

plot20 <- tm_shape(us_states4326) + 
  tm_polygons() + 
  tm_layout(frame = FALSE) +
  tm_shape(coordi_sf) +
  tm_symbols(col = "cluster", scale = 0.5, n = 10,
             legend.col.is.portrait = FALSE,
             palette = colpal[1:nclust]) +
  tm_facets(by = "cluster", scale.factor = 4, ncol = 4) +
  tm_layout(main.title = "s(year) using Ward.D", 
            main.title.size = 1.5,
            legend.title.size = 1.0,
            legend.text.size = 1.5,
            frame = FALSE, bg.color = NA,
            legend.format = "g",
            inner.margins = 0.1) +
  tm_legend(outside = TRUE, outside.position = "bottom")

plot20  
tmap_save(plot20, "map_year8.png", height=8)

temp_norm <- as.data.frame(cbind(station_list$id, temp_norm, station_list$mean_no_b0))
njdate <- dim(b_jdate_norm)[2]
nyear <- dim(b_year_norm)[2]
ntensor <- dim(b_tensor_norm)[2]
cname <- c("id",
  paste0("mjdate_norm_", rep(1:njdate)),
  paste0("myear_norm_", rep(1:nyear)),
  paste0("mtensor_norm_", rep(1:ntensor)),
  "mean_b")
colnames(temp_norm) <-cname

# plot Betas used in the model
#mean_plot <- mean_yr1 %>% pivot_longer(names_to = "knots", values_to = "mean_b", cols = 3:dim(mean_yr1)[2]) %>%
#  mutate(yr_knots = as.factor(knots)) %>%
#  mutate(cluster  = as.factor(mean_b))

#p <- ggplot(mean_plot, aes(x = yr_knots)) + geom_line(aes(y = mean_b, group = id,color = cluster))+
#  scale_colour_manual(values = colpal) +ggtitle("s(year)+s(tensor) using PAM")
#scale_colour_brewer(palette = "Spectral")
#p
#ggsave(filename = "yr_tensors_pam.png", plot = p)

#clust1 <- station_list %>%
#  filter(mean_no_b0 == 1)
#apply(clust1[4:length(clust1)], 2, mean)

########################Get dendrogram########################
#https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/dendrogram
dend1 <- as.dendrogram(tree)
str(dend1, max = 2, last.str =  "'")
plot(dend1)
p <- plot(dend1, ylim = c(0,0.04))

nP <- list(col = 3:2, cex = c(2.0, 0.75), pch =  21:22,
           bg =  c("light blue", "pink"),
           lab.cex = 0.75, lab.col = "tomato")
p <- plot(dend1, nodePar= nP, edgePar = list(col = "gray", lwd = 2), horiz = TRUE)

#############################################################
# s(jdate) + s(tensor) for mean
######################################################
temp_norm <- b_jdate_norm

d0 <- dist(temp_norm) 
tree <- hclustgeo(d0)
sum(tree$height)

station_list <- station_list %>%
  mutate(mean_season = (cutree(tree,nclust)))
coordinates <- as_tibble(station_list)
coordi_sf <- st_as_sf(coordinates,coords = c("Longitude","Latitude"), crs= 4326)

tmap_mode("plot")

plot30 <- tm_shape(us_states4326) +
  tm_polygons() + 
  tm_layout(frame = FALSE) +
  tm_shape(coordi_sf) +
  tm_symbols(col = "mean_season", scale = .5, n = 10,
             legend.col.is.portrait = FALSE,
             palette = colpal[1:nclust]) +
  tm_facets(by = "mean_season", scale.factor = 4, ncol = 3) +
  tm_layout(main.title = "s(jdate)", 
            title.position = c("left", "TOP"),
            main.title.size = 1.0,
            legend.title.size = 2.0,
            frame = FALSE, bg.color = NA,) +
  tm_legend(outside = TRUE, outside.position = "bottom")

plot30  

tmap_save(plot30, "s_jdate.png", height=5)
##########################################################
##                  plot mean_year betas              ####
##########################################################
mean_yr1 <- station_list %>%
  select(c(id, cluster, contains("b_mean_year")))

pattern <- "b_mean_"
mean_plot <- mean_yr1 %>% pivot_longer(names_to = "knots", values_to = "mean_b", cols = 3:9) %>%
  mutate(yr_knots = as.factor(knots)) %>%
  mutate(cluster  = as.factor(cluster)) %>%
  mutate(yr_knots = str_remove(knots,fixed(pattern)))

p <- ggplot(mean_plot, aes(x = yr_knots)) + geom_line(aes(y = mean_b, group = id), color = 'blue') + 
  facet_wrap(.~cluster) +  scale_x_discrete(expand = c(0, 0)) + 
  ggtitle("b(year) of each clusters")

p
ggsave(filename = "byear_betas6.png", plot = p)

#Boxplot

p <- ggplot(mean_plot, aes(x=yr_knots, y=mean_b))  + geom_line(aes(colour = 'blue', group = id), color = "#5A809E")+
  geom_boxplot()+ 
  facet_wrap(.~cluster) + ggtitle("beta (year) of each clusters")

p
ggsave(filename = "box_byear8.png", plot = p)


#####################################################
#Inspect the outlier
#####################################################
row <- which(mean_yr1$`b_mean_year[4]`==min(mean_yr1$`b_mean_year[4]`))
idn <- mean_yr1$id[row]

a <- station_list %>%
  filter (id == idn)

coordinates <- as_tibble(a)
coordi_sf <- st_as_sf(coordinates,coords = c("Longitude","Latitude"), crs= 4326)

plot70 <- tm_shape(us_states4326) +
  tm_polygons() + 
  tm_layout(frame = FALSE) +
  tm_shape(coordi_sf) +
  tm_symbols(col = "mean_yr", scale = .5, n = 10,
             legend.col.is.portrait = FALSE,
             palette = colpal[1:nclust]) +
  tm_layout(main.title = "inspection", 
            title.position = c("left", "TOP"),
            main.title.size = 1.0,
            legend.title.size = 2.0,
            frame = FALSE, bg.color = NA,) +
  tm_legend(outside = TRUE, outside.position = "bottom")

plot70  

######################################################################
#####################
#######################################################################

plot_gauges <- tm_shape(us_states4326) + 
  tm_polygons() + 
  tm_layout(frame = FALSE) +
  tm_shape(coordi_sf) +
  tm_symbols(col = "blue", scale = 0.7, n = 10,
             legend.col.is.portrait = FALSE,
             palette = colpal[1:nclust]) +
  tm_layout(main.title = "gauge location", 
            main.title.size = 1.5,
            legend.title.size = 1.0,
            legend.text.size = 1.5,
            frame = FALSE, bg.color = NA,
            legend.format = "g",
            inner.margins = 0.1) +
  tm_legend(outside = TRUE, outside.position = "bottom")

plot_gauges 
tmap_save(plot_gauges, "map_gauges.png", height=8)

