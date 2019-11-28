library(devtools)
library(adegenet)
library(ggbiplot)
library(factoextra)
library(pca3d)
#function for write a biplot
# create_biplot <- function(pca, p1, p2, prob, group){
#   g<-ggbiplot(pca,
#               # obs.scale = 1,
#               # var.scale = 1,
#               group = group,
#               ellipse = T,
#               circle = F,
#               ellipse.prob = prob,
#               choices = c(p1, p2))
#   
#   print(g)
# }
# create_PCA_only <- function(pca, p1, p2, prob, group){
#   g<-ggbiplot(pca,
#               obs.scale = 1,
#               var.scale = 1,
#               group = group,
#               ellipse = T,
#               circle = T,
#               ellipse.prob = prob,
#               choices = c(p1, p2),
#               var.axes = FALSE
#               )+theme(legend.position = "bottom",axis.ticks=element_line(colour="black"),
#                       panel.border=element_rect(fill=NA))+ggtitle("PCA")
#   print(g)
# }
data <- read.csv("data.csv")
grand <-data
grand_pca <-prcomp(grand[,-1])

# summary(grand_pca)
str(grand_pca)
ggbiplot(grand_pca,ellipse=TRUE, groups=grand$Location
         )
fviz_pca_ind(grand_pca, geom.ind = "point", 
             pointshape = 21,
             pointsize = 2,
             fill.ind = grand$mz_rt,
             col.ind = "black",
             palette = "jco", 
             addEllipses = TRUE,
             ellipse.level = 0.95,
             invisible ="quali",
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = ""
) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom",
        panel.border=element_rect(fill=NA), axis.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20, face = "bold"))
#redo of Loading plot to check the correctness
# PCAloadings <- grand_pca$rotation
# plot(PCAloadings[,1:2],   # x and y data
#      pch=21,              # point shape
#      bg="black",          # point color
#      cex=1,               # point size
#      main="Loadings"      # title of plot
# )
# text(PCAloadings[,1:2],             # sets position of labels
#      labels=rownames(PCAloadings)   # print labels
# )
#just loading
# fviz_pca_var(grand_pca,
#              geom = c("text", "point"),
#              repel = TRUE,     # Avoid text overlapping
#              title = ""
#              
# )

#creating pca biplots
# create_biplot(grand_pca, 1, 2, 0.95, grand$lv)
# create_PCA_only(grand_pca, 1, 2, 0.95, grand$lv)
# create_biplot(grand_pca, 2, 3, 0.95, grand$lv)
# create_biplot(grand_pca, 1, 3, 0.95, grand$lv)

#pcs with individual variance explained
fviz_eig(grand_pca, addlabels = TRUE, main = "", 
         xlab = "Principle Component", ylab = "Variance Contribution Ratio (%)")

#just loading
fviz_pca_var(grand_pca,
             geom = c("text", "point"),
             repel = TRUE,     # Avoid text overlapping
             title = ""
             
)
#PCA
fviz_pca_ind(grand_pca, geom.ind = "point", 
             pointshape = 21,
             pointsize = 2,
             fill.ind = grand$lv,
             col.ind = "black",
             palette = "jco", 
             addEllipses = TRUE,
             ellipse.level = 0.95,
             invisible ="quali",
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = ""
             ) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom",
        panel.border=element_rect(fill=NA), axis.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20, face = "bold"))

#creating 3d pca plot
gr <- factor(grand$lv)
pca3d(grand_pca, group=gr,components = 1:3,
      show.ellipses = FALSE, ellipse.ci = 0.95, legend="bottomleft",
      show.plane = FALSE, title = "3D PCA plots")

# PCA using only 5 elements
grand_top4 <- subset(grand, select = c("lv", "Na", "Al", "Cd", "Rb"))
pca_grand_top4 <-prcomp(grand_top4[,-1])
fviz_pca_ind(pca_grand_top4, geom.ind = "point", 
             axes = c(1,2),
             pointshape = 21,
             pointsize = 2,
             fill.ind = grand_top4$lv,
             col.ind = "black",
             palette = "jco", 
             addEllipses = TRUE,
             ellipse.level = 0.95,
             invisible ="quali",
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "GI Rice") +
  ggtitle("PCA plot of PC 1 and PC 2 with 4 elements only") +
  theme(plot.title = element_text(hjust = 0.5))

fviz_pca_var(pca_grand_top4,
             axes = c(1,2),
             geom = c("text", "point"),
             repel = TRUE,     # Avoid text overlapping
             title = "Loading Plot"
             
)

create_biplot(pca_grand_top4, 1, 2, 0.95, grand_top4$lv)
gr <- factor(grand_top4$lv)
pca3d(pca_grand_top4, group=gr,components = 1:2,
      show.ellipses = FALSE, ellipse.ci = 0.95, legend="bottomleft",
      show.plane = FALSE, title = "3D PCA plots")
