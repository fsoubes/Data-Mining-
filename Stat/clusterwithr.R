# @ authhor : Debras Guillamaury


library(factoextra)
library(NbClust)

database = read.csv(file = "realtable.csv", sep = ";", dec=".",header = TRUE,nrows=1000)
row.names(database) <- database$id_protein
database <- na.omit(database)
database <- database[, -c(1)]
View(database)
database = scale(database)
head(database)

nb <- NbClust(database, distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "kmeans")
fviz_nbclust(nb)







fviz_nbclust(database, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

