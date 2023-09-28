setwd("./")

rm(list = ls())

list.files()
files_ind <- list.files(pattern="_2_.tsv")

##data frame creation

for (i in 1:length(files_ind)) assign(files_ind[i], read.table(files_ind[i],row.names =
1,check.names = F,header = T))

dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))

data_table <- merge(dfs[1],dfs[2],by =0,all =T)
rownames(data_table) <- data_table$Row.names
data_table$Row.names <- NULL

data_table <-merge(data_table, dfs[3],by =0, all= T)
rownames(data_table) <- data_table$Row.names
data_table$Row.names <- NULL

for (i in 4:length(dfs))
{
  data_table <-merge(data_table, dfs[[i]],by =0, all= T)
  rownames(data_table) <- data_table$Row.names
  data_table$Row.names <- NULL
}

data_table[,4] <- NULL

colnames(data_table) <- c('M_DEPTH_A','M_DEPTH_X','F_MISS_X','F_MISS_Y')

##computing of X and Y indexes.

index_x <- data_table[2]/data_table[,1]
colnames(index_x)[1] <- "Index_X"

index_y <- ((1-data_table[3])-(1-data_table[,4]))/(1-data_table[,3])
colnames(index_y)[1] <- "Index_Y"

data_table_2 <- merge(data_table, index_x, by=0, all=TRUE)
rownames(data_table_2) <- data_table_2[,1]
data_table_2[,1] <- NULL

data_table_2 <- merge(data_table_2, index_y, by=0, all=TRUE)
rownames(data_table_2) <- data_table_2[,1]
data_table_2[,1] <- NULL

data_table_2$Sex[data_table_2$Index_Y <= 0.1 ] <- "male"
data_table_2$Sex[data_table_2$Index_Y > 0.1 ] <- "female"

##ploting

pdf("sexing_plots.pdf")
plot(x=data_table_2$Index_X, y=data_table_2$Index_Y, col=factor(data_table_2$Sex), xlab="Index_X", ylab="Index_Y", pch=19)
abline(h=0.1, col="blue")

plot(x=data_table_2$M_DEPTH_A, y=data_table_2$Index_X, col=factor(data_table_2$Sex), xlab="Depth", ylab="Index_X", pch=19)
abline(h=1, col="black")

plot(x=data_table_2$M_DEPTH_A, y=data_table_2$Index_Y, col=factor(data_table_2$Sex), xlab="Depth", ylab="Index_Y", pch=19)
abline(h=0, col="black")
dev.off()

write.csv (data_table_2, "final_sexing.csv")
