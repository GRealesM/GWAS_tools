# Compare Betas pairwise

# Objective: This script aims to compare files of the same trait pairwise and compute if their BETAS are pointing in the same direction, thus allowing us to easily check possible REF/ALT misassignments

traits <- sapply(strsplit(dir("Filtered_datasets/"),"_"), `[`, 1)
traits <- as.data.frame(table(traits))
traits <- traits[traits$Freq > 1,]
Main <- data.frame()

for (i in 1:length(traits$traits)){
  list.of.files <- list()
  files <- dir("Filtered_datasets/", pattern = paste("^", traits$traits[i], "_", sep = ""))
  for (j in 1:length(files)){
    x <- read.table(paste("Filtered_datasets/",files[j], sep = ""), header = T, sep = "\t")
    x <- x[, c("pid38", "BETA")]
    names(x)[2] <- paste(names(x)[2], files[j], sep="_")
    list.of.files[[j]] <- x
  }
  df.of.files <- Reduce(function(x, y) merge(x, y, by = "pid38", all = TRUE), list.of.files)
  dfcopy <- df.of.files
  assign(paste("trait", as.character(traits$traits[i]), sep = "_"), dfcopy)
  df.of.files <- df.of.files[, -1]
  Element1 <- names(df.of.files)[combn(length(df.of.files), 2)[1,]]
  Element2 <- names(df.of.files)[combn(length(df.of.files), 2)[2,]]
  correlation <-combn(length(df.of.files), 2, function(x) cor(df.of.files[,x[1]], df.of.files[,x[2]],use = "complete.obs"))
  newdf <- data.frame(Element1 = Element1, Element2 = Element2, correlation = correlation)
  Main <- rbind(Main, newdf)
}




