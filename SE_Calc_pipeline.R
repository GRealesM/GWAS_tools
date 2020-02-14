## SE EXPLORER ####
## Description: The goal of this script is to calculate SE when missing, using BETA and P.

#install.packages("data.table")
library(data.table)

input <- fread("tmp_schecked.tsv")
input$Z <- sign(input$BETA) * abs(qnorm(input$P/2))
input$SE <- input$BETA/input$Z
fwrite(input, "tmp.tsv", quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")


