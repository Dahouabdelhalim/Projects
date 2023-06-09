#install.packages("~/scripts/related_1.0.tar.gz", repos = NULL, type = "source")
library(related)

input <- readgenotypedata(file.choose())

output <- coancestry(input$gdata, trioml = 2, error.rates = 0.001, ci95.num.bootstrap = 500)

write.csv(output$relatedness, file="trioml_relatedresults.csv")
write.csv(output$relatedness.ci95, file="trioml_relatedCIs.csv")

output$relatedness[output$relatedness$trioml>0.125,] #putative cousins
output$relatedness[output$relatedness$trioml>0.25,] #putative half-siblings/PO
output$relatedness[output$relatedness$trioml>0.5,] #putative full siblings

