x_full <- new("genlight", t(Phenig$G))

indNames(x_full)<- Phenig$id
locNames(x_full)<- Phenig$SNP
pop(x_full)<- Phenig$pop
popNames(x_full)
grp <- find.clusters(x_full, max.n.clust=50, criterion = "diffNgroup")

#Selected the number of clusters to 2 and using 270 PCs

table.value(table(pop(x_full), grp$grp), col.lab=paste("inf", 1:30),row.lab=paste("ori", 1:6))
dapc1 <- dapc(x_full, pop = grp$grp, n.pca = 30, n.da = 2)
optim.a.score(dapc1)
dapc1 <- dapc(x_full, pop = grp$grp, n.pca = 13)
scatter.dapc(dapc1)
grp$grp
dapc1$loadings
