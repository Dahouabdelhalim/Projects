library(statnet)
library(xergm)
library(parallel)

# read in data

law_oc_mat <- as.matrix(read.csv("law_occurence.csv", row.names = 1)) # load occurence of issues in laws
proj_oc_mat <- as.matrix(read.csv("issue_occurence.csv", row.names = 1)) # load occurence of issues in actor activity
attributes <- read.csv("actor_attributes.csv", stringsAsFactors = FALSE)

# Issue association strength based on law overlap

# ochiai similiarity

issue_issue_law_mat <- t(law_oc_mat) %*% law_oc_mat

colsums_x_rowsums_mat_docs <- outer(colSums(issue_issue_law_mat), rowSums(issue_issue_law_mat))
issue_issue_law_mat_norm <- (issue_issue_law_mat / sqrt(colsums_x_rowsums_mat_docs))
diag(issue_issue_law_mat_norm) <- 0 # remove self-ties
issue_issue_lawsim_ochiai <- issue_issue_law_mat_norm

# twomode actor-issue network

#order in same way as issue-issue sim
issue_actor <- t(proj_oc_mat[,colnames(issue_issue_law_mat_norm)])

# check order and write
if(!(identical(rownames(issue_actor), colnames(issue_issue_law_mat_norm)))){
  print("wrong order!")
}

issue_actor_nw <- network(issue_actor, bipartite = T, directed = F)

# assign actor attribute values

set.vertex.attribute(issue_actor_nw, attrname = "mode",
                     value = ifelse(get.vertex.attribute(issue_actor_nw,attrname = "vertex.names") %in% colnames(issue_actor),
                                    "actor","issue"))
set.vertex.attribute(issue_actor_nw, attrname = "type",
                     value = ifelse(get.vertex.attribute(issue_actor_nw,attrname = "vertex.names") %in% colnames(issue_actor),
                                    attributes$type_reduced,NA))
set.vertex.attribute(issue_actor_nw, attrname = "level",
                     value = ifelse(get.vertex.attribute(issue_actor_nw,attrname = "vertex.names") %in% colnames(issue_actor),
                                    attributes$level,NA))

nw_atts_names <- list.vertex.attributes(issue_actor_nw)
nw_atts_df <- data.frame(do.call(cbind,lapply(nw_atts_names[!(nw_atts_names %in% c("na"))],
                                              function(att) get.vertex.attribute(issue_actor_nw, attrname = att))))
colnames(nw_atts_df) <- nw_atts_names[!(nw_atts_names %in% c("na"))]

# create homophily term

# define two-mode homophily function; N = two-mode network; X = covariate input
homophily <- function(N, X, f) {
  mat <- matrix(0, nrow = nrow(N), ncol = ncol(N))
  for (i in 1:nrow(N)) {
    for (j in 1:ncol(N)) {
      for (k in 1:nrow(N)) {
        if (i != k) {
          mat[i, j] <- mat[i, j] + (N[k, j] * f(i, k, X))
        }
      }
    }
  }
  return(mat)
}

# pass issue-issue sim matrix via this
f_lawsim <- function(i, k, X) {
  return(X[i, k])
}

lawsim <- homophily(issue_actor, issue_issue_lawsim_ochiai, f_lawsim)

# ergm setup

np = detectCores() - 2

mcsize <- 1000 #set to 2000 back later

# model

model1 <- ergm(
  issue_actor_nw
  ~ edges                      # endogenous control
  + b2degree(1)                # endogenous control for single-issue actors hist(colSums(issue_actor))
  + b2star(2)                  # endogenous control for clustering around some actors
  # + gwdegree(1, fixed = TRUE)  # endogenous control for degree distribution
  + gwnsp(1, fixed = TRUE)     # endogenous control for general clustering
  + edgecov(lawsim)        # law-induced clustering
  + b2factor("type")
  + b2nodematch("type")
  , control = control.ergm(MCMC.samplesize = mcsize, MCMC.interval = 2000,
                           parallel=np, parallel.type="PSOCK")
)
summary(model1)


# coefficient plot

results <- summary(model1)

res <- as.data.frame(cbind(results$coefs, btergm::confint(model1)))

res <-
  rbind(rep(NA,nrow(res)),
        rep(NA,nrow(res)),
        res[1:which(grepl(pattern = "gwnsp", x = rownames(res))),],
        rep(NA,nrow(res)),
        rep(NA,nrow(res)),
        res[-(1:which(grepl(pattern = "edgecov", x = rownames(res)))),],
        rep(NA,nrow(res)),
        rep(NA,nrow(res)),
        res[which(grepl(pattern = "edgecov", x = rownames(res))):
              which(grepl(pattern = "edgecov", x = rownames(res))),])
res

res$labels <- c("Endogeneous controls"," ",
                "Edges","Actor degree: 1", "Two-stars (centred on actors)",
                "Non-edgewise shared partners (fixed at 1)",
                "Exogeneous controls","  ",
                "Local administration", "Other actors", "Politics",
                "Private sector","Science","Service providers","State and national administration",
                "Actor type homophily",
                "Law-driven integration effect","   ",
                "Issue association in law framework")
res$order <- c(nrow(res):1)#c(1:4,7:9,5:6)

coefficienct_plot_model1 <- ggplot(res, aes(res$Estimate, reorder(labels, order))) +
  geom_errorbarh(aes(xmin=res$`2.5 %`, xmax=res$`97.5 %`),
                 size=2, height=.2, color="darkgrey") +
  geom_point(aes(size = 10)) +
  geom_vline(xintercept=0, linetype=2) + geom_point() +
  labs(title="",
       y="Coefficients",
       x="Value") +
  theme_minimal(25) +
  theme(axis.title.x = element_text(hjust=.5, size = 12),
        axis.title.y = element_text(hjust=.6, size = 12),
        plot.title = element_text(hjust = 0, size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12)) +
  guides(size = F)
coefficienct_plot_model1

# degeneracy check: model 1
mcdiag1 <- mcmc.diagnostics(model1)

# endogenous fit: model 1
gof1 <- gof(model1, statistics = c(b1star, b2star, deg, nsp,
                                   geodesic, rocpr), rocprgof = TRUE, nsim = 500, #500 for proper
            control = control.gof.ergm(parallel=np, parallel.type="PSOCK"))
plot(gof1, roc.rgraph = TRUE, pr.rgraph = TRUE)

# edge prediction: model 1
rocpr1 <- gof(model1, statistics = rocpr, nsim = 500,
              control = control.gof.ergm(parallel=np, parallel.type="PSOCK")) #1000 for proper
plot(rocpr1, roc.rgraph = TRUE, pr.rgraph = TRUE)
