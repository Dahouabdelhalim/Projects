#####Mixed models####
allom_t <- read.csv("BeeColonyScalingData.csv")

lmer_allom_ML_full<-lmer(log_length~log_mass*hive*body_part+(1|temp/file_name),REML=F,
                         na.action = na.omit,data = allom_t)
lmer_allom_REML_full<-lmer(log10(length)~log_mass*hive*body_part+(1|temp/file_name),REML=T,
                           na.action = na.omit,data = allom_t)
anova(lmer_allom_REML_full)
anova(lmer_allom_ML_full)
step_res<-step(lmer_allom_REML_full)
final<-get_model(step_res)
anova(final)
anova(lmer_allom_REML_full,refit=F,ddf="Satterthwaite")

#Comparing with LRT
lmer_allom_ML_red<-lmer(log10(length)~(log_mass+body_part+hive)^2+(1|temp/file_name),REML=F,
                        na.action = na.omit,data = allom_t)
anova(lmer_allom_ML_full,lmer_allom_ML_red)

#Posthoc tests
comps<-emtrends(lmer_allom_REML_full, pairwise ~ body_part*hive, var = "log_mass")
comps2<-emtrends(lmer_allom_REML_full, pairwise ~ body_part|hive, var = "log_mass")
comps3<-emtrends(lmer_allom_REML_full, pairwise ~ body_part, var = "log_mass")
comp_means1<-emmeans(lmer_allom_REML_full,pairwise ~ body_part)
comp_means2<-emmeans(lmer_allom_REML_full,pairwise ~ body_part*hive)

#####FAMD####
ClustData <- read.csv("BeeColonyClusterData.csv")

allom_famd <- FAMD(ClustData)
summary(allom_famd)
famd_info<-get_famd(allom_famd,element = "var")
famd_info$contrib
facto_summarize(allom_famd, "var", axes = 1:2)

#visualise contribution to dimensions:
fviz_screeplot(allom_famd)
fviz_famd_var(allom_famd, "var",repel = TRUE)
fviz_famd_ind(allom_famd, repel = TRUE)

fviz_contrib(allom_famd, choice = "var", axes = 1)
fviz_contrib(allom_famd, choice = "var", axes = 2)
fviz_contrib(allom_famd, choice = "var", axes = 3)

#Cluster affiliations:
clust.hcpc<-HCPC(allom_famd, nb.clust = -1, graph = F)#nb.clust = -1 for HCPC to automatically select cluster number
clust.hcpc$data.clust

