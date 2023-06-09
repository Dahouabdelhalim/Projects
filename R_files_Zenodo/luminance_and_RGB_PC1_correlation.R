#Read in file containing luminance values and file containing body region PC scores
luminance=read.csv("hair_luminance.csv")

RGB_torso_scores=read.csv("RGB_traits_torso_scores.csv")
#OR
RGB_fore_scores=read.csv("RGB_traits_fore_scores.csv")
#OR
RGB_head_scores=read.csv("RGB_traits_head_scores.csv")
#OR
RGB_hind_scores=read.csv("RGB_traits_hind_scores.csv")
#OR
RGB_tail_scores=read.csv("RGB_traits_tail_scores.csv")

#format them properly
rownames(luminance)=luminance$species

rownames(RGB_torso_scores)=RGB_torso_scores$X
#OR
rownames(RGB_fore_scores)=RGB_fore_scores$X
#OR
rownames(RGB_head_scores)=RGB_head_scores$X
#OR
rownames(RGB_hind_scores)=RGB_hind_scores$X
#OR
rownames(RGB_tail_scores)=RGB_tail_scores$X

#organize luminance data in order of RGB_fore_scores data
luminance_order <- luminance[match(RGB_torso_scores$X,rownames(luminance)),]
#OR
luminance_order <- luminance[match(RGB_fore_scores$X,rownames(luminance)),]
#OR
luminance_order <- luminance[match(RGB_head_scores$X,rownames(luminance)),]
#OR
luminance_order <- luminance[match(RGB_hind_scores$X,rownames(luminance)),]
#OR
luminance_order <- luminance[match(RGB_tail_scores$X,rownames(luminance)),]

##run correlation test
luminance_torso_vector=luminance_order[,c(5:6,18:19)]
torso_luminance_RGB_vector<- cbind(RGB_torso_scores[,8],luminance_torso_vector)
colnames(torso_luminance_RGB_vector)[1] <- "PC1"
cor.test(torso_luminance_RGB_vector$PC1, torso_luminance_RGB_vector$m.chest.L, method = c("spearman"), exact = FALSE)
#OR
luminance_fore_vector=luminance_order[,c(7:9,16:17)]
fore_luminance_RGB_vector<- cbind(RGB_fore_scores[,8],luminance_fore_vector)
colnames(fore_luminance_RGB_vector)[1] <- "PC1"
cor.test(fore_luminance_RGB_vector$PC1, fore_luminance_RGB_vector$ven.low.forelimb.L, method = c("spearman"), exact = FALSE)
#OR
luminance_head_vector=luminance_order[,c(3:4,15)]
head_luminance_RGB_vector<- cbind(RGB_head_scores[,8],luminance_head_vector)
colnames(head_luminance_RGB_vector)[1] <- "PC1"
cor.test(head_luminance_RGB_vector$PC1, head_luminance_RGB_vector$ven.neck.L, method = c("spearman"), exact = FALSE)
#OR
luminance_hind_vector=luminance_order[,c(10:12,20:21)]
hind_luminance_RGB_vector<- cbind(RGB_hind_scores[,8],luminance_hind_vector)
colnames(hind_luminance_RGB_vector)[1] <- "PC1"
cor.test(hind_luminance_RGB_vector$PC1, hind_luminance_RGB_vector$ven.low.hindlimb.L, method = c("spearman"), exact = FALSE)
#OR
luminance_tail_vector=luminance_order[,c(13:14,22)]
tail_luminance_RGB_vector<- cbind(RGB_tail_scores[,8],luminance_tail_vector)
colnames(tail_luminance_RGB_vector)[1] <- "PC1"
cor.test(tail_luminance_RGB_vector$PC1, tail_luminance_RGB_vector$ven.dist.tail.L, method = c("spearman"), exact = FALSE)


