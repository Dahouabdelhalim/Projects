### Phylogenetic contrast for the Ammophila-Paraxenos project (A Double-Edged Sword: Parental care increases risk of offspring infection by a maternally-vectored parasite)
## Authors: Rebecca Jean Millena, Jay Rosenheim

## Using "ape" to load the tree from Field et al. 2011, using just the Newick tree from the .nexorg file (names are edited to match the tip label names to my data names)
library(ape)
ammophila_tree <- read.tree(text = "(Lasioglossum:0.200882,Apis:0.228084,(Hesperapis:0.210663,(Macropis:0.122921,((Anacrabro:0.254768,Tachysphex:0.218863):0.080965,((((((((((((((aberti:0.023605,((azteca:0.003716,silvestris:0.005236):0.025090,breviceps:0.026590):0.015199):0.006636,((dysmica:0.028393,kennedyi:0.021044):0.007019,(gracilis:0.038645,urnaria:0.029724):0.012881):0.013958):0.005826,(mediata:0.026353,(sabulosa:0.025416,strenua:0.025737):0.008225):0.013754):0.006499,(((ferrugineipes:0.030917,poecilocnemis:0.040361):0.005554,heydeni:0.027392):0.002669,laevicollis:0.037536):0.009747):0.013401,((femurrubra:0.024520,pictipennis:0.026916):0.020993,(fernaldi:0.051831,nigricans:0.026111):0.019499):0.019642):0.003999,pubescens:0.050420):0.009108,((((marshi:0.030364,stangei:0.033871):0.010587,wrightii:0.039319):0.006451,procera:0.039236):0.012789,vulcania:0.067664):0.006470):0.006310,(((Eremnophila_aureonotata:0.017420,Eremnophila_binodis:0.022032):0.031733,(Hoplarmata:0.024032,Hoplclypeata:0.030850):0.030640):0.011383,(Podalonia_affinis:0.041041,(((Podalonia_canescens:0.044766,(Podalonia_luffii:0.022506,Podalonia_mauritanica:0.042731):0.016642):0.010582,((Podalonia_melaena:0.030018,Podalonia_pubescens:0.042306):0.008928,Podalonia_valida:0.036910):0.011473):0.005385,Podalonia_hirsuta:0.033723):0.015872):0.003067):0.005181):0.009358,(Parapsherero:0.043806,Parapsturanica:0.051906):0.019279):0.017127,Eremochares_dives:0.054061):0.040428,Prionyx_kirbii:0.106347):0.026420,(Isodontia_mexicana:0.068969,Sphex_ichneumoneus:0.093587):0.055576):0.028526,Stangeella_cyaniventris:0.080489):0.031672,(((Chalybion_zimmermanni_aztecum:0.080428,(Penepodium_sp.:0.068750,Podium_rufipes:0.041077):0.029841):0.015119,Sceliphron_caementarium:0.100646):0.014302,Chlorion_aerarium:0.070248):0.031630):0.091007):0.137451):0.048223):0.052700);")

## Step after this is to prune the tree to just the 10 taxa in question, then create a CSV of the tip.labels and the two traits I'm comparing (provisioning data and parasitism rates)
plot(ammophila_tree)
## Interactively removing unnecessary taxa from the tree by clicking on them
stylotree <- drop.tip(ammophila_tree, interactive = TRUE)
plot(stylotree)
## One more round of removing unnecessary taxa
stylotree <- drop.tip(stylotree, interactive = TRUE)
plot(stylotree, font = 2, edge.width = 2, align.tip.label = TRUE)
## Save the tree with export

## Here I'm essentially just following the outline of "https://www.r-phylo.org/wiki/HowTo/Phylogenetic_Independent_Contrasts" to familiarize myself with it
### Creating a numeric vector for the traits of interest (provisioning data and parasitism rates) for use in the contrast, named "styloptree.txt"
stylodata <- read.table("styloptree.txt")
prov <- stylodata$mean_prey
para <- stylodata$percent_parasitized
## Associating the data with the correct tips
names(prov) <- row.names(stylodata)
names(para) <- row.names(stylodata)

### Calculate contrasts that have been scaled using the expected variances
contrastprov <- pic(prov, stylotree)
contrastpara <- pic(para, stylotree)

### Extract expected variances and contrasts
contrastprov.var <- pic(prov, stylotree, var.contrasts = TRUE)
contrastpara.var <- pic(para, stylotree, var.contrasts = TRUE)

### Calculating contrasts for both variables
contrasts_stylodata <- apply(stylodata, 2, pic, stylotree)

### Displaying contrasts
contrastprov.var
contrastprov
contrastpara.var
contrastpara

### Plotting the contrasts at the appropriate nodes with variances
plot(stylotree, font = 2, edge.width = 2, align.tip.label = TRUE)
nodelabels(round(contrastprov.var [,1],3), adj = c(0, -0.5), frame = "n", col = "black", family = "Times")
nodelabels(round(contrastpara.var [,1],3), adj = c(0, 1), frame = "n", col = "red", family = "Times")

# Here, at the end, is the most important part
## Correlated trait evolution; running a linear regression on the contrasts
regressprovpara <- lm(contrastpara~contrastprov)
summary.lm(regressprovpara)

## Plot the contrasts
plot(contrastprov, contrastpara)
abline(regressprovpara)

## Plotting a prettier image with ggplot2
library(ggplot2)
library(ggpubr)
ggplot(mapping = aes(x = contrastprov, y = contrastpara)) +
  geom_jitter(alpha = 0.8) +
  ## method = 'lm' means that I am fitting a linear regression model; shaded area is the 95% confidence interval that signifies there is 95% confidence that the true regression line lies within the shaded region
  ## formula = y ~ x - 1 forces the regression through the origin
  geom_smooth(method ='lm',color = "black", formula = y ~ x - 1) +
  stat_regline_equation(formula = y ~ x - 1, label.x = -6.3, label.y = 25, family = "Times") +
  xlab("Provisioning Rate Contrast") +
  ylab("Parasitism Rate Contrast") +
  theme_gray(base_family = "Times") +
  stat_cor(color = "black", family = "Times")
  
ggsave("../ammophila_strepsiptera/phylocontrast_regression_0_paper.pdf", dpi = 300)
