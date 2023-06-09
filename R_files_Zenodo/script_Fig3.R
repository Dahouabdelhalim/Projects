
# Script corresponding to Fig. 3 in Beinsteiner et al., accepted in Plos Genetics, 2021

# loading the ape and phytools libraries (R version 3.2.2)
library(ape)
library(phytools)

# setting working directory. Warning: the path has to be adjusted on your own computer!
setwd( "/User/yourname/yourworkingdirectory")

# reading a tree without branch lengths
tree<-read.nexus("tree_Fig3.nex"); 

# reading the table with RxxxE presence/absence data (panel A)
traitA <- read.csv("Fig3A_RxxxE.csv", header =T, row.names = 1);
# coding: R = RxxxE motif present; Q = RxxxE motif absent; S = variable pattern

# transforming the data table into a vector format.
vec_traitA<-traitA[,1]
names(vec_traitA)<-rownames(traitA)
head(vec_traitA)

# stochastic mapping under the symmetrical model 
simmaptrees <- make.simmap(tree,vec_traitA, model="SYM", nsim= 10000, Q="empirical", pi="estimated")

# Below is the expected output:
Using pi estimated from the stationary distribution of Q assuming a flat prior.
pi =
       Q        R        S 
0.333333 0.333333 0.333333 

make.simmap is sampling character histories conditioned on the transition matrix
Q =
           Q          R          S
Q -0.1590832  0.1590832  0.0000000
R  0.1590832 -0.7056928  0.5466095
S  0.0000000  0.5466095 -0.5466095
(estimated using likelihood);
and (mean) root node prior probabilities
pi =
        Q         R         S 
0.3333333 0.3333333 0.3333333 

# end of expected output


# computing the state frequencies from the stochastic maps for each internal nodes. Posterior probabilities are illustrated by pies. 
simmap_summary<-describe.simmap(simmaptrees, plot=TRUE);
simmap_summary

10000 trees with a mapped discrete character with states:
 Q, R, S 

trees have 5.6168 changes between states on average

changes are of the following types:
        Q,R Q,S   R,Q    R,S S,Q    S,R
x->y 0.2277   0 3.126 1.6915   0 0.5716

mean total time spent in each state is:
              Q         R          S    total
raw  17.8942285 3.2822715 0.88050993 22.05701
prop  0.8112717 0.1488085 0.03991973  1.00000

simmap_summary$ace

       Q      R      S
29 0.0151 0.9738 0.0111
30 0.0257 0.9639 0.0104
31 0.0144 0.9758 0.0098
32 0.0152 0.9772 0.0076
33 0.0133 0.7713 0.2154
34 0.0005 0.3881 0.6114
35 0.0237 0.9734 0.0029
36 0.9872 0.0127 0.0001
37 0.9945 0.0052 0.0003
38 0.9997 0.0003 0.0000
39 0.9994 0.0006 0.0000
40 1.0000 0.0000 0.0000
41 0.9995 0.0005 0.0000
42 1.0000 0.0000 0.0000
43 0.9997 0.0003 0.0000
44 1.0000 0.0000 0.0000
45 0.9999 0.0001 0.0000
46 0.9996 0.0004 0.0000
47 0.9998 0.0002 0.0000
48 1.0000 0.0000 0.0000
49 1.0000 0.0000 0.0000
50 0.9999 0.0001 0.0000
51 0.9998 0.0002 0.0000
52 0.9998 0.0002 0.0000
53 1.0000 0.0000 0.0000
54 0.9999 0.0001 0.0000
55 1.0000 0.0000 0.0000

# plotting posterior density of stochastic mapping on the tree and saving it into a pdf file.
pdf("Simmap_RxxxE.pdf", height=7, width=6)
plot(simmap_summary, cex=0.6, fsize=1);
dev.off()
# note: final color for the variable state was further adjusted using Inkscape (color code: d4aa00ff)



# reading the table with dimerization data (panel B)
traitB <- read.csv("Fig3B_dimerization.csv", header =T, row.names = 1);
# coding: A = unknown; B = DR homodimer; C = IR homodimer; D = monomer; E = heterodimer 

# transforming the data table into a vector format.
vec_traitB<-traitA[,1]
names(vec_traitB)<-rownames(traitB)
head(vec_traitB)

# stochastic mapping under the symmetrical model 
simmaptrees <- make.simmap(tree,vec_traitB, model="SYM", nsim= 10000, Q="empirical", pi="estimated")

# Below is the expected output:
Using pi estimated from the stationary distribution of Q assuming a flat prior.
pi =
  A   B   C   D   E 
0.2 0.2 0.2 0.2 0.2 

make.simmap is sampling character histories conditioned on the transition matrix
Q =
           A          B         C          D          E
A -0.4964306  0.4964306  0.000000  0.0000000  0.0000000
B  0.4964306 -1.0833986  0.113973  0.2126396  0.2603555
C  0.0000000  0.1139730 -0.113973  0.0000000  0.0000000
D  0.0000000  0.2126396  0.000000 -0.4703572  0.2577176
E  0.0000000  0.2603555  0.000000  0.2577176 -0.5180730
(estimated using likelihood);
and (mean) root node prior probabilities
pi =
  A   B   C   D   E 
0.2 0.2 0.2 0.2 0.2 
# end of expected output


# computing the state frequencies from the stochastic maps for each internal nodes. Posterior probabilities are illustrated by pies. 
simmap_summary<-describe.simmap(simmaptrees, plot=TRUE);
simmap_summary

10000 trees with a mapped discrete character with states:
 A, B, C, D, E 

trees have 13.6732 changes between states on average

changes are of the following types:
        A,B A,C A,D A,E    B,A    B,C    B,D    B,E C,A    C,B C,D C,E D,A    D,B D,C    D,E E,A    E,B E,C
x->y 0.8709   0   0   0 3.3923 1.0652 1.3538 2.5725   0 0.0522   0   0   0 0.3815   0 0.5401   0 1.3382   0
        E,D
x->y 2.1065

mean total time spent in each state is:
              A         B         C          D         E    total
raw  2.00590310 6.4981531 3.2961030 1.72767631 8.5291745 22.05701
prop 0.09094175 0.2946072 0.1494356 0.07832777 0.3866877  1.00000

simmap_summary$ace
       A      B      C      D      E
29 0.1036 0.8945 0.0000 0.0005 0.0014
30 0.1261 0.8714 0.0001 0.0006 0.0018
31 0.1019 0.8964 0.0000 0.0004 0.0013
32 0.0097 0.9833 0.0000 0.0005 0.0065
33 0.0043 0.9946 0.0000 0.0005 0.0006
34 0.0036 0.9955 0.0000 0.0003 0.0006
35 0.0040 0.9689 0.0000 0.0023 0.0248
36 0.0016 0.9684 0.0002 0.0109 0.0189
37 0.0058 0.8731 0.0001 0.0900 0.0310
38 0.0000 0.0223 0.0000 0.9704 0.0073
39 0.0021 0.9723 0.0045 0.0026 0.0185
40 0.0002 0.0158 0.9840 0.0000 0.0000
41 0.0000 0.0023 0.9977 0.0000 0.0000
42 0.0000 0.0001 0.9999 0.0000 0.0000
43 0.0102 0.9345 0.0005 0.0044 0.0504
44 0.0012 0.1533 0.0000 0.0295 0.8160
45 0.0000 0.0200 0.0000 0.0106 0.9694
46 0.0000 0.0925 0.0000 0.1336 0.7739
47 0.0007 0.1616 0.0000 0.1010 0.7367
48 0.0000 0.0082 0.0000 0.0053 0.9865
49 0.0000 0.0033 0.0000 0.0032 0.9935
50 0.0000 0.0004 0.0000 0.0010 0.9986
51 0.0000 0.0005 0.0000 0.0014 0.9981
52 0.0000 0.0007 0.0000 0.0006 0.9987
53 0.0000 0.0002 0.0000 0.0005 0.9993
54 0.0000 0.0004 0.0000 0.0010 0.9986
55 0.0000 0.0120 0.0000 0.0897 0.8983

# plotting posterior density of stochastic mapping on the tree and saving it into a pdf file.
pdf("Simmap_binding.pdf", height=7, width=6)
plot(simmap_summary, cex=0.6, fsize=1);
dev.off()
# note: final color for the unknown state was further adjusted using Inkscape version 0.91 (color code: 808080ff)
