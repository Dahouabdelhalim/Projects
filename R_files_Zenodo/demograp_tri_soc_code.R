data_tri_soc<-read.csv("~/Downloads/demography_sse.csv")
			attach(data_tri_soc)
			names(data_tri_soc)
			
			fQ1<-factor(Q1)
			levels(fQ1)

##########################
#########################
#FIGURE 1, VENN DIAGRAM
##########################
##########################
pdf(file="Figure1_VD_SSE.pdf", height=6, width=6)
					vp <- venn.diagram(list(), 
 					                   fill = c(2,5,7), alpha = 0.3, filename = NULL);
					grid.draw(vp)
dev.off()

#COMPOSITION OF THE SURVEY

table(data_tri_soc $Q6)

                                                         Academic professional (staff / administration)                   Graduate student (Masters or PhD) 
                                                  2                                                  19                                                 326 
Non-academic professional (publisher, vendor, etc.)                            Non-tenure track faculty                                               Other 
                                                 10                                                  16                                                   7 
                                            Postdoc                                  Pre-tenure faculty                                     Tenured faculty 
                                                169                                                 109                                                 160 
                                      Undergraduate 
                                                 34 

##########################
##########################
#FIGURE 2, Phd year      #
##########################
##########################
			
library(gplots)
phdyear<-read.csv("~/Downloads/phd_degree.csv")
		names(phdyear)
				attach(Q7)
				
Q2_f<-factor(phdyear$Q2)

pdf(file="Figure2_SSE.pdf", height=5, width=10)

				par(mfrow=c(1,2))
				males<-(subset(phdyear, Q2_f =='Male'))
				hist(males $Q7, xlim=c(1950, 2020), ylim=c(0,100), breaks=10, xlab="Year of Ph.D. completion", main="A. Men, \\n all evolution societies")
				error_males <- qt(0.975,df=220-1)*(sd(males $Q7))/sqrt(220)
				abline(v=mean(males $Q7))
				abline(v=mean(males $Q7)+ error_males, lty=2)
				abline(v=mean(males $Q7)- error_males, lty=2)

				women<-(subset(phdyear, Q2_f =='Female'))
				hist(women $Q7, xlim=c(1950, 2020), ylim=c(0,100),  breaks=10, xlab="Year of Ph.D. completion", main="B. Women, \\n all evolution societies")
				error_ women <- qt(0.975,df=234-1)*(sd(women $Q7))/sqrt(234)
				abline(v=mean(women $Q7))
				abline(v=mean(women $Q7)+ error_, lty=2)
				abline(v=mean(women $Q7)- error_, lty=2)


dev.off()

##########################
##########################
# Question 1: Gender
##########################
##########################

#################
#TABLE 1.       #
#################
table(data_tri_soc $Q1, data_tri_soc $Q2)
                                                                                                                       
                                                                                                                            Choose not to answer Female Male
                                                                                                                          1                    0      0    2
  American Society of Naturalists (ASN)                                                                                   0                    2     33   22
  American Society of Naturalists (ASN),Society for the Study of Evolution (SSE)                                          0                    0     90   68
  American Society of Naturalists (ASN),Society of Systematic Biologists (SSB)                                            0                    0      3    2
  American Society of Naturalists (ASN),Society of Systematic Biologists (SSB),Society for the Study of Evolution (SSE)   0                    1     22   19
  None                                                                                                                    0                    1     28   19
  Society for the Study of Evolution (SSE)                                                                                0                    8    222  160
  Society of Systematic Biologists (SSB)                                                                                  0                    0     30   29
  Society of Systematic Biologists (SSB),Society for the Study of Evolution (SSE)                                         0                    0     41   37


#TABLE S1#
##########
##################################
###TABLE S1. SOCIETIES VS CENSUS
##################################

CENSUS TOTAL: 328,239,523
Female persons, percent

50.8%
Total females: 166745678
ASN
			#allowing for members in two or more societies
						females  <- c( (148), (166745678))
						totals <- c( 259, 328239523)
						prop.test(females, totals)

			#NOT allowing for members in two or more societies
						females  <- c( (33), (166745678))
						totals <- c( 55, 328239523)
						prop.test(females, totals)

SSB
			#allowing for members in two or more societies
						females  <- c( (96), (166745678))
						totals <- c( 183, 328239523)
						prop.test(females, totals)

			#NOT allowing for members in two or more societies
						females  <- c( (30), (166745678))
						totals <- c( 59, 328239523)
						prop.test(females, totals)


SSE
			#allowing for members in two or more societies
						females  <- c( (373), (166745678))
						totals <- c( 659, 328239523)
						prop.test(females, totals)

			#NOT allowing for members in two or more societies
						females  <- c( (160), (166745678))
						totals <- c( 382, 328239523)
						prop.test(females, totals)

						X-squared = 11.794, df = 1, p-value = 0.0005943

Total
			#allowing for members in two or more societies
						females  <- c( (469), (166745678))
						totals <- c( 825, 328239523)
						prop.test(females, totals)

#########################################
#total of the societies vs phds awarded in LIFE SCIENCES (2019)
#########################################

ASN
			#allowing for members in two or more societies
						females  <- c( (148), (4501))
						totals <- c( 259, 8702)
						prop.test(females, totals)

						X-squared = 2.7456, df = 1, p-value = 0.09752

			#NOT allowing for members in two or more societies
						females  <- c( (33), (4501))
						totals <- c( 55, 8702)
						prop.test(females, totals)
						
						X-squared = 1.1862, df = 1, p-value = 0.2761

SSB	
			#allowing for members in two or more societies
						females  <- c( (96), (4501))
						totals <- c( 183, 8702)
						prop.test(females, totals)
						X-squared = 0.014945, df = 1, p-value = 0.9027

			#NOT allowing for members in two or more societies
						females  <- c( (30), (4501))
						totals <- c( 59, 8702)
						prop.test(females, totals)

						X-squared = 1.2502e-05, df = 1, p-value = 0.9972

SSE
			#allowing for members in two or more societies
						females  <- c( (373), (4501))
						totals <- c( 659, 8702)

						prop.test(females, totals)

						X-squared = 5.645, df = 1, p-value = 0.01751

			#NOT allowing for members in two or more societies
						females  <- c( (160), (4501))
						totals <- c( 382, 8702)
						prop.test(females, totals)

						X-squared = 13.788, df = 1, p-value = 0.0002046



TOTAL
						female  <- c(469, 4501)
						totals <- c( (825), (8702))
						prop.test(female, totals)

						X-squared = 7.7271, df = 1, p-value = 0.00544

#########################################
#total of the societies vs phds awarded in ev bio (https://ncses.nsf.gov/pubs/nsf21308/data-tables TABLE 16)
#########################################

Total females: 166745678
ASN
		#allowing for members in two or more societies
						females  <- c( (148), (112))
						totals <- c( 259, 238)
						prop.test(females, totals)

						X-squared = 4.6595, df = 1, p-value = 0.03088

		#NOT allowing for members in two or more societies
						females  <- c( (33), (112))
						totals <- c( 55, 238)
						prop.test(females, totals)

						X-squared = 2.4978, df = 1, p-value = 0.114

SSB	
		#allowing for members in two or more societies
						females  <- c( (96), (112))
						totals <- c( 183, 238)
						prop.test(females, totals)

						X-squared = 1.0006, df = 1, p-value = 0.3172
			
		#NOT allowing for members in two or more societies
						females  <- c( (30), (112))
						totals <- c( 59, 238)
						prop.test(females, totals)

						X-squared = 0.14133, df = 1, p-value = 0.707

SSE
		#allowing for members in two or more societies
						females  <- c( (373), (112))
						totals <- c( 659, 238)
						prop.test(females, totals)

						X-squared = 6.0322, df = 1, p-value = 0.01405

		#NOT allowing for members in two or more societies
						females  <- c( (160), (112))
						totals <- c( 382, 238)
						prop.test(females, totals)

						X-squared = 1.391, df = 1, p-value = 0.2382



TOTAL
						female  <- c(469, 112)
						totals <- c( (825), (112+126))
						prop.test(female, totals)

						2-sample test for equality of proportions with continuity correction

						data:  female out of totals
						X-squared = 5.4722, df = 1, p-value = 0.01932
						alternative hypothesis: two.sided
						95 percent confidence interval:
						0.01367209 0.16187030
						sample estimates:
						   prop 1    prop 2 
						0.5671100 0.4793388 

#################
#TABLE S2       #
#################

table(data_tri_soc $Q2, data_tri_soc $Q6)
                                        Choose not to answer Female Male 
                                                       0      0    1                                   
  Academic professional (staff / administration)       0     11    8                                   
  Graduate student (Masters or PhD)                    2    201  114                                   
  Non-academic professional (publisher, vendor, etc.)  0      3    7                                   
  Non-tenure track faculty                             1      7    8                                   
  Other                                                0      5    2                                   
  Postdoc                                              5     87   75                                   
  Pre-tenure faculty                                   2     64   43                                   
  Tenured faculty                                      2     71   87                                   
  Undergraduate                                        0     20   13                                   

###############
# TABLE S3.   #
###############
                                                                                                                           Choose not to answer Female Male
                                                                                                                          1                    0      0    2
  American Society of Naturalists (ASN)                                                                                   0                    2     33   22
  American Society of Naturalists (ASN),Society for the Study of Evolution (SSE)                                          0                    0     90   68
  American Society of Naturalists (ASN),Society of Systematic Biologists (SSB)                                            0                    0      3    2
  American Society of Naturalists (ASN),Society of Systematic Biologists (SSB),Society for the Study of Evolution (SSE)   0                    1     22   19
  None                                                                                                                    0                    1     28   19
  Society for the Study of Evolution (SSE)                                                                                0                    8    222  160
  Society of Systematic Biologists (SSB)                                                                                  0                    0     30   29
  Society of Systematic Biologists (SSB),Society for the Study of Evolution (SSE)                                         0                    0     41   37


#################################
#Professional stage GENDER, TABLE S3
#################################

> table(SSE_data$Q6, SSE_data$Q2)
                                                     
                                                      Choose not to answer Female Male
  Academic professional (staff / administration)                         0      8    4                                            
  Graduate student (Masters or PhD)                                      2    143   85                                           
  Non-academic professional (publisher, vendor, etc.)                    0      0    4                                            
  Non-tenure track faculty                                               1      7    5                                            
  Other                                                                  0      3    1                                            
  Postdoc                                                                2     73   57                                            
  Pre-tenure faculty                                                     2     50   29                                            
  Tenured faculty                                                        1     59   71                                            
  Undergraduate                                                          0     10    9                                            

#grad vs postdoc
								female  <- c(201, 87)
								totals <- c( (114+201), (87+75))
								prop.test(female, totals)

								2-sample test for equality of proportions with continuity correction

								data:  female out of totals
								X-squared = 4.1544, df = 1, p-value = 0.04153
								alternative hypothesis: two.sided
								95 percent confidence interval:
								0.003047124 0.199069278
								 estimates:
								   prop 1    prop 2 
								0.6380952 0.5370370 

#grad vs. untenured
								female  <- c(201, 64)
								totals <- c( (114+201), (64+43))
								prop.test(female, totals)

								>  prop.test(female, totals)

								2-sample test for equality of proportions with continuity correction

								data:  female out of totals
								X-squared = 0.38836, df = 1, p-value = 0.5332
								alternative hypothesis: two.sided
								95 percent confidence interval:
								 -0.07328115  0.15320994
								sample estimates:
								   prop 1    prop 2 
								0.6380952 0.5981308 

#grad vs. tenured

								> female  <- c(201, 71)
								>  totals <- c( (114+201), (71+87))
								>  prop.test(female, totals)

								2-sample test for equality of proportions with continuity correction

								data:  female out of totals
								X-squared = 14.574, df = 1, p-value = 0.0001347
								alternative hypothesis: two.sided
								95 percent confidence interval:
								0.08999683 0.28745947
								sample estimates:
								   prop 1    prop 2 
								0.6380952 0.4493671 

#postdoc vs untenured
								>  female  <- c(87,64)
								>  totals <- c((87+75), (64+43))
								>  prop.test(female, totals)

								2-sample test for equality of proportions with continuity correction
								data:  female out of totals
								X-squared = 0.7444, df = 1, p-value = 0.3883
								alternative hypothesis: two.sided
								95 percent confidence interval:
								 -0.18937408  0.06718647
								sample estimates:
								   prop 1    prop 2 
								0.5370370 0.5981308 

#postdoc vs tenured
								> female  <- c(87,71)
								>  totals <- c((87+75), (71+87))
								>  prop.test(female, totals)

								2-sample test for equality of proportions with continuity correction

								data:  female out of totals
								X-squared = 2.1213, df = 1, p-value = 0.1453
								alternative hypothesis: two.sided
								95 percent confidence interval:
								-0.02772119  0.20306109
								sample estimates:
								   prop 1    prop 2 
								0.5370370 0.4493671 

#untenured vs tenured
								>  prop.test(female, totals)

								2-sample test for equality of proportions with continuity correction

								data:  female out of totals
								X-squared = 5.0698, df = 1, p-value = 0.02435
								alternative hypothesis: two.sided
								95 percent confidence interval:
								0.01990738 0.27762012
								sample estimates:
								   prop 1    prop 2 
								0.5981308 0.4493671 


##########################
##########################
# Question 2: Orientation
##########################
##########################

#####################
# TABLE S4.   #
###############
> table(data_tri_soc $Q1, data_tri_soc $Q3)
                                                                                                                       
                                                                                                                            Choose not to answer Gay, lesbian, bisexual, pansexual or asexual
                                                                                                                          1                    1                                            0
  American Society of Naturalists (ASN)                                                                                   0                    3                                            9
  American Society of Naturalists (ASN),Society for the Study of Evolution (SSE)                                          0                    3                                           25
  American Society of Naturalists (ASN),Society of Systematic Biologists (SSB)                                            0                    0                                            1
  American Society of Naturalists (ASN),Society of Systematic Biologists (SSB),Society for the Study of Evolution (SSE)   0                    1                                            5
  None                                                                                                                    0                    1                                            8
  Society for the Study of Evolution (SSE)                                                                                1                   17                                           68
  Society of Systematic Biologists (SSB)                                                                                  0                    0                                            9
  Society of Systematic Biologists (SSB),Society for the Study of Evolution (SSE)                                         0                    1                                           12
                                                                                                                       
                                                                                                                        Straight or heterosexual
                                                                                                                                               1
  American Society of Naturalists (ASN)                                                                                                       45
  American Society of Naturalists (ASN),Society for the Study of Evolution (SSE)                                                             131
  American Society of Naturalists (ASN),Society of Systematic Biologists (SSB)                                                                 4
  American Society of Naturalists (ASN),Society of Systematic Biologists (SSB),Society for the Study of Evolution (SSE)                       37
  None                                                                                                                                        39
  Society for the Study of Evolution (SSE)                                                                                                   310
  Society of Systematic Biologists (SSB)                                                                                                      50
  Society of Systematic Biologists (SSB),Society for the Study of Evolution (SSE)                                                             69




############################
#####LGBT vs census#########
#############################
LGBT  <- c(137, 328239523*0.049)
totals <- c( (823), (328239523))
					 prop.test(LGBT, totals)
					
						2-sample test for equality of proportions with continuity correction
					
					data:  lgbt out of totals
					X-squared = 241.17, df = 1, p-value < 2.2e-16
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 0.09140761 0.14352070
					sample estimates:
					   prop 1    prop 2 
					0.1664642 0.0490000 


> lgbt  <- c(16, 15)
> totals <- c( (158), (107))
> prop.test(lgbt, totals)

					2-sample test for equality of proportions with continuity correction
					
					data:  lgbt out of totals
					X-squared = 0.59672, df = 1, p-value = 0.4398
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.12762958  0.04978739
					sample estimates:
					   prop 1    prop 2 
					0.1012658 0.1401869 

> lgbt  <- c(16, 15)
> totals <- c( (158), (156))
> prop.test(lgbt, totals)

					2-sample test for equality of proportions with continuity correction
					data:  lgbt out of totals
					X-squared = 2.2626e-31, df = 1, p-value = 1
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.06597612  0.07620007
					sample estimates:
					    prop 1     prop 2 
					0.10126582 0.09615385 
					
> 79+238
[1] 317
> lgbt  <- c(79, 15)
> totals <- c( (238), (156))
> prop.test(lgbt, totals)

					2-sample test for equality of proportions with continuity correction
					data:  lgbt out of totals
					X-squared = 27.554, df = 1, p-value = 1.527e-07
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 0.1548468 0.3167111
					sample estimates:
					    prop 1     prop 2 
					0.33193277 0.09615385 

> 
> lgbt  <- c(79, 16)
> totals <- c( (238), (158))
> prop.test(lgbt, totals)

					2-sample test for equality of proportions with continuity correction
					data:  lgbt out of totals
					X-squared = 26.458, df = 1, p-value = 2.694e-07
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 0.1492965 0.3120374
					sample estimates:
					   prop 1    prop 2 
					0.3319328 0.1012658 

> lgbt  <- c(79, 15)
> totals <- c( (238), (107))
> prop.test(lgbt, totals)

					2-sample test for equality of proportions with continuity correction
					data:  lgbt out of totals
					X-squared = 12.741, df = 1, p-value = 0.0003578
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 0.09605306 0.28743865
					sample estimates:
					   prop 1    prop 2 
					0.3319328 0.1401869 

> 
> lgbt  <- c(79, 11)
> totals <- c( (238), (156))
> prop.test(lgbt, totals)

					2-sample test for equality of proportions with continuity correction
					data:  lgbt out of totals
					X-squared = 35.071, df = 1, p-value = 3.179e-09
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 0.1840504 0.3387895
					sample estimates:
					    prop 1     prop 2 
					0.33193277 0.07051282 

##############
##############		
#FIGURE 3#####
##############
##############
binom.bayes(201, 326, conf.level = 0.95)
binom.bayes(87, 169, conf.level = 0.95)
binom.bayes(64, 107, conf.level = 0.95)
binom.bayes(71, 158, conf.level = 0.95)



pdf(file="Figure3_SSEAv02.pdf", height=5, width=10)

		par(mfrow=c(1,2))

						means_women<-c(0.616208, 0.5147059, 0.5972222, 0.4496855)		
						upper_women<-c(0.6685923, 0.5894719, 0.6885701, 0.5268334)		
						lower_women<-c(0.5634259, 0.4398434, 0.504873, 0.3728905)

		plotCI(means_women, li= lower_women, ui= upper_women, ylim=c(0,1), col="gray19", lwd=1, pch=1,lty=2, cex=2, xlab="Professional stage", ylab="Proportion")
		legend("topleft",legend=c("Women", "Men"), col=c("gray43", "black"), lty=1:1, cex=1)

						means_LBGT<-c(0.2431193, 0.09705882, 0.1409091, 0.07142857)		
						upper_LBGT<-c(0.2898894, 0.1421273, 0.206581, 0.1117308)		
						lower_LBGT<-c(0.19723, 0.05470582, 0.07895607, 0.03422507)

		plotCI(means_LBGT, li= lower_LBGT, ui= upper_LBGT, ylim=c(0,1), col="tan2", lwd=1, pch=15, cex=2)
		legend("topleft",legend=c("LGBTQ+", "Heterosexual"), col=c("tan2", "turquoise4"), lty=1:1, cex=1)

dev.off()

#################
# TABLE S5      #
#################

###########################
# LGBT professional stages
###########################
> table(data_tri_soc $Q6, data_tri_soc $Q3)
                                                     
                                                          Choose not to answer Gay, lesbian, bisexual, pansexual or asexual Straight or heterosexual
                                                        1                    1                                            0                        0
  Academic professional (staff / administration)        0                    1                                            1                       17
  Graduate student (Masters or PhD)                     1                    8                                           79                      238
  Non-academic professional (publisher, vendor, etc.)   0                    0                                            1                        9
  Non-tenure track faculty                              0                    0                                            3                       13
  Other                                                 0                    0                                            1                        6
  Postdoc                                               0                   11                                           16                      142
  Pre-tenure faculty                                    0                    2                                           15                       92
  Tenured faculty                                       0                    4                                           11                      145
  Undergraduate                                         0                    0                                           10                       24


#############################################
#
#LGBT: BETWEEN SOCIETIES COMPARISONS 
#
#############################################
FEMALES
TOTAL ASN: 33+90+3+22 = 148
TOTAL SSB: 3+22+30+41 = 96
TOTAL SSE: 90+22+222+41 = 375

MEMBERS:
TOTAL ASN: 55+158+5+41 = 259
TOTAL SSB: 5+41+59+78 = 183
TOTAL SSE: 158+41+382+78 = 659

SSE vs. ASN
###WITH OVERLAP###
female  <- c(375, 148)
totals <- c( 659, 259)
prop.test(female, totals)

				X-squared = 2.048e-29, df = 1, p-value = 1

SSE vs. SSB
###WITH OVERLAP###
female  <- c(375, 96)
totals <- c( 659, 183)
prop.test(female, totals)
			
				X-squared = 0.97507, df = 1, p-value = 0.3234


SSB vs. ASN
###WITH OVERLAP###
female  <- c(96, 148)
totals <- c(183, 259)
prop.test(female, totals)

				X-squared = 0.77133, df = 1, p-value = 0.3798

###SSE ASN; WITH NO OVERLAP###
female  <- c(222, 33)
totals <- c( 382, 55)
prop.test(female, totals)

				X-squared = 0.01412, df = 1, p-value = 0.9054


###SSE SSB; WITH NO OVERLAP###
female  <- c(222, 30)
totals <- c( 382, 59)
prop.test(female, totals)

				X-squared = 0.82548, df = 1, p-value = 0.3636


###ASN SSB; WITH NO OVERLAP###
female  <- c(33, 30)
totals <- c( 55, 59)
prop.test(female, totals)

				X-squared = 0.6298, df = 1, p-value = 0.4274

########################
#####LGBQT##############
#########################

###SSE ASN; WITH OVERLAP###
					lgbt  <- c(110, 40)
					totals <- c( (110+547), (40+217))
					prop.test(lgbt, totals)

					X-squared = 0.11101, df = 1, p-value = 0.739


###SSE SSB; WITH OVERLAP###
					lgbt  <- c(110, 27)
					totals <- c( (110+547), (27+160))
					prop.test(lgbt, totals)

					X-squared = 0.41159, df = 1, p-value = 0.5212


###ASN SSB; WITH OVERLAP###
					lgbt  <- c(40, 27)
					totals <- c( (40+217),(27+160))
					prop.test(lgbt, totals)

					X-squared = 0.03722, df = 1, p-value = 0.847

###SSE ASN; WITH NO OVERLAP###
					lgbt  <- c(68, 9)
					totals <- c( (68+310), (9+45))
					prop.test(lgbt, totals)

					X-squared = 0.0022577, df = 1, p-value = 0.9621


###SSE SSB; WITH NO OVERLAP###
					lgbt  <- c(68, 9)
					totals <- c( (68+310), (9+50))
					prop.test(lgbt, totals)

					X-squared = 0.10834, df = 1, p-value = 0.742

###ASN SSB; WITH NO OVERLAP###
					lgbt  <- c(9, 9)
					totals <- c( (9+45), (9+50))
					prop.test(lgbt, totals)
					
					X-squared = 1.1879e-31, df = 1, p-value = 1

############################
############################
#QUestion 3: Underrepresented
############################
############################
#TABLE S6: SOC comparisons ALLOWING FOR OVERLAP
                                                        None ASN    SSE    SSB SSE/SSB	   ASN/SSE ASN/SSB    ASM/SSB/SSE
Black or African American                                  *   0     4       3   *	 	    0		  0				  *
Hispanic or Latinx                                        10   9    63       5 	18		    5		  *				  3
Multi-racial                                               3   2    13       3 	 0 		    8		  *				  0
Native Hawaiian or other Pacific islander                  0   0     0       0 	 0 		    0		  0				  *
Other                                                      3   *     7       1 	 2 		    3		  *				  0
South Asian, East Asian, or Southeast Asian                3   6    47       6 	 5 		    8		  0				  *
White, non-Hispanic                                       26   37  249      41 	54		  133		  2				  37
American Indian, Alaskan Native, Indigenous			       *   0     *       0	 1		    0		  0				  0
White, non-Hispanic,Other                                  0   0	 *       0 	 0	 		0		  0				  0
South Asian, East Asian, or Southeast Asian                0   0	 0       0 	 0			*		  0				  0
South Asian, East Asian, or Southeast Asian,Multi-racial   0   0	 *       0 	 0 			0		  0				  0
Total													  47  55   386      59  81        158         5              43                

white ASN:37+133+2+37 = 209
white SSE:249+54+133+37 = 473
white SSB: 41+54+2+37 = 134

asian ASN: 6+8+1+1+1= 17
asian SSE: 47+5+1+1+1= 55
asian SSB: 6+8+1= 15

BLACK ASN 1
BLACK SSE 4+1+1 = 6
BLACK SSB 3+1+1 = 5


Hispanic ASN: 9+5+1+3 = 18
Hispanic SSB: 5+18+1+3 = 27
Hispanic SSE: 63+18+5+3 = 89

Multiracial ASN: 2+8+1 = 11
Multiracial SSB: 3+1 =4
Multiracial SSE: 13+8 =21

Total ASN: 55+158+5+43 = 261
Total SSB: 59+81+5+43 = 188
Total SSE: 386+81+ 158 +43 = 668

###########
#TABLE S7##
###########

##NO OVERLAP
##HISPANICS

#ASN vs SSE
					hispanics  <- c( (8), (63))
					totals <- c( 55, 396)
					prop.test(hispanics, totals)

					X-squared = 0.028571, df = 1, p-value = 0.8658

#ASN vs SSB
					hispanics  <- c( (8), (10))
					totals <- c( 57, 59)
					prop.test(hispanics, totals)

				X-squared = 0.031286, df = 1, p-value = 0.8596

#SSB vs SSE
					hispanics  <- c( (10), (63))
					totals <- c( 59, 396)
					prop.test(hispanics, totals)

					X-squared = 0.00016778, df = 1, p-value = 0.9897

##BLACKS

#ASN vs SSE
					blacks  <- c( (0), (4))
					totals <- c( 55, 385)
					prop.test(blacks, totals)

					X-squared = 0, df = 1, p-value = 1

#ASN vs SSB
					blacks  <- c( (0), (3))
					totals <- c( 55, 59)
					prop.test(blacks, totals)

					X-squared = 1.2305, df = 1, p-value = 0.2673


#SSB vs SSE
					blacks  <- c( (3), (4))
					totals <- c( 59, 385)
					prop.test(blacks, totals)

					X-squared = 3.1042, df = 1, p-value = 0.07809

##MULTIRACIAL

#ASN vs SSE
					multi  <- c( (2), (14))
					totals <- c( 55, 385)
					prop.test(multi, totals)

					X-squared = 0, df = 1, p-value = 1

#ASN vs SSB
					multi  <- c( (2), (3))
					totals <- c( 55, 59)
					prop.test(multi, totals)

					X-squared = 5.6891e-31, df = 1, p-value = 1


#SSB vs SSE
					multi  <- c( (3), (13))
					totals <- c( 59, 385)
					prop.test(multi, totals)
					X-squared = 0.078654, df = 1, p-value = 0.7791



#####HISPANICS, OVERLAP###
#ASN vs SSB
					hisp  <- c(18, 27)
					totals <- c( (261), (188))
					prop.test(hisp, totals)
					X-squared = 5.951, df = 1, p-value = 0.01471

#ASN vs SSE
					hisp  <- c(18, 89)
					totals <- c( (261), (668))
					prop.test(hisp, totals)
					X-squared = 6.9886, df = 1, p-value = 0.008203

#SSB vs SSE
					vhisp  <- c(27, 89)
					totals <- c( (188), (668))
					prop.test(hisp, totals)
					X-squared = 0.060934, df = 1, p-value = 0.805


#####BLACKS, OVERLAP###
#ASN vs SSB
					blacks  <- c(1, 5)
					totals <- c( (261), (188))
					prop.test(blacks, totals)
					X-squared = 2.7423, df = 1, p-value = 0.09773

#ASN vs SSE
					blacks  <- c(1, 6)
					totals <- c( (261), (668))
					prop.test(blacks, totals)
					X-squared = 0.15646, df = 1, p-value = 0.6924

#SSB vs SSE
					blacks  <- c(5, 6)
					totals <- c( (188), (668))
					prop.test(blacks, totals)
					X-squared = 2.3339, df = 1, p-value = 0.1266

#####Multiracial, OVERLAP###
#ASN vs SSB
					multi  <- c(11, 4)
					totals <- c( (261), (188))
					prop.test(multi, totals)
					X-squared = 0.89847, df = 1, p-value = 0.3432

#ASN vs SSE
					multi  <- c(11, 21)
					totals <- c( (261), (668))
					prop.test(multi, totals)
					X-squared = 0.36514, df = 1, p-value = 0.5457


#SSB vs SSE
					multi  <- c(4, 21)
					totals <- c( (188), (668))
					prop.test(multi, totals)
					X-squared = 0.23593, df = 1, p-value = 0.6272

####LATINO in trisociety VS NSF

NSF from https://ncses.nsf.gov/pubs/nsf21308/data-tables table 22
total hispanics: 9+63+5+18+5+1+3 = 104
total members: 47 + 55 +  386 +     59 + 81    +    158  +       5   +           43 = 834

#vs. Life sciences recipeints

latino  <- c(600, 104)
totals <- c(6380, 834)
prop.test(latino, totals)


					2-sample test for equality of proportions with continuity correction
					data:  latino out of totals
					X-squared = 7.5271, df = 1, p-value = 0.006078
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.054872552 -0.006440153
					sample estimates:
					    prop 1     prop 2 
					0.09404389 0.12470024 

#vs. ev bio recipeints

 latino  <- c(15, 104)
 totals <- c( 192, 834)
 prop.test(latino, totals)


				2-sample test for equality of proportions with continuity correction
				data:  latino out of totals
				X-squared = 2.8633, df = 1, p-value = 0.09062
				alternative hypothesis: two.sided
				95 percent confidence interval:
				 -0.0938666734  0.0007161938
				sample estimates:
				   prop 1    prop 2 
				0.0781250 0.1247002 


library(pwr)

		h=2*asin(sqrt(104/834))-2*asin(sqrt(600/6380))= 0.09845135
		h=2*asin(sqrt(104/834))-2*asin(sqrt(15/192))= 0.1552631	
		pwr.2p2n.test(h = 0.09845135, n1 =834 , n2 =6380 , sig.level = , power = ) 
         				 power = 0.7623338
		pwr.2p2n.test(h = 0.1552631, n1 =834 , n2 =192 , sig.level = , power = )
         				 power = 0.4919532

####BLACK in trisociety VS NSF

NSF from https://ncses.nsf.gov/pubs/nsf21308/data-tables table 22
total blacks: 4    +   3  + 1	 +	    0		+  0	+			  1 = 9
total members: 47 + 55 +  386 +     59 + 81    +    158  +       5   +           43 = 834

#vs. Life sciences recipeints

 black  <- c(276, 9)
 totals <- c(6380, 834)
 prop.test(black, totals)
					2-sample test for equality of proportions with continuity correction
					data:  black out of totals
					X-squared = 19.645, df = 1, p-value = 9.324e-06
					alternative hypothesis: two.sided
					95 percent confidence interval:
					0.02318338 0.04175426
					sample estimates:
					    prop 1     prop 2 
					0.04326019 0.01079137 

#vs. ev bio recipients

 black  <- c(3, 9)
 totals <- c(192, 834)
 prop.test(black, totals)

					2-sample test for equality of proportions with continuity correction
					data:  black out of totals
					X-squared = 0.035871, df = 1, p-value = 0.8498
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.01726194  0.02692921
					sample estimates:
					    prop 1     prop 2 
					0.01562500 0.01079137 



#################################
#################################
#Societies comparisons : URM
#################################
#################################
BLACK ASN 1
BLACK SSE 4+1+1 = 6
BLACK SSB 3+1+1 = 5


Hispanic ASN: 9+5+1+3 = 18
Hispanic SSB: 5+18+1+3 = 27
Hispanic SSE: 63+18+5+3 = 89

Multiracial ASN: 2+8+1 = 11
Multiracial SSB: 3+1 =4
Multiracial SSE: 13+8 =21

Total ASN: 55+158+5+43 = 261
Total SSB: 59+81+5+43 = 188
Total SSE: 386+81+ 158 +43 = 668

#####HISPANICS, OVERLAP###
			#ASN vs SSB
			hisp  <- c(18, 27)
			totals <- c( (261), (188))
			prop.test(hisp, totals)
			X-squared = 5.951, df = 1, p-value = 0.01471
			
#ASN vs SSE
			hisp  <- c(18, 89)
			totals <- c( (261), (668))
			prop.test(hisp, totals)
			
			X-squared = 6.9886, df = 1, p-value = 0.008203

#SSB vs SSE
			hisp  <- c(27, 89)
			totals <- c( (188), (668))
			prop.test(hisp, totals)
			
			X-squared = 0.060934, df = 1, p-value = 0.805


#####BLACKS, OVERLAP###
#ASN vs SSB
			blacks  <- c(1, 5)
			totals <- c( (261), (188))
			prop.test(blacks, totals)
			X-squared = 2.7423, df = 1, p-value = 0.09773

#ASN vs SSE
			blacks  <- c(1, 6)
			totals <- c( (261), (668))
			prop.test(blacks, totals)
			X-squared = 0.15646, df = 1, p-value = 0.6924

#SSB vs SSE
			blacks  <- c(5, 6)
			totals <- c( (188), (668))
			prop.test(blacks, totals)
			
			X-squared = 2.3339, df = 1, p-value = 0.1266

#####Multiracial, OVERLAP###
#ASN vs SSB
			multi  <- c(11, 4)
			totals <- c( (261), (188))
			prop.test(multi, totals)
			X-squared = 0.89847, df = 1, p-value = 0.3432
			
#ASN vs SSE
			multi  <- c(11, 21)
			totals <- c( (261), (668))
			prop.test(multi, totals)
			
			X-squared = 0.36514, df = 1, p-value = 0.5457


#SSB vs SSE
			multi  <- c(4, 21)
			totals <- c( (188), (668))
			prop.test(multi, totals)
			X-squared = 0.23593, df = 1, p-value = 0.6272


########################################
#URM Vs census and NSF PH.D. information
########################################
####Hispanic in trisociety VS NSF

NSF from https://ncses.nsf.gov/pubs/nsf21308/data-tables table 22
total hispanics: 9+63+5+18+5+1+3 = 104
total members: 47 + 55 +  386 +     59 + 81    +    158  +       5   +           43 = 834

		#vs. Life sciences recipeints

				 hispanics  <- c(600, 104)
				 totals <- c(6380, 834)
				 prop.test(hispanics, totals)
				
				
					2-sample test for equality of proportions with continuity correction
				
				data:  hispanics out of totals
				X-squared = 7.5271, df = 1, p-value = 0.006078
				alternative hypothesis: two.sided
				95 percent confidence interval:
				 -0.054872552 -0.006440153
				sample estimates:
				    prop 1     prop 2 
				0.09404389 0.12470024 

		#vs. ev bio recipients

				  latino  <- c(15, 104)
				  totals <- c( 192, 834)
				  prop.test(latino, totals)
				 
				 	2-sample test for equality of proportions with continuity correction
				 
				 data:  latino out of totals
				 X-squared = 2.8633, df = 1, p-value = 0.09062
				 alternative hypothesis: two.sided
				 95 percent confidence interval:
				  -0.0938666734  0.0007161938
				 sample estimates:
				    prop 1    prop 2 
				 0.0781250 0.1247002 

####BLACK in trisociety VS NSF

NSF from https://ncses.nsf.gov/pubs/nsf21308/data-tables table 22
total blacks: 4    +   3  + 1	 +	    0		+  0	+			  1 = 9
total members: 47 + 55 +  386 +     59 + 81    +    158  +       5   +           43 = 834

		#vs. Life sciences recipeints

				  black  <- c(276, 9)
				  totals <- c(6380, 834)
				  prop.test(black, totals)
				 	2-sample test for equality of proportions with continuity correction
				 
				 data:  black out of totals
				 X-squared = 19.645, df = 1, p-value = 9.324e-06
				 alternative hypothesis: two.sided
				 95 percent confidence interval:
				  0.02318338 0.04175426
				 sample estimates:
				     prop 1     prop 2 
				 0.04326019 0.01079137 

		#vs. ev bio recipeints

				  black  <- c(3, 9)
				  totals <- c(192, 834)
				  prop.test(black, totals)
				 
				 	2-sample test for equality of proportions with continuity correction
				 
				 data:  black out of totals
				 X-squared = 0.035871, df = 1, p-value = 0.8498
				 alternative hypothesis: two.sided
				 95 percent confidence interval:
				  -0.01726194  0.02692921
				 sample estimates:
				     prop 1     prop 2 
				 0.01562500 0.01079137 

###########################################
#PROFESSIONAL STAGES COMPARISON HISPANICS
###########################################
#grad vs postdoc
					hispa  <- c(33, 20)
					 totals <- c(234, 128)
					 prop.test(hispa, totals)
						2-sample test for equality of proportions with continuity correction

					data:  hispa out of totals
					X-squared = 0.05581, df = 1, p-value = 0.8132
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.09837269  0.06792397
					sample estimates:
					   prop 1    prop 2 
					0.1410256 0.1562500 

#grad vs untenured
					hispa  <- c(33, 9)
					 totals <- c(234, 77)
					 prop.test(hispa, totals)

					2-sample test for equality of proportions with continuity correction
					data:  hispa out of totals
					X-squared = 0.11935, df = 1, p-value = 0.7297
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.06897604  0.11726109
					sample estimates:
					   prop 1    prop 2 
					0.1410256 0.1168831 

#grad vs tenured
					hispa  <- c(33, 14)
					 totals <- c(234, 132)
					 prop.test(hispa, totals)

					2-sample test for equality of proportions with continuity correction
					data:  hispa out of totals
					X-squared = 0.6359, df = 1, p-value = 0.4252
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.03986431  0.10979438
					sample estimates:
					   prop 1    prop 2 
					0.1410256 0.1060606 

# postdoc vs untenured
					hispa  <- c(20, 9)
					 totals <- c(128, 77)
					 prop.test(hispa, totals)

						2-sample test for equality of proportions with continuity correction
					data:  hispa out of totals
					X-squared = 0.33217, df = 1, p-value = 0.5644
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.0664593  0.1451931
					sample estimates:
					   prop 1    prop 2 
					0.1562500 0.1168831 

# postdoc vs tenured
					hispa  <- c(20, 14)
					 totals <- c(128, 132)
					 prop.test(hispa, totals)

						2-sample test for equality of proportions with continuity correction
					data:  hispa out of totals
					X-squared = 1.0324, df = 1, p-value = 0.3096
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.03945465  0.13983344
					sample estimates:
					   prop 1    prop 2 
					0.1562500 0.1060606 

# untenured vs tenured
					hispa  <- c(9, 14)
					 totals <- c(77, 132)
					 prop.test(hispa, totals)

					2-sample test for equality of proportions with continuity correction
					data:  hispa out of totals
					X-squared = 0.0001454, df = 1, p-value = 0.9904
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.08839048  0.11003550
					sample estimates:
					   prop 1    prop 2 
					0.1168831 0.1060606 

####################################
####################################
#QUESTION 4: DISSABILITIES
####################################
####################################


> (31+27+22)/(892+637+546)
[1] 0.03855422
> 3930/(44209+2930)
[1] 0.08337046
> 3930/(44209+3930)
[1] 0.08163859
> 567/(567+7354)
[1] 0.07158187
> disa  <- c(91, 567)
> totals <- c( (91+752), (567+7354))
> prop.test(disa, totals)

					2-sample test for equality of proportions with continuity correction
					data:  female out of totals
					X-squared = 13.991, df = 1, p-value = 0.0001837
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 0.01400632 0.05872555
					sample estimates:
					    prop 1     prop 2 
					0.10794781 0.07158187 

 
> 40747411/(40747411+277428456)
[1] 0.1280657
> disa  <- c(91, 40747411)
> totals <- c( (91+752), (40747411+277428456))
> prop.test(disa, totals)

					2-sample test for equality of proportions with continuity correction
					data:  disa out of totals
					X-squared = 2.8779, df = 1, p-value = 0.0898
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.041658735  0.001422972
					sample estimates:
					   prop 1    prop 2 
					0.1079478 0.1280657 

> 
> female  <- c(567, 40747411)
> totals <- c( (567+7354), (40747411+277428456))
> prop.test(female, totals)

					2-sample test for equality of proportions with continuity correction
					data:  female out of totals
					X-squared = 225.8, df = 1, p-value < 2.2e-16
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.06222422 -0.05074341
					sample estimates:
					    prop 1     prop 2 
					0.07158187 0.12806569 

#DISSABILITY BREAKDOWN.
#CEnsus data from file:///Users/danielmatute/Downloads/wmpd19-sr-tab01-003.pdf
#vs census
with disabibility: 40,747,411
total: 318,175,867

					disa_total  <- c(91, 40747411)
					totals <- c( (843), (318175867))
					prop.test(disa_total, totals)
					data:  disa_total out of totals
					data:  disa_total out of totals
					X-squared = 2.8779, df = 1, p-value = 0.0898
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.041658735  0.001422972
					sample estimates:
					   prop 1    prop 2 
					0.1079478 0.1280657 

#census vs biological sciences phd

 
> female  <- c(567, 40747411)
> totals <- c( (567+7354), (40747411+277428456))
> prop.test(female, totals)

					2-sample test for equality of proportions with continuity correction
					data:  female out of totals
					X-squared = 225.8, df = 1, p-value < 2.2e-16
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.06222422 -0.05074341
					sample estimates:
					    prop 1     prop 2 
					0.07158187 0.12806569 


#survey vs biological sciences

> 567/(567+7354)
[1] 0.07158187
> disa  <- c(91, 567)
> totals <- c( (91+752), (567+7354))
> prop.test(disa, totals)

					2-sample test for equality of proportions with continuity correction
					data:  female out of totals
					X-squared = 13.991, df = 1, p-value = 0.0001837
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 0.01400632 0.05872555
					sample estimates:
					    prop 1     prop 2 
					0.10794781 0.07158187 


#TABLE S8#
                                             Do not wish to provide                Hearing impairment Mobility or orthopedic Impairment                              None 
                                9                                25                                15                                 5                               752 
           Other (please specify)                     Physiological                 Visual impairment 
                               23                                10                                13 

25+15+5+752+23+10+13= 843

> hearing  <- c(15, 79)
> totals <- c( (843), (7921))
> prop.test(hearing, totals)

					2-sample test for equality of proportions with continuity correction
					data:  female out of totals
					X-squared = 3.6852, df = 1, p-value = 0.0549
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.00202468  0.01766489
					sample estimates:
					     prop 1      prop 2 
					0.017793594 0.009973488 

> hearing  <- c(15, 11515283)
> totals <- c( (843), (320775014))
> prop.test(hearing, totals)

					2-sample test for equality of proportions with continuity correction
					data:  hearing out of totals
					X-squared = 7.4693, df = 1, p-value = 0.006276
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.027622031 -0.008587414
					sample estimates:
					    prop 1     prop 2 
					0.01779359 0.03589832 

#Visual

visual  <- c(13, 214)
totals <- c( (843), (7921))
> prop.test(visual, totals)

					2-sample test for equality of proportions with continuity correction
					data:  female out of totals
					X-squared = 3.6138, df = 1, p-value = 0.0573
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.021303827 -0.001887524
					sample estimates:
					    prop 1     prop 2 
					0.01542112 0.02701679 

> visual  <- c(13, 7555551)
> totals <- c( (843), (320775014))
> prop.test(visual, totals)

					2-sample test for equality of proportions with continuity correction
					data:  female out of totals
					X-squared = 2.0837, df = 1, p-value = 0.1489
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.0170440457  0.0007781741
					sample estimates:
					    prop 1     prop 2 
					0.01542112 0.02355405 

#MOBILITY

#Vs Life sciences
> mobility  <- c(5, 32)
> totals <- c( (843), (7921))
> prop.test(mobility, totals)

					2-sample test for equality of proportions with continuity correction
					data:  female out of totals
					X-squared = 0.27645, df = 1, p-value = 0.599
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.004133256  0.007915865
					sample estimates:
					     prop 1      prop 2 
					0.005931198 0.004039894 

#Vs census

> mobility  <- c(5, 20903105)
> totals <- c( (843), (320775014))

> prop.test(mobility, totals)

					2-sample test for equality of proportions with continuity correction
					data:  female out of totals
					X-squared = 47.585, df = 1, p-value = 5.267e-12
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.06500976 -0.05345660
					sample estimates:
					     prop 1      prop 2 
					0.005931198 0.065164380 
					
#DIsability Representation across professional stages
#grad vs postdoc
					all_disa  <- c(40, 21)
					totals <- c( (40+284), (21+146))
					prop.test(all_disa, totals)

					data:  all_disa out of totals
					X-squared = 2.1036e-30, df = 1, p-value = 1
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.06632372  0.06174029
					sample estimates:
					   prop 1    prop 2 
					0.1234568 0.1257485 

#grad vs untenured

					> all_disa  <- c(40, 8)
					> totals <- c( (40+284), (8+98))
					> prop.test(all_disa, totals)
					2-sample test for equality of proportions with continuity correction

					data:  all_disa out of totals
					X-squared = 1.4022, df = 1, p-value = 0.2364
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.02001425  0.11598443
					sample estimates:
					   prop 1    prop 2 
					0.1234568 0.0754717 

#grad vs tenure
					> all_disa  <- c(40, 14)
					> totals <- c( (40+284), (14+145))
					> prop.test(all_disa, totals)

					2-sample test for equality of proportions with continuity correction
					data:  all_disa out of totals
					X-squared = 1.0135, df = 1, p-value = 0.3141
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.02605316  0.09686611
					sample estimates:
					    prop 1     prop 2 
					0.12345679 0.08805031 

#postdoc vs untenured
					all_disa  <- c(21, 8)
					totals <- c( (21+146), (8+98))
					prop.test(all_disa, totals)

					2-sample test for equality of proportions with continuity correction
					data:  all_disa out of totals
					X-squared = 1.2374, df = 1, p-value = 0.266
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.02855036  0.12910397
					sample estimates:
					   prop 1    prop 2 
					0.1257485 0.0754717 

#postdoc vs tenured
					all_disa  <- c(21, 14)
					totals <- c( (21+146), (14+145))
					prop.test(all_disa, totals)

					2-sample test for equality of proportions with continuity correction
					data:  all_disa out of totals
					X-squared = 0.84651, df = 1, p-value = 0.3575
					alternative hypothesis: two.sided
					95 percent confidence interval:
					 -0.03528973  0.11068610
					sample estimates:
					    prop 1     prop 2 
					0.12574850 0.08805031 

#untenured vs tenured
					all_disa  <- c(8, 14)
					totals <- c( (8+98), (14+145))
					prop.test(all_disa, totals)
					2-sample test for equality of proportions with continuity correction

					data:  all_disa out of totals
					X-squared = 0.018589, df = 1, p-value = 0.8916
					alternative hypothesis: two.sided
					95 percent confidence interval:
					-0.08728844  0.06213121
					sample estimates:
					    prop 1     prop 2 
					0.07547170 0.08805031 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
