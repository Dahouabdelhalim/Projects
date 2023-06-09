
 library("smatr")

 # The following tests are conducted on the dataset "Ontogenetic.Allometry.Data". 
 # For each figure and sub-figure, the script generates the relevent dataset for slope comparisons.
 # The slope comparison test is provided at the end of the script.
 # The uncorrected P-value is provided after each dataset. 
 # For significant P-values for data used in multiple tests, a Bonferonni-corrected P-value is also provided.
 # The Bonferonni correction is the uncorrected P-value multiplied by the highest number of hypotheses tested with one genotype.
 # The number of comparisons for each genotype is: 
 # ci-WT: 4
 # ci-Minute(3R): 3
 # ci-Mosaic(3R): 5
 # ci-Mosaic(3R)EtOH: 3
 # ci-Inr.DN-Mosaic(3R)EtOH: 2
 # All other genotypes used in only one comparison.
 
 #Fig 1B
 testdata<-subset(Ontogenetic.Allometry.Data, Genotype=="ci-WT"|Genotype=="ci-Minute(3R)")
      #uncorrected P=0.8486027 
 
 #Fig 1B'
 testdata<-subset(Ontogenetic.Allometry.Data, Genotype=="ci-WT"|Genotype=="ci-Mosaic(3R)")
      #uncorrected P=0.0004406042 corrected P=0.002203021
 testdata<-subset(Ontogenetic.Allometry.Data, Genotype=="ci-Minute(3R)"|Genotype=="ci-Mosaic(3R)")
      #uncorrected P=0.0007549755 corrected P=0.003774878
 
 #Fig 2A
 testdata<-subset(Ontogenetic.Allometry.Data, Genotype=="ci-Mosaic(3R)EtOH"|Genotype=="ci-Mosaic(3R)20E")
      #uncorrected P=0.004444931 corrected P=0.01333479
 
 #Fig 2B
 testdata<-subset(Ontogenetic.Allometry.Data, Genotype=="ci-Minute(3R)EtOH"|Genotype=="ci-Minute(3R)20E")
      #uncorrected P=0.9850054
 
 #Fig 3A
 testdata<-subset(Ontogenetic.Allometry.Data, Genotype=="ci-EcR.RNAi(II)"|Genotype=="ci-WT17oC")
      #uncorrected P=0.02425538
 
 #Fig 4A
 testdata<-subset(Ontogenetic.Allometry.Data, Genotype=="ci-WT"|Genotype=="ci-Inr.del"|Genotype=="ci-Inr.DN")
      #uncorrected P=0.4898537
 
 #Fig 4B
 testdata<-subset(Ontogenetic.Allometry.Data, Genotype=="ci-Mosaic(3R)"|Genotype=="ci-Inr.del-Mosaic(3R)"|Genotype=="ci-Inr.DN-Mosaic(3R)")
      #uncorrected P=2.009307e-07 corrected P= 1.004653e-06
 
 #Fig 4C
 testdata<-subset(Ontogenetic.Allometry.Data, Genotype=="ci-Inr.DN-Mosaic(3R)EtOH"|Genotype=="ci-Inr.DN-Mosaic(3R)20E")
      #uncorrected P=0.3942004
 
 #Sup Fig 2
 testdata<-subset(Ontogenetic.Allometry.Data, Genotype=="en-Mosaic(3L)"|Genotype=="Minute(3L)"|Genotype=="en-RFP")
      #uncorrected P=0.003418632
 
 #Sup Fig 3
 testdata<-subset(Ontogenetic.Allometry.Data, Genotype=="ci-Mosaic(3R)"|Genotype=="ci-Mosaic(3R)-Replicate 2")
      #uncorrected P=0.9619436
 testdata<-subset(Ontogenetic.Allometry.Data, Genotype=="ci-WT-Replicate 2"|Genotype=="ci-WT")
      #uncorrected P=0.9508046
 
 #Compare effect of EtOH 
 testdata<-subset(Ontogenetic.Allometry.Data, Genotype=="ci-Inr.DN-Mosaic(3R)EtOH"|Genotype=="ci-Mosaic(3R)EtOH")
      #uncorrected P=0.2632683
 testdata<-subset(Ontogenetic.Allometry.Data, Genotype=="ci-Mosaic(3R)"|Genotype=="ci-Mosaic(3R)EtOH")
      #uncorrected P=0.382984
 testdata<-subset(Ontogenetic.Allometry.Data, Genotype=="ci-Inr.DN-Mosaic(3R)"|Genotype=="ci-Inr.DN-Mosaic(3R)EtOH")
      #uncorrected P=0.6586119
 
 #Fig 3A
 testdata<-subset(Ontogenetic.Allometry.Data, Genotype=="ci-EcR.RNAi(II)"|Genotype=="ci-WT17oC")
 #uncorrected P=0.02425538
 
 # commonslope test on test data-table comparing slopes across genotypes
 
 with(testdata, slope.com(lnAnterior, lnPosterior, Genotype, method="SMA"))
 
 