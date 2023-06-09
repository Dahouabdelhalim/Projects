## This R script calculates nucleotide diversity using fasta files 
## that are provided by Stacks populations function with --fasta-samples option 
## FASTA should not hold more than one chromosomes

## Author: ITO Tsuyoshi
## Date: 2021-01-24
## E-mail: ito.tsuyoshi.3a@kyoto-u.ac.jp


library(ape)
library(tidyverse)
library(purrr)
library(pegas)
library(snowfall)

# Path to the directories of fasta files 
# (the outputs of Stacks populations function with --fasta-samples option)
dir_fasta <- ""
  
# Tab-delimited text file of id, sex (female or male), population, and cluster (if any)
file_list <- ""

# Name of Y and X chromosomes
chrname_y <- "NC_027914.1"
chrname_x <- "NC_041774.1"

# Chromosomes
vec_chromosome <- list.files(dir_fasta)

# Read fasta
list_fasta <- list()
for (i in 1:length(vec_chromosome)){
  list_fasta[[i]] <- read.FASTA(file = paste(dir_fasta, 
                                             vec_chromosome[i], 
                                             "/populations.samples.fa", 
                                             sep = ""))
}
rm(i)

# Population
vec_population <- c(
  "Shimokita", 
  "Yamagata", 
  "Gunma",
  "Shiga",
  "Arashiyama",
  "Takahama",
  "Wakasa",
  "Kochi",
  "Kojima",
  "Yakushima",
  "Taiwanese",
  "ChineseRhesus",
  "IndianRhesus",
  "Cynomolgus"
)

## cluster
#vec_cluster <- c(
#  "West",
#  "East",
#  "Yakushima"
#)

# Sample list
tb_sample <- read_delim(file = file_list, 
                        delim="\\t",
                        col_names = c("sample_id", "sex", "population")
                        #col_names = c("sample_id", "sex", "population", "cluster")
                        )

# Function of calculating pi
calc_pi <- function(fasta){
  
  # Extract sequence information
  fasta_names <- str_split(names(fasta), "_Sample_|_Locus_|_Allele_|;| \\\\[|,| ")
  
  tb_seq <- tibble(
    locus_id = map_chr(fasta_names, 3),
    sample_id = map_chr(fasta_names, 5),
    allele = map_chr(fasta_names, 4),
    chromosome = map_chr(fasta_names, 7),
    position = map_chr(fasta_names, 9)
  ) %>% 
    left_join(
      tb_sample,
      by = "sample_id")
  
  # Extract chromosome name
  chr <- unique(tb_seq$chromosome)
  
  # Check whether fasta contains more than one chromosomes
  if (length(chr) > 1){
    stop("Fasta contains more than one chromosomes")
  }
  
  # Check whether Y chromosome fasta contains female samples
  if (chr == chrname_y && "female" %in% tb_seq$sex){
    stop("Y chromosome fasta contains female samples")
  }
  
  n_seq_het <- NaN
  n_seq_male <- NaN
  n_seq_female <- NaN
  # For sex chromosomes, check the heterogeneity of sequence pairs
  if (chr == chrname_y || chr == chrname_x){
    # Split by sex
    index_female <- which(tb_seq$sex == "female")
    index_male <- which(tb_seq$sex == "male")
    
    n_seq_female <- length(index_female)
    n_seq_male <- length(index_male)
    
    fasta_f <- fasta[index_female]
    fasta_m <- fasta[index_male]
    
    tb_seq_f <- tb_seq[index_female, ]
    tb_seq_m <- tb_seq[index_male, ]
    
    # Index of selected data
    # Check the heterogeneity of sequence pairs
    n_seq_het <- 0
    index_sel <- NULL
    
    for (i in seq(1, length(fasta_m), 2)) {
      
      fasta_pair <- fasta_m[c(i, i+1)]
      names(fasta_pair) <- NULL # For all.equal()
      fasta_pair <- as.matrix(fasta_pair)
      
      if (isTRUE(all.equal(fasta_pair[1, ], fasta_pair[2, ]))) {
        index_sel <- c(index_sel, i)
      }
      else{
        # Count heterogeneous sequence pairs
        n_seq_het <- n_seq_het + 1
      }
    }
  }
  
  # Selected loci
  # X chromosome
  if (chr == chrname_x){
    tb_seq_new <- rbind(tb_seq_m[index_sel, ], tb_seq_f)
    fasta_new <- c(fasta_m[index_sel], fasta_f)
  }
  # Y chromosome
  else if(chr == chrname_y){
    tb_seq_new <- tb_seq_m[index_sel,]
    fasta_new <- fasta_m[index_sel]
  }
  # Autosome
  else{
    tb_seq_new <- tb_seq
    fasta_new <- fasta
  }
  
  # Calculate pi
  list_pi <- list()
  n_loci_use <- 0
  for (locus in unique(tb_seq_new$locus_id)) {
    
    # Index of a target locus
    index_sub <- which(tb_seq_new$locus_id == locus)
    
    # Subset of tibble of sequence information
    tb_seq_sub <- tb_seq_new[index_sub, ]
    
    # Check missing rate for each population
    tb_missing <- tb_sample %>%
      # For Y chromosome, only males are extracted
      {if(chr == chrname_y) filter(., sex == "male") else .} %>%
      left_join(
        tibble(
          sample_id = tb_seq_sub$sample_id,
          allele = tb_seq_sub$allele,
          fasta = 1
        ),
        by = "sample_id"
      ) %>%
      replace_na(list(fasta = 0)) %>%
      group_by(population) %>%
      summarize(n = n(), n_fasta = sum(fasta), r_fasta = sum(fasta)/n())
    
    # Skip loci with # of sequences < 2 at least one population
    if (min(tb_missing[tb_missing$n > 1, ]$n_fasta) > 1){
      n_loci_use <- n_loci_use + 1
      
      # Selected sequence
      # Subset of fasta
      fasta_sub <- as.matrix(fasta_new[index_sub])
      
      vec_pi <- NULL
      for (population in vec_population){
        
        vec_pi <- c(vec_pi, 
                    nuc.div(fasta_sub[which(tb_seq_sub$population == population), ], 
                            variance = FALSE, 
                            # Delete the sites with missing data in a pairwise way
                            pairwise.deletion = TRUE)
        )
      }
      names(vec_pi) <- vec_population
      
      #for (cluster in vec_cluster){
      #  
      #  vec_pi <- c(vec_pi, 
      #              nuc.div(fasta_sub[which(tb_seq_sub$cluster == cluster), ], 
      #                      variance = FALSE, 
      #                      # Delete the sites with missing data in a pairwise way
      #                      pairwise.deletion = TRUE)
      #  )
      #}
      #names(vec_pi) <- vec_cluster
      
      list_pi <- append(list_pi, list(bind_rows(vec_pi)))
    }
  }
  
  # To tibble
  tb_pi <- bind_rows(list_pi)
  
  # Summary table
  tb_summary <- tibble(
    n_seq_male = n_seq_male,
    n_seq_female = n_seq_female,
    n_seq_het = n_seq_het,
    r_seq_het = n_seq_het/n_seq_male,
    n_loci = length(unique(tb_seq$locus_id)),
    n_loci_new = length(unique(tb_seq_new$locus_id)),
    n_loci_use = n_loci_use
  )
  
  print(paste(chr, " done", sep = ""))
  
  return(list(
    pi = tb_pi,
    summary = tb_summary))
  
}


# Run in parallel
sfInit(parallel = TRUE, cpus = 6)
sfLibrary(tidyverse)
sfLibrary(pegas)
sfExport(list = c("tb_sample", 
                  "vec_population",
                  #"vec_cluster",
                  "vec_chromosome", 
                  "calc_pi",
                  "chrname_y", 
                  "chrname_x")
         )

t <- proc.time()
results <- sfLapply(list_fasta, calc_pi)
proc.time() - t
rm(t)

sfRemoveAll()
sfStop()


# Summary
tb_pi_mean <- purrr::map_df(lapply(results, function(x) x[[1]]), 
                            function(x){summarize_all(x, mean, na.rm = TRUE)}) %>%
  bind_cols(
    tibble(
      chromosome_id = vec_chromosome,
      chromosome_type = c("Y", rep("A", 20), "X")
    ) 
  )
tb_pi_summary <- purrr::map_df(lapply(results, function(x) x[[2]]), 
                               function(x){summarize_all(x, mean, na.rm = TRUE)}) %>%
  bind_cols(
    tibble(
      chromosome_id = vec_chromosome,
      chromosome_type = c("Y", rep("A", 20), "X")
    ) 
  )
write_csv(tb_pi_mean, "tb_pi_mean.csv")
write_csv(tb_pi_summary, "tb_pi_summary.csv")

## Mean pi of autosome
#tb_pi_a <- tibble()
#for (i in 2:21){
#  tb_pi_a <- tb_pi_a %>% 
#    bind_rows(results[[i]]$pi)
#}
#
#mean_pi_cluster <- tb_pi_a %>%
#  summarize_all(mean)
