##### R script to reproduce all analyses in Phillips et al. (2023) #####
##### Note: DNA sequence alignments must in uppercase font and be exported in FASTA format from MEGA6 #####

library(VLF) # load VLF package

vlfFun <- function(x, p = 0.001, seqlength = 648, own = NULL) 
{
  species.names <- x[, 2]
  specimen.Number <- nrow(x)
  rownames(x) <- species.names
  Nuc.count <- count.function(x, specimen.Number, seqlength)
  frequency.matrix <- ffrequency.matrix.function(Nuc.count, 
                                                 seqlength)
  spec.freq <- specimen.frequencies(frequency.matrix, x, specimen.Number, 
                                    species.names, seqlength)
  nucleotide.modalSequence <- MODE(frequency.matrix, seqlength)
  first.modal.frequencies <- MODE.freq(frequency.matrix, seqlength)
  second.modal.frequencies <- MODE.second.freq(frequency.matrix, 
                                               seqlength)
  First_conserved_100 <- conservation_first(first.modal.frequencies, 
                                            1, seqlength)
  First_conserved_99.9 <- conservation_first(first.modal.frequencies, 
                                             (1 - p), seqlength)
  FirstAndSecond_conserved_99.9 <- conservation_two(first.modal.frequencies, 
                                                    second.modal.frequencies, (1 - p), seqlength)
  specimen_VLFcount <- VLF.count.spec(spec.freq, p, seqlength)
  position_VLFcount <- VLF.count.pos(spec.freq, p, seqlength)
  VLFconvert <- VLF.convert.matrix(x, spec.freq, p, seqlength)
  VLFnuc <- VLF.nucleotides(VLFconvert, x, seqlength)
  VLFreduced <- VLF.reduced(VLFnuc, specimen_VLFcount, seqlength)
  species <- separate(VLFreduced)
  singleAndShared <- find.singles(species, seqlength)
  if (is.null(own)) {
    foo <- list(modal = nucleotide.modalSequence, con100 = First_conserved_100, 
                conp = First_conserved_99.9, combine = FirstAndSecond_conserved_99.9, 
                specimen = specimen_VLFcount, position = position_VLFcount, 
                sas = singleAndShared, VLFmatrix = VLFreduced)
    class(foo) <- "vlf"
    foo
  }
  else {
    own.species.names <- own[, 2]
    own.specimen.Number <- nrow(own)
    rownames(own) <- own.species.names
    own.Nuc.count <- count.function(own, own.specimen.Number, seqlength)
    own.frequency.matrix <- ffrequency.matrix.function(own.Nuc.count, 
                                                       seqlength)
    ownspec.freq <- specimen.frequencies(frequency.matrix, 
                                         own, nrow(own), own[, 2], seqlength)
    own.nucleotide.modalSequence <- MODE(own.frequency.matrix, seqlength)
    own.first.modal.frequencies <- MODE.freq(own.frequency.matrix, seqlength)
    own.second.modal.frequencies <- MODE.second.freq(own.frequency.matrix, 
                                                     seqlength)
    own.First_conserved_100 <- conservation_first(own.first.modal.frequencies, 
                                                  1, seqlength)
    own.First_conserved_99.9 <- conservation_first(own.first.modal.frequencies, 
                                                   (1 - p), seqlength)
    own.FirstAndSecond_conserved_99.9 <- conservation_two(own.first.modal.frequencies, 
                                                          own.second.modal.frequencies, (1 - p), seqlength)
    ownspec.VLFcount <- VLF.count.spec(ownspec.freq, p, seqlength)
    ownpos.VLFcount <- VLF.count.pos(ownspec.freq, p, seqlength)
    own.VLFconvert <- VLF.convert.matrix(own, ownspec.freq, 
                                         p, seqlength)
    own.VLFnuc <- VLF.nucleotides(own.VLFconvert, own, seqlength)
    own.VLFreduced <- VLF.reduced(own.VLFnuc, ownspec.VLFcount, 
                                  seqlength)
    own.species <- separate(own.VLFreduced)
    own.singleAndShared <- find.singles(own.species, seqlength)
    foo <- list(own.modal = own.nucleotide.modalSequence, own.con100 = own.First_conserved_100, 
                own.conp = own.First_conserved_99.9, own.combine = own.FirstAndSecond_conserved_99.9, 
                own.sas = own.singleAndShared, own.VLFmatrix = own.VLFreduced, ownSpecCount = ownspec.VLFcount, 
                ownPosCount = ownpos.VLFcount, ownVLFMatrix = own.VLFnuc, 
                ownVLFreduced = own.VLFreduced)
    class(foo) <- "vlf"
    foo
  }
}

### VLF analysis on full birds dataset, N = 11333 barcodes ###

data(birds) # load birds nucleotide dataset

out_birds <- vlfFun(birds)

sum(out_birds$specimen) # total number of specimen VLFs - 771
sum(out_birds$position) # total number of positionsl VLFs - 771

length(which(out_birds$specimen != 0)) # total number of specimens with VLFs - 552
length(which(out_birds$position != 0)) # # total number of positions with VLFs - 241

which(out_birds$specimen != 0) # specimens with VLFs
which(out_birds$position != 0) # sites containing VLFs

rowSums(out_birds$sas) # number of singleton and shared VLFs

species.names <- birds[,2]
specimen.Number <- nrow(birds)
rownames(birds) <- species.names
Nuc.count <- count.function(birds, specimen.Number, seqlength = 648)
frequency.matrix <- ffrequency.matrix.function(Nuc.count, seqlength = 648)
birdSpec.freq <- specimen.frequencies(frequency.matrix, birds, specimen.Number, species.names, seqlength = 648)
Bird_specimen_VLFcount <- VLF.count.spec(birdSpec.freq, p = 0.001, seqlength = 648)
Bird_position_VLFcount <- VLF.count.pos(birdSpec.freq, p = 0.001, seqlength = 648)
bird_VLFconvert <- VLF.convert.matrix(birds, birdSpec.freq, p = 0.001, seqlength = 648)
bird_VLFnuc <- VLF.nucleotides(bird_VLFconvert, birds, seqlength = 648)
bird_VLFreduced <- VLF.reduced(bird_VLFnuc, Bird_specimen_VLFcount, seqlength = 648)
bird_species <- separate(bird_VLFreduced)
birds_singleAndShared <- find.singles(bird_species, seqlength = 648)

# singleton error rate
birds_single_error <- Error.Rate(birds_singleAndShared[1,], spec = specimen.Number, seqlength = 648) 
# shared error rate
birds_shared_error <- Error.Rate(birds_singleAndShared[2,], spec = specimen.Number, seqlength = 648)
# total error rate
birds_error_total <- Error.Rate(birds_singleAndShared[1,] + birds_singleAndShared[2,], spec = specimen.Number, seqlength = 648)

Decile.Plot(birds_singleAndShared, seqlength = 648) # plot decile plot
Sliding.Window(birds_singleAndShared, seqlength = 648, n = 30) # plot sliding window


### VLF analysis for Canada goose (Branta canadensis), N =  125 barcodes ###

x <- separate(birds) # partition into lists by species name

y <- x[[317]]

out_goose <- vlfFun(birds, own = y) # B. canadensis is list element 317

sum(out_goose$ownSpecCount) # total number of specimen VLFs - 27
sum(out_goose$ownPosCount) # total number of positionsl VLFs - 27

length(which(out_goose$ownSpecCount != 0)) # total number of specimens with VLFs - 18
length(which(out_goose$ownPosCount != 0)) # total number of positions with VLFs - 10

which(out_goose$ownSpecCount != 0) # specimens with VLFs 
which(out_goose$ownPosCount != 0) # sites containing VLFs

table(out_goose$ownSpecCount) # frequency of specimen VLFs
table(out_goose$ownPosCount) # frequency of positional VLFs

rowSums(out_goose$own.sas)

goose_single_error <- Error.Rate(out_goose$sas[1,], out_goose$sas[2,], spec = 125, seqlength = 648) 
goose_shared_error <- Error.Rate(out_goose$sas[2,], out_goose$sas[1,], spec = 125, seqlength = 648)
goose_totalerror <- Error.Rate(out_goose$sas[1,] + out_goose$sas[1,], spec = 125, seqlength = 648)

Decile.Plot(out_goose$own.sas, seqlength = 648) # plot decile plot
Sliding.Window(out_goose$own.sas, seqlength = 648, n = 30) # plot sliding window


### VLF analysis for fishes, N = 2371 barcodes ###

sequences <- fasta.read("/Users/jarrettphillips/desktop/VLF/fish_aligned.fas", seqlength = 652, pos1 = 1, pos2 = 2)

species.names <- sequences[,2]
specimen.Number <- nrow(sequences)
rownames(sequences) <- species.names
Nuc.count <- count.function(sequences, specimen.Number, seqlength = 652)
frequency.matrix <- ffrequency.matrix.function(Nuc.count, seqlength = 652)
spec.freq <- specimen.frequencies(frequency.matrix, sequences, specimen.Number, species.names, seqlength = 652)
specimen_VLFcount <- VLF.count.spec(spec.freq, p = 0.001, seqlength = 652)
position_VLFcount <- VLF.count.pos(spec.freq, p = 0.001, seqlength = 652)
VLFconvert <- VLF.convert.matrix(sequences, spec.freq, p = 0.001, seqlength = 652)
VLFnuc <- VLF.nucleotides(VLFconvert, sequences, seqlength = 652)
VLFnuc_reduced <- VLF.reduced(VLFnuc, specimen_VLFcount, seqlength = 652)
VLFsingletons <- find.singles(separate(VLFnuc_reduced), seqlength = 652)

rowSums(VLFsingletons) 

sum(specimen_VLFcount) # total number of specimen VLFs - 117
sum(position_VLFcount) # total number of positional VLFs - 117

length(which(specimen_VLFcount != 0)) # total number of specimens with vlfS - 58
length(which(position_VLFcount != 0)) # total number of positions with VLFs - 84

which(specimen_VLFcount != 0) # specimens with VLFs
which(position_VLFcount != 0) # sites containing VLFs

table(specimen_VLFcount) # frequency of specimen VLFs
table(position_VLFcount) # frequency of positional VLFs

# singleton error rate
single <- Error.Rate(VLFsingletons[1,], VLFsingletons[2,], spec = 2371, seqlength = 652)
# shared error rate
shared <- Error.Rate(VLFsingletons[2,], VLFsingletons[1,], spec = 2371, seqlength = 652)
# total error rate
total <- Error.Rate(VLFsingletons[1,] + VLFsingletons[2,], spec = 2371, seqlength = 652)

Decile.Plot(VLFsingletons, seqlength = 652) # plot decile plot
Sliding.Window(VLFsingletons, seqlength = 652) # plot sliding window


