# A new implementation of network dissimilarity (betalink and related) --------------------

# idea: use vegdist to allow any dissimilarity index (implemented there, with the drawback of unusual naming)
  # but see below for which vegdist call to use for the poisot and bartomeus defaults
# implement poisots partitioning (with my modification) plus bartomeus partitioning (quantitative dissimilarity anyways)

# should explain my alternative partitioning concept in detail in the help file
# mention the conceptual problem of quantitative dissimilarity vs assumption "shared species subweb means no species dissimilarity"

#!!! webs2array is already present in bipartite, even with a vegdist example and Poisot2012 reference
  # if possible, change that function that it ALSO accepts a single argument that is already a list
  # if that's done, replace all weblist2array by webs2array
  # also, need to move the warning about only 2 webs to within the betalinkr-Function

# for bipartite
  # add test data to testfile
  # write rd files (use a stump from a simple function)
  # one file per function (use aliases for rd files)

# more Notes:
# can I handle multisite-comparisons in the partitioning? Probably not (or too complicated); make multiple comparisons by a metafunction
# allow for array input (default), list input (wishlist) and two single webs [or should that be a list, so that order of arguments is easier?]
# add Nacho's 2nd partitioning!?
  # not so sure what it means here / why it would be needed; and how it really fits to the rest of the framework
  # note that for the quantitative case, he still calls it BRich, but it should be BAbun (dissimilarity explained by abundance differences)
  # (nb: how does this partition differ from dissimilarity based on proportions, that also removes abundance differences?)
  # if I include this, then it's actually easier if index is newly developed (see next note); b_replacement = 2*min(B,C) / denom, b_rich or b_abun = abs(B-C) / denom
# develop new approach, calculating dissimilarities not with external function but internally
  # following Legendre 2014 (App), use only Jaccard/Ruzicka and Sorensen/"bray"
  # partition only numerator parts (number of different links, amounts of link differences), as consistent with set theory (cf. discussion with Benoit)
  # then use same denominator for all components
  # note that here (how I implemented it), frequency-changes of shared interactions are also included in rewiring! thus, quantitative version is conceptually different from function.dist="vegdist"! (and Bartomeus, I assume)
# for helper functions
  # "zero/forbidden links" should be removed from the output again to reduce storage size, for many large networks, the full array / linkmatrix might have to be avoided in the first place
  # consider adding an option "edgelist" to frame2webs (I guess that's only to handle large datasets efficiently; not sure right now what it's for otherwise; igraph can also convert a one-mode (square) adjacency matrix)



library(bipartite)  # includes vegan

# FUNCTION to make the needed array from a weblist -----
weblist2array <- function(weblist){
  if (length(weblist)!=2) warning("currently only accepting a list of two webs; function may break")
  rows.both <- sort(unique(c(rownames(weblist[[1]]), rownames(weblist[[2]]))))
  cols.both <- sort(unique(c(colnames(weblist[[1]]), colnames(weblist[[2]]))))
  webarray <- array(0, dim=c(length(rows.both), length(cols.both), 2))
  dimnames(webarray) <- list(rows.both, cols.both, names(weblist))
  for (i in 1:2) webarray[rownames(weblist[[i]]), colnames(weblist[[i]]), i] <- weblist[[i]]
  return(webarray)
}  
# weblist2array(list(Safariland, vazarr))  
# testarray <- weblist2array(list(Safariland, vazarr))  
    
# FUNCTION to turn an array with sites as third dimension into a siteXlink matrix (for dissim. calculations) ----------
  # note that this conversion is not needed when using my own 2-vector-function bcdist (and related)
  # also note that in the CA2002 project I create the linkmx directly from the dataframe (similar to igraph conversion used by Poisot)
  # names are optional (already taken care of in frame2webs, and new weblist2array)
array2linkmx <- function(webarray){
  # assumes the third dimension is the webID in the array, which will be the first dimension in the linkmx (vegan community matrix)
  # to easily fill the array into the matrix, first creating the transpose of linkmx
  linkmx.transp <- matrix(nrow = dim(webarray)[1] *  dim(webarray)[2], ncol = dim(webarray)[3])
    colnames(linkmx.transp) <- dimnames(webarray)[[3]] # use webIDs as names
  linkmx.transp[] <- webarray
    rownames(linkmx.transp) <- as.vector(outer(dimnames(webarray)[[1]], dimnames(webarray)[[2]], paste, sep="__")) # give names to links
  return(t(linkmx.transp))
}
# linkmx <- array2linkmx(webarray)  # example how to use


# main FUNCTION betalinkr --------------------
betalinkr <- function(webarray, method.dist = "bray", binary=TRUE, partition.osst="poisot", proportions=FALSE, function.dist="vegdist", distofempty="na"){
  # first, setting all the defaults to Poisot (for comparisons), later change to improved values
  # function arguments
    # webarray  input data, an array with 3dim, third dimension has length 2 and separates the two webs to compare
    # method.dist  the dissimilarity index, passed to "method" of either vegdist or betadiver (see there for naming)
    # binary  should binary data or quantitative data be used (i.e., quantitative or binary versions of dissimilarity indices)
    # partition.osst  "poisot" for the original ST=WN-OS, and "corrected" for the revised approach that calculates ST.raw directly and then adjusts OS and ST (OR the new alternative partitioning approach)
    # proportions  should data be standardized to proportions before calculating quantitative dissimilarity metrics?
    # function.dist  vegan package function to use for calculating dissimilarity (beta diversity), either "vegdist" or "betadiver" (the last gives compatibility to the 24 numbered indices of Koleff et al. 2003, but only binary)
      # a third option now also available! "new", which is a newly developed partitioning method, implemented incl. the dissimilarity index
    # distofempty "zero" or "na" how should dissimilarity be defined when there are either no links to use for b_os (i.e. only links involving species just found in one of the 2 webs) or no links to use for b_st (i.e. only "rewiring links" present); "zero" is appropriate when interested in the contribution of components b_os and b_st, whereas "na" is appropriate when b_os should be interpreted separately as dissimilarity of shared species subwebs

  # alternative subsets of linkmx (partitioning ST and OS) --
    # my approach here is to set the species/links to zero instead of excluding them (which will be done in the index calculation anyways)
    # this makes it easier to match species / links (even without names)
  
  if (class(webarray)=="list") {webarray <- weblist2array(webarray)}
  
  if (partition.osst=="poisot" & function.dist=="new") warning("poisot-style partitioning not compatible with function.dist='new'")
  
  # for "shared species subweb", set non-shared species links to zero
  array.sharedsp <- webarray
  array.sharedsp[rowSums(apply(webarray, MARGIN=c(1,3), sum)>0)!=2, , ] <- 0
  array.sharedsp[, rowSums(apply(webarray, MARGIN=c(2,3), sum)>0)!=2, ] <- 0

  # all links
  linkmx <- array2linkmx(webarray)

  # only links of shared species
  linkmx.sharedsp <- array2linkmx(array.sharedsp)  

  # the following block is NOT needed for function.dist=="new" (but cannot leave out unless changing the structure of "proportions" option below)
  # shared links of shared species (only LINKS occurring in both sites)
  linkmx.sharedli <- linkmx   
  linkmx.sharedli[, colSums(linkmx.sharedli>0)==1] <- 0
  # varying links of shared species
  linkmx.rewiring <- linkmx.sharedsp - linkmx.sharedli
  linkmx.RewSha <- linkmx.rewiring + linkmx.sharedli  # all links excluding those from unique species
  # links of non-shared species
  linkmx.uniquesp <- linkmx - linkmx.sharedsp
  linkmx.UniSha <- linkmx.uniquesp + linkmx.sharedli # all links excluding rewiring links
  # species community matrix (combining upper and lower of bipartite web)
  specmx.lower <- apply(webarray, c(3,1), sum)
  specmx.higher <- apply(webarray, c(3,2), sum)
  specmx.all <- cbind(specmx.lower, specmx.higher)  # e.g. sites X (plants, pollinators)
  
  # standardizing to proportions if wanted
  if (proportions){
    if (binary){warning("standardizing to proportions for binary index; do you really want to do this?!?")}
    specmx.all <- decostand(specmx.all, method="total")
    linkmx <- decostand(linkmx, method="total")
    linkmx.sharedsp <- decostand(linkmx.sharedsp, method="total")
    linkmx.RewSha <- decostand(linkmx.RewSha, method="total")
    linkmx.UniSha <- decostand(linkmx.UniSha, method="total")
  }
  
  # calculating dissimilarity / the betalink components --
  if (function.dist=="vegdist"){
    b_s <- vegdist(specmx.all, method=method.dist, binary=binary) # "S"
    b_wn <- vegdist(linkmx, method=method.dist, binary=binary) # "WN"
    if (distofempty=="zero" & any(rowSums(linkmx.RewSha)==0)){  # adjustment (set to the conceptually correct value); avoids warning
      b_os.raw <- b_wn
      b_os.raw[] <- 0
    } else {
      b_os.raw <- vegdist(linkmx.RewSha, method=method.dist, binary=binary) # "OS"
    }
    if (distofempty=="zero" & any(rowSums(linkmx.UniSha)==0)){ # adjustment (set to the conceptually correct value); avoids warning
      b_st.raw <- b_wn
      b_st.raw[] <- 0
    } else {
      b_st.raw <- vegdist(linkmx.UniSha, method=method.dist, binary=binary) # "ST"
    }
  } 
  if (function.dist=="betadiver"){
    if (binary==FALSE) {
      warning("betadiver only uses binary data; for quantitative indices use vegdist")
    } else {
      b_s <- betadiver(specmx.all, method=method.dist) # "S"
      b_wn <- betadiver(linkmx, method=method.dist) # "WN"
      if (distofempty=="zero" & any(rowSums(linkmx.RewSha)==0)){  # adjustment (set to the conceptually correct value); avoids warning
        b_os.raw <- b_wn
        b_os.raw[] <- 0
      } else {
        b_os.raw <- betadiver(linkmx.RewSha, method=method.dist) # "OS"
      }
      if (distofempty=="zero" & any(rowSums(linkmx.UniSha)==0)){ # adjustment (set to the conceptually correct value); avoids warning
        b_st.raw <- b_wn
        b_st.raw[] <- 0
      } else {
        b_st.raw <- betadiver(linkmx.UniSha, method=method.dist) # "ST"
      }
    }
  }
  if (function.dist=="new"){
    # here, method must be explicitly specified as one of Sorensen or Jaccard
      # quantitative equivalents available with binary=F; I actually use the quantitative formulae
    if (binary==TRUE){
      linkmx <- (linkmx > 0)
      linkmx.sharedsp <- (linkmx.sharedsp >0)
    }
    # A, B, C follows Legendre 2014
    # alle Komponenten berechnen aus linkmx, linkmx.sharedli, linkmx.rewiring and linkmx.uniquesp 
      # (do I need all 4 matrices?; maybe linkmx.sharedsp instead?)
    pmins <- pmin(linkmx[1,], linkmx[2,])
    A <- sum(pmins)
    B.tot <- sum(linkmx[1,] - pmins)
    C.tot <- sum(linkmx[2,] - pmins)
    B.rew <- sum(linkmx.sharedsp[1,] - pmin(linkmx.sharedsp[1,], linkmx.sharedsp[2,])) # note frequency-changes of shared interactions also included in rewiring here!
    B.uni <- B.tot - B.rew  # here it works with subtraction
    C.rew <- sum(linkmx.sharedsp[2,] - pmin(linkmx.sharedsp[1,], linkmx.sharedsp[2,]))
    C.uni <- C.tot - C.rew
    if (method.dist == "Sorensen"){denominator <- 2*A + B.tot + C.tot}
    if (method.dist == "Jaccard"){denominator <- A + B.tot + C.tot}
    b_wn <- (B.tot + C.tot) / denominator
    b_os <- (B.rew + C.rew) / denominator
    b_st <- (B.uni + C.uni) / denominator 
    b_s <- vegdist(specmx.all, method=switch(method.dist, "Jaccard"="jaccard", "Sorensen"="bray"), binary=binary) # "S"
  }
  
  # further calculations for the partitioning
  if(function.dist!="new"){
    b_st.minus <- b_wn - b_os.raw
    b_os <- b_os.raw * b_wn / (b_os.raw + b_st.raw)  # the new correction I propose to apply
    b_st <- b_st.raw * b_wn / (b_os.raw + b_st.raw)  # equivivalently for st
  }

  # output
  #!! probably change to a vector, e.g. using unlist!!
  if (partition.osst=="poisot") return(list(S=b_s, OS=b_os.raw, WN=b_wn, ST=b_st.minus))
  if (partition.osst=="corrected") return(list(S=b_s, OS=b_os, WN=b_wn, ST=b_st))
}


