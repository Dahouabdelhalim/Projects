## Identify References - Revised

# Identify references by assigning them to known sources

## configuration
#setwd("/home/fm/GlobalPathways/")
setwd("c:/Users/Felix/Desktop/GlobalPathways/")
data_path <- "./data/"
image_path <- "./images/"


library(stringdist)


############## Functions for reference identification
# (See next section for running the identification)



## Function OLA_dist

# Adapted Smith-Waterman algorithm for 'O'ptimal 'L'ocal sequence 'A'lignment.
# Calculates the score for the best "local" (end gaps in either sequence are free!) sequence alignment of two strings.
# A perfect local alignment means the shorter string aligns without error at some position in the longer string.
# In this case the optimal local alignment score is equal to the length of the shorter string.

# x, y: two strings to be compared

# Returns a distance, computed by subtracting the number of matched characters in the best local alignment from the length of the shorter string.

OLA_dist <- function(x,y){
  # If x or y are empty strings or NA, return NA
  if((is.na(x) | x == "") | (is.na(y) | y == "")){
    return(NA)
  }
  nx <- nchar(x)
  ny <- nchar(y)
  match_credit <- 1
  gap_cost <- .95 # cost for internal gaps
  mismatch_cost <- .9 # cost for internal mismatched characters
  score_max <- 0 # initialize alignment score
  score_max_matches <- 0 # initialize number of matched characters
  score_max_ij <- c(0,0)
  # initialize the dynamic programming matrices with zero (zero for free end gaps)
  M <- matrix(0, nrow = nx+1, ncol = ny+1)
  M2 <- matrix(0, nrow = nx+1, ncol = ny+1)
  # Run the Smith-Waterman algorithm and remember the best local alignment score so far.
  # Also, remember the number of matched characters in the best local alignment so far.
  for(i in 1:nx){
    for(j in 1:ny){
      # calculate alignment scores in the dynamic programming matrix (see Smith-Waterman alg.)
      if(substr(x,i,i) == substr(y,j,j)){ match <- match_credit }else{ match <- -mismatch_cost }
      M[i+1,j+1] <- max(
        M[i,j] + match,
        M[i,j+1] - gap_cost,
        M[i+1,j] - gap_cost,
        0
      )
      # use a second matrix to store the number of perfectly matched characters in each possible alignment
      if(M[i+1,j+1] == 0){
        M2[i+1,j+1] <- 0
      }else if(M[i+1,j+1] == M[i,j] + match_credit){
        M2[i+1,j+1] <- max(M2[i,j] + 1, M2[i+1,j], M2[i,j+1])
      }else{
        M2[i+1,j+1] <- max(M2[i,j], M2[i+1,j], M2[i,j+1])
      }
      # store the the maximum score and the number of matches in the corresponding alignment
      if(M[i+1,j+1] > score_max){
        score_max <- M[i+1,j+1]
        score_max_matches <- M2[i+1,j+1]
        score_max_ij <- c(i+1,j+1)
      }
    }
  }
  # Traceback to find number of internal gaps for each string in the alignment
  s <- score_max
  i <- score_max_ij[1]; j <- score_max_ij[2]
  gx <- 0; gy <- 0
  while(s > 0){
    if(M[i-1,j-1] >= M[i-1,j] & M[i-1,j-1] >= M[i,j-1]){
      i <- i-1; j <- j-1
    }else if(M[i-1,j] >= M[i-1,j-1] & M[i-1,j] >= M[i,j-1]){
      i <- i-1; gy <- gy+1 
    }else{
      j <- j-1; gx <- gx+1
    }
    s <- M[i,j]
  }
  # Return length of the shorter string, minus the number of matches in the optimal alignment, as a distance.
  return(min(nx+gx,ny+gy) - score_max_matches)
}


## Function findBestCandidate

# Performs different forms of approximate string matching for a reference with a number of source candidates.
# Returns a string with the source ID of the best matching source candidate together with keywords
# which signify whether the reference could be assigned to a candidate and the distance between the reference
# and the best candidate.
# Returns a string with multiple source IDs, if the best candidate is ambiguous due to small distances for 
# multiple source candidates.

# "candidates" is a data frame of up to "n.candidates" many source candidates (rows) 
# that should be compared to the reference.

# "reference" is a single row of the references data frame.

# "threshold" is the maximum distance threshold up to which a reference is assigned to a source candidate.

# "missing_field_dists" is a named vector of default distances to apply if source candidates have no entry
#  or NA in a specific field.

# "field_weights" is a named vector of multiplicative weights which are applied to the calculated distances for each field.

# "min_separation_threshold" is the minimum difference in source candidate - reference distances 
#  for which one of two similarly well fitting candidates can be selected as the real source. 

# If multiple accepted candidates have distances closer to each other than "min_separation_threshold",
# the reference is marked as ambiguous with the source IDs of all respective candidates.
# Otherwise the best accepted candidate is assigned to the reference. 

# All parameters are passed from the main function "identifyReferences". See there for details.

findBestCandidate <- function(candidates, reference, threshold, 
                              missing_field_dists, field_weights,
                              min_separation_threshold=2, 
                              f.page_allowed_numeric_distance=10, 
                              f.page_extra_distance=3){
  
  # transform reference from data frame (single row) to character vector
  fieldnames <- names(reference)
  reference <- as.character(reference)
  names(reference) <- fieldnames
  
  nc <- nrow(candidates)
  
  # check if any candidates have at least (number of fields - 1) non-empty identical fields
  # to the reference. These will be assigned a distance of 0 by default (see end of function).
  enough_id_fields <- rep(FALSE, nc)
  ref_filled <- !(is.na(reference) | reference == "")
  if(sum(! ref_filled) <= 1){
    nonempty_ref_fields <- fieldnames[ref_filled]
    candidate_id_fields <- apply(candidates, 1, function(x){
      r <- reference[nonempty_ref_fields] == x[nonempty_ref_fields]
      r[is.na(r)] <- FALSE
      return(sum(r))
    })
    enough_id_fields <- candidate_id_fields >= (length(fieldnames) - 1)
  }
  
  # Distances for fields other than the specific fields in our reference data 
  # get the standard Dammerau-Levenshtein-distance
  extra_dist <- rep(0, nc)
  if(any(! fieldnames %in% c("DOI","aut.1","year","pub.name","vol","f.page"))){
    
    extra_fields <- fieldnames[! fieldnames %in% c("DOI","aut.1","year","pub.name","vol","f.page")]
    
    for(field in extra_fields){
      dist <- rep(0, nc)
      has_field <- !is.na(candidates[,field]) & candidates[,field] != ""
      if(!is.na(reference[field]) & reference[field] != ""){
        dist[has_field] <- stringdist(reference[field], candidates[has_field,field], method="dl") * field_weights
        dist[!has_field] <- missing_field_dists[field]
      }
      extra_dist <- extra_dist + dist
      
    }
  }
  
  # DOI field distance
  # Assign distance 0 if the shorter DOI between reference and a source candidate is matches exectly
  # in the longer DOI. Assign the field weight as distance otherwise.
  # If either reference or source candidate have no DOI, also assign distance 0.
  doi_dist <- rep(0, nc)
  if("DOI" %in% fieldnames){
    has_doi <- !is.na(candidates[,"DOI"]) & (candidates[,"DOI"] != "")
    if(!is.na(reference["DOI"]) & reference["DOI"] != "" & any(has_doi)){
      n_ref_DOI <- nchar(reference["DOI"])
      doi_dist[has_doi] <- sapply(candidates[has_doi,"DOI"], function(x){
        if(nchar(x) < n_ref_DOI){
          r <- as.numeric(! grepl(x, reference["DOI"], fixed=TRUE)) * field_weights["DOI"]
        }else{
          r <- as.numeric(! grepl(reference["DOI"], x, fixed=TRUE)) * field_weights["DOI"]
        }
        return(r)
      })
    }
  }
  # aut.1 field distance
  # also use OLA distance (see "OLA_dist" function above)
  # Apply missing aut.1 distance if only the reference has aut.1 or author.1 and the candidate has neither. 
  aut.1_dist <- rep(0, nc)
  if("aut.1" %in% fieldnames){
    has_aut.1 <- !is.na(candidates[,"aut.1"]) & candidates[,"aut.1"] != ""
    has_author.1 <- !is.na(candidates[,"author.1"]) & candidates[,"author.1"] != ""
    has_aut <- has_aut.1 | has_author.1
    if(!is.na(reference["aut.1"]) & reference["aut.1"] != ""){
      
      aut.1_c_dist <- sapply(candidates[has_aut,"aut.1"], function(x){OLA_dist(reference["aut.1"], x)})
      
      author.1_c_dist <- sapply(candidates[has_aut,"author.1"], function(x){OLA_dist(reference["aut.1"], x)})
      
      aut.1_c_dist[is.na(aut.1_c_dist)] <- 1000
      author.1_c_dist[is.na(author.1_c_dist)] <- 1000
      
      nchar_aut.1 <- nchar(candidates[has_aut,"aut.1"])
      nchar_aut.1[is.na(nchar_aut.1)] <- 0
      nchar_author.1 <- nchar(candidates[has_aut,"author.1"])
      nchar_author.1[is.na(nchar_author.1)] <- 0
      
      weights_aut.1 <- as.numeric(nchar_aut.1 <= 6) + 1
      weights_author.1 <- as.numeric(nchar_author.1 <= 6) + 1
      weights_ref <- rep(as.numeric(nchar(reference["aut.1"]) <= 6) + 1, sum(has_aut))
      
      aut.1_c_dist <- aut.1_c_dist * apply(cbind(weights_aut.1, weights_ref), 1, max)
      author.1_c_dist <- author.1_c_dist * apply(cbind(weights_author.1, weights_ref), 1, max)
      
      aut.1_dist[has_aut] <- apply(cbind(aut.1_c_dist, author.1_c_dist), 1, min) * field_weights["aut.1"]
      aut.1_dist[!has_aut] <- missing_field_dists["aut.1"]
    }
  }
  # year field distance
  # use Dammerau-Levenshtein distance
  # Apply missing dist. if only the reference has year
  year_dist <- rep(0, nc)
  if("year" %in% fieldnames){
    has_year <- !is.na(candidates[,"year"]) & candidates[,"year"] != ""
    if(!is.na(reference["year"]) & reference["year"] != ""){
      year_dist[has_year] <- stringdist(reference["year"], candidates[has_year,"year"], method="dl") * field_weights["year"]
      year_dist[!has_year] <- missing_field_dists["year"]
    }
  }
  # pub.name field distance
  # Use OSA distance and substract the difference in string length between reference and candidate 
  # if both have length >= 5. Otherwise use unmodified OSA distance.
  # Apply missing dist. if only the reference has pub.name/full.pub
  pub.name_dist <- rep(0, nc)
  if("pub.name" %in% fieldnames){
    has_pub.name <- (!is.na(candidates[,"pub.name"]) & candidates[,"pub.name"] != "") | (!is.na(candidates[,"full.pub"]) & candidates[,"full.pub"] != "")
    if(!is.na(reference["pub.name"]) & reference["pub.name"] != ""){
      if(sum(has_pub.name) > 0){
        pub_dists <- rbind(stringdist(reference["pub.name"], candidates[has_pub.name,"pub.name"], method="osa"), 
                           stringdist(reference["pub.name"], candidates[has_pub.name,"full.pub"], method="osa")) 
        nr <- nchar(reference["pub.name"])
        if(nr >= 5){
          pub.name_length <- sapply(candidates[has_pub.name,"pub.name"], nchar)
          pub.name_length[is.na(pub.name_length)] <- 0
          pub_dists[1, pub.name_length >= 5] <- pub_dists[1, pub.name_length >= 5] - abs(pub.name_length[pub.name_length >= 5] - nr)
          full.pub_length <- sapply(candidates[has_pub.name,"full.pub"], nchar)
          full.pub_length[is.na(full.pub_length)] <- 0
          pub_dists[2, full.pub_length >= 5] <- pub_dists[2, full.pub_length >= 5] - abs(full.pub_length[full.pub_length >= 5] - nr)
        }
        pub_dists <- apply(pub_dists, 2, function(c){c[is.na(c)] <- c[!is.na(c)]; return(c)})
        pub.name_dist[has_pub.name] <- apply(pub_dists, 2, min) * field_weights["pub.name"]
      }
      pub.name_dist[!has_pub.name] <- missing_field_dists["pub.name"]
    }
  }
  
  # vol field distance
  # also use Dammerau-Levenshtein
  # Apply missing dist. if only reference or only source candidate have vol
  vol_dist <- rep(0, nc)
  if("vol" %in% fieldnames){
    has_vol <- !is.na(candidates[,"vol"]) & candidates[,"vol"] != ""
    if(!is.na(reference["vol"]) & reference["vol"] != ""){
      vol_dist[has_vol] <- stringdist(reference["vol"], candidates[has_vol,"vol"], method="dl") * field_weights["vol"]
      vol_dist[!has_vol] <- missing_field_dists["vol"]
    }
  }
  # f.page field distance
  # also use Dammerau-Levenshtein
  # Apply missing dist. if only reference or only source candidate have f.page
  f.page_dist <- rep(0, nc)
  if("f.page" %in% fieldnames){
    has_f.page <- !is.na(candidates[,"f.page"]) & candidates[,"f.page"] != ""
    if(!is.na(reference["f.page"]) & reference["f.page"] != ""){
      f.page_dist[!has_f.page] <- missing_field_dists["f.page"]
      # If entries can be converted to numeric, compute the numeric page number difference.
      # If it is higher than "f.page_allowed_numeric_distance", add extra distance.
      # The extra difference is NOT multiplied by the field weight (hopefully more intuitive since it is given separately to the field_weights vector)
      if(!is.na(as.numeric(reference["f.page"]))){
        ref_num_page <- as.numeric(reference["f.page"])
        num_page <- suppressWarnings(as.numeric(candidates[has_f.page,"f.page"]))
        #set entries that cannot be converted to numeric to -(allowed distance + 1)
        #so that they always get the extra distance applied. 
        num_page[is.na(num_page)] <- -(f.page_allowed_numeric_distance + 1) 
        small_num_diff <- abs(num_page - ref_num_page) != 0 & abs(num_page - ref_num_page) <= f.page_allowed_numeric_distance
        large_num_diff <- abs(num_page - ref_num_page) > f.page_allowed_numeric_distance
        f.page_dist[has_f.page][small_num_diff] <- field_weights["f.page"]
        f.page_dist[has_f.page][large_num_diff] <- field_weights["f.page"] + f.page_extra_distance
      }else{
        f.page_dist[has_f.page & !(candidates[has_f.page,"f.page"] == reference["f.page"])] <- field_weights["f.page"] + f.page_extra_distance
      }
      # for f.page we also add the missing field distance, if the reference has no entry but the source candidate has an entry.
      # This is different from the other fields, where a missing reference entry is never punished!
    }else{ 
      f.page_dist[has_f.page] <- missing_field_dists["f.page"]
    }
  }
  # construct an overall distance
  all_dist <- extra_dist + doi_dist + aut.1_dist + year_dist + pub.name_dist + vol_dist + f.page_dist
  
  # modifications
  # If volume and page BOTH have differences for a candidate, increase the distance by 1.
  all_dist[vol_dist != 0 & f.page_dist != 0] <- all_dist[vol_dist != 0 & f.page_dist != 0] + 1
  
  if(any(is.na(all_dist))){
    my_names <- c("extra_dist","doi_dist","aut.1_dist","year_dist","pub.name_dist","vol_dist","f.page_dist")
    my_dists <- list(extra_dist,doi_dist,aut.1_dist,year_dist,pub.name_dist,vol_dist,f.page_dist)
    for(i in 1:length(my_dists)){
      print(paste("NA in dist: ", my_names[i], " ", any(is.na(my_dists[[i]]))))
    }
  }
  
  # add candidates with enough identical fields as candidates with distance zero
  all_dist[enough_id_fields] <- 0
  
  # valid candidates
  valid <- which(all_dist <= threshold)
  
  # Select candidate with minimum distance to the reference
  # "best_candidate" should hold a string with either the distance to the best match and its source id  (i.e. "d12 17493")
  # or "overTresh" with distance to the best match and the new source id of the newly created source, if no match is under the distance threshold.
  # If the best candidates are too close to each other in distance, it should hold "amb" with the minimum distance 
  # and the source idsall ambiguous source candidates (i.e. "amb d5 123987 83476 812312").
  
  # initialize for the case that no candidate is under the distance threshold (a new source is and source id is added in the main function)
  closest <- which(all_dist == min(all_dist))[1]
  best_candidate <- paste(c("overThresh ", "d", min(all_dist), " ", candidates[closest,"sou.id"]), collapse="") 
  
  # if there is exactly one valid candidate, give its distance and source id.
  if(length(valid) == 1){
    best_candidate <- paste(c("d", min(all_dist)," ", candidates[valid,"sou.id"]), collapse="")
    # if there are multiple valid candidates:
  }else if(length(valid) > 1){
    # if the best candidates are closer to each other than "min separation threshold", give the smallest distance and all ambiguous source ids
    if( any( (sort(all_dist[valid])[-1] - min(all_dist)) < min_separation_threshold ) ){
      amb_valid_dist_indices <- which((all_dist[valid] - min(all_dist[valid])) < min_separation_threshold)
      amb_sou.ids <- candidates$sou.id[valid][order(all_dist[valid][amb_valid_dist_indices])]
      amb_dist <- paste(c("amb ","d", min(all_dist)), collapse="")
      best_candidate <- paste(c(amb_dist, amb_sou.ids), collapse=" ")
      # if not, just give the best distance with corresponding source id
    }else{
      best_candidate <- paste(c("d", min(all_dist)," ", candidates[which(all_dist == min(all_dist)),"sou.id"]), collapse="")
    }
  }
  return(best_candidate)
}


## Function buildSourceRefMatrix

# "sources" is the source table.
# "reference" is a character vector of all relevant reference fields for one reference.

# Builds a numerical matrix of dimensons nrow(sources) x length(fields)+1.
# The first column contains an index
# Other matrix fields [i,j] are either 0 or 1, depending on whether source i matches 
# reference field j-1 exactly (-1 because of index column). 

buildSourceRefMatrix <- function(sources, reference){
  
  fields <- names(reference)
  ns <- nrow(sources)
  nf <- length(fields)
  s.matrix <- matrix(FALSE, nrow=ns, ncol=nf+1)
  s.matrix[,1] <- 1:ns
  
  for(i in 1:nf){
    
    field <- fields[i]
    rf <- as.character(reference[field])
    
    if(!is.na(rf) & rf != ""){
      
      if(field == "DOI"){
        s.matrix[, i+1] <- grepl(rf, sources[,field], fixed=TRUE)
      }else if(field == "aut.1"){
        s.matrix[, i+1] <- grepl(rf, sources[,field], fixed=TRUE) | grepl(rf, sources$author.1, fixed=TRUE)
      }else if(field == "pub.name"){
        s.matrix[, i+1] <- grepl(rf, sources[,field], fixed=TRUE) | grepl(rf, sources$full.pub, fixed=TRUE)
      }else{
        s.matrix[, i+1] <- rf == sources[,field]
      }  
      
    }else{
      s.matrix[, i+1] <- FALSE
    }
  }
  
  return(s.matrix)
  
}# end buildSourceRefMatrix




## Function identifyReferences:

# This function takes 2 data frames "references" and "sources" and finds up to
# "n.candidates" many candidates for each row of references in sources by exact string matching (fast).
# It then gives these source candidates and the reference to function "findBestCandidate"
# to compare them to the reference with approximate string matching (slow).

# If a good enough source candidate is found, the reference is marked with the candidate's sou.id (source identifier).
# If no valid candidate is found, the reference is added to sources with a new sou.id (and marked with that id in references).
# If there are multiple good candidates, too close in distance, the reference is marked as ambiguous with all possible sou.ids

# "ordered.subsets" is a matrix of boolean vectors, where each component (row) 
# corresponds to a column (e.g. DOI, aut.1, year) of both dataframes, such that
# each vector (column of ordered.subsets) corresponds to a subset of data frame columns.
# The rownames of "ordered.subsets" must be column names of "references" and "sources".
# The purpose of these field sets is to have an ordered list of field combinations to
# find the best source candidates: Candidates that exactly match the reference in all fields of the field combination,
# e.g. a candidate that matches the reference in author, year and DOI is better than one that matches in 
# author, volume and page.
# Going through "ordered.subsets" from left to right, the function selects candidates
# from "sources" which are equal to the processed "references" row in the specified fields/columns,
# until up to "n.candidates" candidates are found.

# "ref_set" allows the user to pass a selection of reference indices (if not the whole "references" data frame
# should be used). Passing the selection as "ref_set" is preferable to restricting "references" directly,
# because the result candidates will show the correct "ref.nr" (row index of the unrestricted "references" table).
# This allows easier comparison of the results with the original "references" table.

# "ref_ids" allows the user to pass a vector of reference IDs that will be used as IDs for newly added sources
# (references that cannot be assigned to known source in the "sources" table). 
# "ref_ids" should hold an ID for every row of the "references" table. 
# If a "ref_set" is given, the appropriate IDs will be selected from "ref_ids" automatically.

# "n.candidates" lets the user specify how many candidates should at most be passed to the aproximate string matching
# function "findBestCandidates".

# "dist_threshold" lets the user specify the maximum weighted Dammerau-Levenshtein-distance threshold 
# for a reference to be assigned to a source in the aproximate string matching part.

identifyReferences <- function(references, sources, ordered.subsets, dist_threshold=3, n.candidates=20,
                               ref_set=NULL, ref_ids=NULL, add_as_source=NULL,
                               missing_field_dists=NULL, field_weights=NULL,
                               min_separation_threshold=2, 
                               f.page_allowed_numeric_distance=10, 
                               f.page_extra_distance=3, self_run=FALSE, verbose=FALSE){
  
  # some constants for quicker computation in the main loop
  os_sums <- apply(ordered.subsets,2,sum)
  minfields <- min(os_sums)
  # get ordered field names
  fields <- rownames(ordered.subsets) 
  nfields <- nrow(ordered.subsets)
  
  if(is.null(ref_set)){
    ref_set <- 1:nrow(references)
  }else if(!all(ref_set %in% 1:nrow(references))){
    cat("\\nERROR: ref_set indices out of bounds. Check if all indices in ref_set are in the range 1:nrow(references).\\n")
    return(1)
  }
  if(is.null(ref_ids)){
    ref_ids <- 1:nrow(references)
  }
  if(is.null(add_as_source)){
    add_as_source <- rep(TRUE, nrow(references))
  }
  if(is.null(missing_field_dists)){
    missing_field_dists <- rep(1,length(fields))
    names(missing_field_dists) <- fields
  }else if(! all(fields %in% names(missing_field_dists))){
    cat("\\nERROR: not all field names present in names(missing_field_dists)\\nfield names: ", 
        fields,"\\nnames(missing_field_dists):",names(missing_field_dists),"\\n")
    return(1)
  }
  if(is.null(field_weights)){
    field_weights <- rep(1,length(fields))
    names(field_weights) <- fields
  }else if(! all(fields %in% names(field_weights))){
    cat("\\nERROR: not all field names present in names(field_weights)\\nfield names: ", 
        fields,"\\nnames(field_weights):",names(field_weights),"\\n")
    return(1)
  }
  # restrict the "references" data frame to the selection passed in "ref_set"
  references <- references[ref_set,]
  # restrict the "ref_ids" vector to the selection passed in "ref_set"
  ref_ids <- ref_ids[ref_set]
  # restrict the "add_as_source" vector to the selection passed in "ref_set"
  add_as_source <- add_as_source[ref_set]
  
  # TEMP (change later): A vector to store the max nr. of identical fields for each reference
  maxfields_vector <- numeric(nrow(references))
  
  
  print("Running reference identification:")
  if(!verbose){print(paste("Reference set size:", length(ref_set), "comparing reference:"))}
  
  # For each row of references (= for each reference). Main loop of the function.
  for(ref.nr in 1:nrow(references)){  
    
    # number of sources at this time
    ns <- nrow(sources)
    
    # initialize the data frame for source candidates 
    # (select "nothing" from sources so that candidates has the correct columns)
    candidates <- data.frame(cbind(ref.nr=numeric(0), n.identical.fields=numeric(0), sources[rep(FALSE, ns), ]), stringsAsFactors=FALSE)

    # build source == reference matrix "s.matrix", see function description above for details
    self_souid <- references$sou.id[ref.nr]
    s.matrix <- buildSourceRefMatrix(sources, references[ref.nr, fields])
    if(self_run){
      s.matrix <- s.matrix[!(sources$sou.id == self_souid),]
    }
    
    # make sure that s.matrix is a logical (0/1 or FALSE/TRUE) matrix
    s.matrix[is.na(s.matrix) | is.nan(s.matrix)] <- FALSE
    
    # number of found source candidates for this reference
    found.candidates <- 0
    
    # preselections for faster computation: 
    
    s.matrix_sums <- rowSums(s.matrix[,2:(nfields+1)])
    maxfields <- max(s.matrix_sums)
    # we can ignore sources with less corresponding values than the smallest field subset in "ordered.subsets" 
    s.matrix <- s.matrix[s.matrix_sums >= minfields,]
    
    # If only one source remains, R converts s.matrix it to a vector. It must then be transposed to work as a matrix again.
    if(!is.matrix(s.matrix)){ s.matrix <- t(s.matrix) }
    # we can ignore field subsets with reference fields that don't match any source
    missing_fields <- matrix( !apply(matrix(s.matrix[,2:(nfields+1)],ncol=6), 2, any), nrow=1 )
    # if there are no sources with enough equal field values for this reference, give up
    if(! (maxfields < minfields | all(missing_fields)) ){
      
      # 
      relevant.subset.indices <- which( (missing_fields %*% ordered.subsets == 0)[1,] )
      # resulting number of field subsets
      nss <- length(relevant.subset.indices)
      
      # end preselections
      
      # main loop for selecting source candidates:
      
      # index for subsets
      s <- 1
      # if s.matrix has only one row, it is converted to a vector by R and must be transposed to work as a matrix in R matrix multiplication
      if(!is.matrix(s.matrix)){ s.matrix <- t(s.matrix) }
      used.subsets <- c()
      # while not all subsets were checked AND s.matrix is not empty AND fewer than "n.candidates" source candidates have been found 
      while(s <= nss & nrow(s.matrix) >= 1 & found.candidates < n.candidates){
        # get the next relevant field subset (sbs) as a 1 column matrix
        sbs <- matrix(ordered.subsets[, relevant.subset.indices[s]], ncol=1)
        # Find rows of s.matrix which have TRUE entries for each field in the subset sbs
        # by matrix multiplication of s.matrix with sbs.
        # These rows correspond to sources which are equal to the reference in all fields for this field subset
        sum.sbs <- sum(sbs) 
        is_candidate <- as.vector((s.matrix[,2:(nfields+1)] %*% sbs) == sum.sbs)
        # number of source candidates for this subset
        sbs.n.cand <- sum(is_candidate)
        # if there are any candidates
        if(sbs.n.cand >= 1){
          # add up to n.candidates new source candidates to the return data frame "candidates"
          # the right source is found by the first column of s.matrix which has the corresponding row number in "sources"
          new_source_candidates <- data.frame(cbind(ref.nr=ref_set[ref.nr], n.identical.fields=sum(sbs), sources[ s.matrix[is_candidate ,1] ,] )[ 1:min(sbs.n.cand, n.candidates - found.candidates) , ], stringsAsFactors = FALSE)
          candidates <- data.frame(rbind(candidates, new_source_candidates), stringsAsFactors = FALSE) 
          # exclude these candidates from s.matrix to save computation time for further field subsets
          s.matrix <- s.matrix[ !is_candidate , ]
          # add the number of found candidates for this subset to the total number of found candidates
          found.candidates <- found.candidates + sbs.n.cand
        }
        # if we found a candidate with that matches the reference exactly in all fields, stop looking for more candidates
        if(sum(sbs) == length(fields) & sbs.n.cand >= 1){break}
        # increase field subset index
        s <- s+1
        
        # if s.matrix has only one row left, it is converted to a vector by R and must be transposed to work as a matrix again.
        if(!is.matrix(s.matrix)){
          s.matrix <- t(s.matrix)
        }else if(nrow(s.matrix) == 0){ break } # if it has zero rows left, just end the while loop
        
        
      }# end while
      
    }# end if
    
    # end find candidates loop
    
    # Call findBestCandidate for the source candidates and assign it to the reference
    # OR add the reference as a new source. 
    
    # create a new source identifier in case we need it
    new_sou.id <- ref_ids[ref.nr]
    
    # if no candidates were found, directly add the reference as a new source (if allowed).
    if(nrow(candidates) == 0){
      if(add_as_source[ref.nr]){
        new_source <- character(ncol(sources))
        names(new_source) <- names(sources)
        new_source[fields] <- references[ref.nr, fields]
        new_source["sou.id"] <- new_sou.id
        sources[nrow(sources) + 1,] <- new_source
        references[ref.nr, "sou.id"] <- paste("noCand", new_sou.id)
      }else{
        # if the reference is not allowed as a new source, mark with "noCand".
        references[ref.nr, "sou.id"] <- "noCand"
      }
      
      # otherwise call findBestCandidate to find the best source candidate 
    }else{
      
      # findBestCandidate returns either the sou.id of the best candidate or "overThresh" if all candidates are over the
      # distance threshold and "amb sou.id 1 sou.id 2... sou.id n" if the best candidates are too close in distance. 
      best_candidate <- findBestCandidate(candidates, references[ref.nr, fields], dist_threshold, 
                                          missing_field_dists=missing_field_dists,
                                          field_weights=field_weights,
                                          min_separation_threshold=min_separation_threshold,
                                          f.page_allowed_numeric_distance=f.page_allowed_numeric_distance,
                                          f.page_extra_distance=f.page_extra_distance)
      
      # if the candidate assignment is ambiguous, mark the reference with "amb " + all ambiguous sou.ids
      if(length(grep("^amb", best_candidate)) > 0){
        references[ref.nr, "sou.id"] <- best_candidate
        
        # if no valid candidate was found, add the reference as a new source and mark it with "overThresh " + its new source identifier in the references table
      }else if(length(grep("^overThresh", best_candidate)) > 0){
        if(add_as_source[ref.nr]){
          new_source <- character(ncol(sources))
          names(new_source) <- names(sources)
          new_source[fields] <- references[ref.nr, fields]
          new_source["sou.id"] <- new_sou.id
          sources[nrow(sources) + 1,] <- new_source
          references[ref.nr, "sou.id"] <- paste(c(best_candidate, " new:",new_sou.id), collapse="")
        }else{
          references[ref.nr, "sou.id"] <- paste(best_candidate, collapse="")
        }
        # if a valid best candidate was found, mark the reference with its sou.id 
      }else{
        references[ref.nr, "sou.id"] <- best_candidate
      } 
      
    }
    
    if(verbose){
      print(paste("ref:" , ref.nr, ", found:", found.candidates, ", max id fields:", maxfields))
    }else if(ref.nr %% 20 == 1){
      cat(" ",ref.nr)
    }
    maxfields_vector[ref.nr] <- maxfields
    
  }# end reference for-loop
  
  # return a list with 2 elements: the updated sources- and references- data frames
  return(list(sources, references, maxfields_vector))
  
}











########## Run Identification

### Load data

# Usually the algorithm requires a "sources" data frame, a "references" data frame
# and a vector of "reference ids" (generated in script "make.ref.ids.R").
# Reference ids can be omitted, but they are useful for mapping assigned references back to the original reference data frames,
# if some references were filtered out before the identification (e.g. only scholarly references).

#load(file=paste(data_path, "sources.RObj", sep=""))
load(file=paste(data_path, "sources.preprocessed.RObj", sep=""))
#load(file=paste(data_path, "sources.after.rest.RObj", sep=""))
#load(file=paste(data_path, "references-sel.1.RObj", sep=""))
#load(file=paste(data_path, "references-sel.2.RObj", sep=""))
#load(file=paste(data_path, "references-rest.RObj", sep=""))
#load(file=paste(data_path, "ref-rest-edited.preprocessed.RObj", sep=""))
load(file=paste(data_path, "ref-sel.1-edited.preprocessed.RObj", sep=""))
#load(file=paste(data_path, "ref-sel.2-edited.preprocessed.RObj", sep=""))
#load(file=paste(data_path, "ref.ids-sel.1.RObj", sep=""))
#load(file=paste(data_path, "ref.ids-sel.2.RObj", sep=""))
#load(file=paste(data_path, "ref.ids-rest.RObj", sep=""))
#load(file=paste(data_path, "ref.ids-rest.preprocessed.RObj", sep=""))
load(file=paste(data_path, "ref.ids-sel.1.preprocessed.RObj", sep=""))

### Data preprocessing

# convert author and publication name fields to upper case (not necessary for preprocessed)
sou$author.1 <- toupper(sou$author.1)
sou$aut.1 <- toupper(sou$aut.1)
sou$pub.name <- toupper(sou$pub.name)
sou$full.pub <- toupper(sou$full.pub)
sou$vol <- toupper(sou$vol)
sou$f.page <- toupper(sou$f.page)
sou$DOI <- toupper(sou$DOI)
ref$aut.1 <- toupper(ref$aut.1)
ref$pub.name <- toupper(ref$pub.name)
ref$vol <- toupper(ref$vol)
ref$f.page <- toupper(ref$f.page)
ref$DOI <- toupper(ref$DOI)

# STRANGE TABLE IN Y6
ref$y6 <- as.numeric(ref$y6)

# function for converting data frame columns in "factor" format to character columns
factorToChar <- function(v){
  if(is.factor(v)){
    return(as.character(v))
  }else{return(v)}
}

# Convert factor- to character-columns.
for(i in 1:ncol(sou)){
  sou[,i] <- factorToChar(sou[,i])
}

for(i in 1:ncol(ref)){
  ref[,i] <- factorToChar(ref[,i])
}

### Set algorithm inputs for the identification process

# Build a matrix of all field subsets. Change this if more or fewer fields should be compared!
all.subsets <- t(as.matrix(expand.grid(0:1, 0:1, 0:1, 0:1, 0:1, 0:1)))

# IMPORTANT: define the field names: The function "identifyReferences" will take field names from this definition!
# Both the references- and the sources-data frame must have all of these field names as column names.
# Columns with names not defined here will be ignored in the identification.
# Columns with names other than "DOI","aut.1","year","pub.name","vol","f.page" will be compared with a "standard" 
# procedure (using the unmodified Dammerau-Levenshtein distance), if they are specified here.
rownames(all.subsets) <- c("DOI","aut.1","year","pub.name","vol","f.page")

# Define field weights for the source candidate selection. Source candidates (usually up to 20 per reference) will be selected for
# a more thorough comparison to the reference, according to the sum of weights of the fields, which are identical
# to the reference.  
# (these are not the weights for distance values between reference and source candidate fields, 
# see "field_weights" below)
w.custom <- c(6,3,3,2,3,4) # weights of "importance" for the fields "DOI","aut.1","year","pub.name","vol","f.page".
# order the field-subset matrix by the sum of weights of the fields they include. 
ordered.subsets <- all.subsets[, order(w.custom %*% all.subsets, decreasing=TRUE)]
# Define a selection of field subsets whose combined field weights are larger than a constant cutoff value
# (higher cutoffs will decrease the number of field combinations tested for exact matching, but speed up the algorithm)
ssel <- w.custom %*% ordered.subsets >= 5
# Show the final, ordered, field-subset matrix used for source candidate selection.
ordered.subsets[,ssel]


# Which references are allowed to be added as a new source?
n <- nrow(ref) # total number of refrences in this dataset
tc <- apply(ref[, paste("y", 6:15, sep="")], 1, sum) # citation count over all relevant years (2006-2015) 
add_as_source <- tc >= 2 # only allow references with a minimum of 2 citations to be added as a new source
#tc_yearly_max <- apply(ref[, paste("y", 6:15, sep="")], 1, max) # max. yearly citation count
#add_as_source <- tc_yearly_max >= 2 # only allow references with a minimum of 2 citations in a single year to be added as a new source
#add_as_source <- rep(TRUE, n)
#l_60 <- as.numeric(ref$year) < 1960 | is.na(ref$year) # which references are published before 1960
#add_as_source[l_60] <- add_as_source[l_60] & tc[l_60] > 10 # only allow references published before 1960 as a new source, if the were cited at least 10 times in total

# In what order should references be processed?
ref_ordered_indices <- order(tc, decreasing=TRUE) # order reference indices by citation count
#ref_ordered_indices <- order(tc)
#ref_ordered_indices <- n:1

# assign a set size for each call of function "identifyReferences" 
# (to save intermediate results when identifying a large number of references)
nb <- 1000
ref_blocks <- c(seq(1, n, nb), n+1) # starting indices for each block

# The input variable "ref_set" should contain a set of indices for rows of the "references"-table.
# Only these references (rows) will be identified, using the order defined by this set.
ref_set <- ref_ordered_indices[(ref_blocks[1]):(ref_blocks[2]-1)]

# If you want to identify the whole "ref" table in one go or only one specific set of rows, 
# define "ref_set" accordingly and dont use the for loop below! Examples:
#ref_set <- ref_ordered_indices[1:100]
# ref_set <- 1:nrow(ref) # whole table
# ref_set <- which(ref$year %in% seq(2009:2011)) # specific selection
# ref_set <- ref_ordered_indices # whole table ordered by citation count

# Assign reference ids (if not assigned, the numbers 1 to length(ref_set) will be used by default)
ref_ids <- ref.ids_sel.1 # ref.ids_sel.1 are the loaded reference ids which were generated in the "make.ref.ids.R" script
#ref_ids <- ref.ids_sel.2
#ref_ids <- ref.ids_rest
#ref_ids <- NULL

# Define distance weights for compared fields and distance "punishments" for missing fields
field_weights <- c(20, 1, 3, 1, 2, 1)
names(field_weights) <- rownames(ordered.subsets)
missing_field_dists <- c(1, 2, 1, 2, 1, 1)
names(missing_field_dists) <- rownames(ordered.subsets)

# Set the maximum allowed numeric difference of page numbers
f.page_allowed_numeric_distance=10 
# what extra distance should be added if the page number of reference and source candidate differ by more? 
f.page_extra_distance=3

# maximum approximate string matching distance threshold to assign a reference to a source.
threshold <- 2
# maximum number of source candidates that are checked with approximate string matching after
# the initial exact string matching step.
#n.candidates <- 20
n.candidates <- 200


# run "identifyReferences" on the first block
# The result is a list with 3 elements:
#   1) The updated sources table (including newly added sources)
#   2) The table of identified references (with source id in column "sou.id", 
#       only includes references defined in "ref_set", 
#       sou.id is debug version!! See function below "analysis" section to convert debug sou.id to final sou.id). 
#   3) A vector with the highest number of identical fields (author, year, etc.) that each identified reference
#       has with any source candidate (for debug/analysis).
L <- identifyReferences(ref, sou, ordered.subsets[, ssel], ref_set=ref_set, ref_ids=ref_ids, 
                        n.candidates=n.candidates, dist_threshold=threshold, add_as_source=add_as_source,
                        missing_field_dists=missing_field_dists, field_weights=field_weights,
                        f.page_allowed_numeric_distance=f.page_allowed_numeric_distance, 
                        f.page_extra_distance=f.page_extra_distance, self_run=TRUE, verbose=FALSE)
# save intermediate result list (comment out when testing to avoid overwriting previous results!)
save(L, file=paste(data_path, "ref-sel.1_identified_List.RObj", sep="")) # save intermediate result

# STOP here when not identifying in multiple blocks with saved intermediate results!

# If identifying in blocks, continue...
# indices of the blocks in "ref_blocks" you want to identify (change to continue paused/aborted runs)
initial_i <- 2
end_i <- (length(ref_blocks)-1)

# continue identification for specified reference blocks, saving intermediate results after each block
for(i in initial_i:end_i){
  cat("\\nIdentification of reference set", i, "of", end_i, "\\n")
  sou <- L[[1]]
  ref_set <- ref_ordered_indices[(ref_blocks[i]):(ref_blocks[i+1]-1)]
  L_tmp <- identifyReferences(ref, sou, ordered.subsets[, ssel], ref_set=ref_set, ref_ids=ref_ids, 
                              n.candidates=n.candidates, dist_threshold=threshold, add_as_source=add_as_source,
                              missing_field_dists=missing_field_dists, field_weights=field_weights,
                              f.page_allowed_numeric_distance=f.page_allowed_numeric_distance, 
                              f.page_extra_distance=f.page_extra_distance, self_run=TRUE, verbose=FALSE)
  L[[1]] <- L_tmp[[1]]  # replace the updated source table with the newest version
  L[[2]] <- rbind(L[[2]], L_tmp[[2]]) # add the identified references of this block to the table of identified references
  L[[3]] <- c(L[[3]], L_tmp[[3]]) # add the maximum identical field numbers for the references of this block
  rm(L_tmp)
  save(L, file=paste(data_path, "ref-sel.1_identified_List.RObj", sep="")) # save intermediate result
}


########## Post processing


### Extract final source ids from L[[2]] and set in references 
finalSourceID <- function(ref){
  simpleSID <- gsub("^d[[:digit:][:punct:]]+ (\\\\d+)", "\\\\1", ref$sou.id)
  simpleSID <- gsub("^overThresh d[[:digit:][:punct:]]+ \\\\d+ new:(\\\\d+)", "\\\\1", simpleSID)
  simpleSID <- gsub("^overThresh d[[:digit:][:punct:]]+ \\\\d+", "-1", simpleSID)
  simpleSID <- gsub("^noCand (\\\\d+)", "\\\\1", simpleSID)
  simpleSID <- gsub("^noCand", "-1", simpleSID)
  simpleSID <- gsub("^amb d[[:digit:][:punct:]]+ (\\\\d+) [[:digit:][:space:]]+", "\\\\1", simpleSID)
  return(as.numeric(simpleSID))
}

finalDistance <- function(ref){
  dist <- gsub("^d([[:digit:][:punct:]]+) \\\\d+", "\\\\1", ref$sou.id)
  dist <- gsub("^overThresh d([[:digit:][:punct:]]+) .*", "\\\\1", dist)
  dist <- gsub("^amb d([[:digit:][:punct:]]+) [[:digit:][:space:]]+", "\\\\1", dist)
  dist[grep("noCand", dist)] <- -1
  dist <- as.numeric(dist)
  return(dist)
}

# IMPORTANT: Do not overwrite sou.ids, if you want to test identified results, 
# since meta-data is written in sou.id and will be lost if overwritten by the final "cleaned" sou.id
references$sou.id[ref_ordered_indices] <- finalSourceID(L[[2]])
#View(references)



################## Testing

# Load identification results
load(paste(data_path, "ref-sel.1_identified_List.RObj", sep=""))

# Extract data frames from results
sou2 <- L[[1]]
ref2 <- L[[2]]
max_identical_fields <- L[[3]]
simpleSID <- finalSourceID(ref2)
ref2 <- cbind(ref2, simpleSID, max_identical_fields)
View(ref2[1:10,])


# Function for manually comparing reference-source associations in a given distance interval
# Shows best associations in the interval, regardless of how the distance threshold was actually set.
# Returns: 
# 1) a table of all distances between each reference and its best source candidate.
# 2) a list of all references within the specified distance interval to their best source candidate(s)
#    and the according source candidate for comparison.
# You can set a maximum number of shown results (since this procedure is slow)
# You can set the function to return a random sample of results instead of starting from the top of "ref".
# You can show only ambiguous results (with multiple different source candidates)
compareSR <- function(ref, sou, distInterval=c(0,2), maxResults=NULL, randomSample=FALSE, onlyAmb=FALSE, origSou=FALSE){
  
  I <- seq(distInterval[1], distInterval[2], 1)
  dist <- gsub("^d([[:digit:][:punct:]]+) \\\\d+", "\\\\1", ref$sou.id)
  dist <- gsub("^overThresh d([[:digit:][:punct:]]+) .*", "\\\\1", dist)
  dist <- gsub("^amb d([[:digit:][:punct:]]+) [[:digit:][:space:]]+", "\\\\1", dist)
  dist[grep("noCand", dist)] <- -1
  dist <- as.numeric(dist)
  simpleSID <- gsub("^d[[:digit:][:punct:]]+ (\\\\d+)", "\\\\1", ref$sou.id)
  simpleSID <- gsub("^overThresh d[[:digit:][:punct:]]+ (\\\\d+).*", "\\\\1", simpleSID)
  simpleSID <- gsub("^amb d[[:digit:][:punct:]]+ ([[:digit:][:space:]]+)", "\\\\1", simpleSID)
  simpleSID <- gsub("noCand", "-1", simpleSID)
  if(onlyAmb){
    amb_indices <- grep("amb", ref$sou.id)
    u <- which(dist %in% I)[order(dist[dist %in% I])]
    u <- intersect(u, amb_indices)
  }else{
    u <- which(dist %in% I)[order(dist[dist %in% I])]
  }
  if(randomSample){u <- sample(u)}
  uRefs <- ref[u,]
  if(is.null(maxResults)){maxResults <- nrow(uRefs)}
  resultList <- list()
  
  # Table of distances to best source candidate. References without a source candidate have distance -1.
  roundDist <- dist
  resultList[[1]] <- table(roundDist)
  
  if(nrow(uRefs) == 0){
    resultList[[2]] <- paste("No references within the specified distance interval:",distInterval[1],"-",distInterval[2])
  }else{
    relevant_sou_fields <- c("sou.id","author.1","aut.1","year","pub.name","full.pub","title",
                             "vol","iss","f.page","l.page","DOI")
    if(origSou){relevant_sou_fields <- c("origSou.id", relevant_sou_fields)}
    common_fields_ref <- names(uRefs)[names(uRefs) %in% relevant_sou_fields]
    common_fields_sou <- names(sou)[names(sou) %in% relevant_sou_fields]
    for(i in 1:min(nrow(uRefs), maxResults)){
      cat(i," ")
      if(i %% 100 == 0){
        print("")
      }
      SID <- strsplit(simpleSID[u[1:maxResults]][i], " ")[[1]]
      #res_frame <- sou[NULL,relevant_sou_fields]
      res_frame <- setNames(data.frame(matrix(ncol=length(relevant_sou_fields), nrow = 0)), relevant_sou_fields)
      res_frame[1,common_fields_ref] <- uRefs[i,common_fields_ref]
      candidate_indices <- which(sou$sou.id %in% SID)
      #res_frame <- rbind(res_frame, sou[sou$sou.id %in% SID, relevant_sou_fields])
      res_frame[2:(length(candidate_indices)+1),common_fields_sou] <- sou[candidate_indices,common_fields_sou]
      resultList[[i+1]] <- res_frame
    }
  }
  
  return(resultList)
}

# Examples

# IMPORTANT: not selecting a maximum number of results will lead to a very long processing time!
positives <- compareSR(ref2, sou2, distInterval=c(0,2)) 
positives

# This will show a randomn selection of 30 ambiguous reference-source assignments
positives_selection <- compareSR(ref2, sou2, distInterval=c(0,2), maxResults = 30, randomSample=TRUE, onlyAmb=TRUE)
positives_selection

# For shown results, the first row always contains the reference and following rows contain the source candidate(s)



