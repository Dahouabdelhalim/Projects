#### Modified null modelling function, April 2021 ####

generate_null_net <- function(consumers, resources, sims = 100,
                              data.type = "names", maintain.d = NULL,
                              summary.type = "sum", c.samples = NULL,
                              r.samples = NULL, r.weights = NULL,
                              prog.count = TRUE){
  
  # Ensure input data are in data frame format (important when data are stored
  #   in tidyverse 'tibble' format)
  consumers <- as.data.frame(consumers)
  consumers[, 1] <- as.factor(consumers[, 1])
  resources <- as.data.frame(resources)
  if(!is.null(c.samples)) c.samples <- as.data.frame(c.samples)[, 1]
  if(!is.null(r.samples)) r.samples <- as.data.frame(r.samples)[, 1]
  if(!is.null(r.weights)) r.weights <- as.data.frame(r.weights)
  
  # --------------------------------------
  # Initial error handling:
  # --------------------------------------
  # 1. Column names are identical (names and order):
  if (!identical(colnames(consumers)[-1], colnames(resources))) {
    stop("Resource names do not match in 'consumers' and 'resources'")}
  if(!is.null(r.weights) && (!identical(colnames(r.weights)[-1],
                                        colnames(resources)))) {
    stop("Names of 'r.weights' do not match names of 'resources'")}
  
  # 2. Either both or neither of r.samples and c.samples are present
  if (sum(is.null(r.samples), is.null(c.samples)) == 1) {
    stop("Only one of 'r.samples' and 'c.samples' supplied")}
  
  # 3. Check resource abundance data: length and sample codes
  # When sample codes are not present, 'resources' should be a single row
  if (is.null(r.samples) && (nrow(resources) != 1)) {
    stop("Resource abundances should have one row")}
  
  # When sample codes are present, the the length of 'c.samples' and 'r.samples'
  #   should match the number of rows in 'consumers' and 'resources' respectively,
  #   and the factor levels should be the same
  if(!is.null(r.samples) && {
    (nrow(resources) != length(r.samples)) ||
      (nrow(consumers) != length(c.samples)) ||
      (!identical(sort(unique(r.samples)), sort(unique(c.samples))))}) {
    stop("There is a problem with the sample codes:
'c.samples' and/or 'r.samples'
       are of incorrect length or their sample codes do not match")}
  
  # 4. Abundance weights (typically forbidden link values) are in the range 0 - 1
  if(!is.null(r.weights)) {
    if(max(r.weights[, -1]) > 1 || min(r.weights[, -1]) < 0) {
      stop("Abundance weights ('r.weights') must be between 0 and 1")}}
  
  # 5. Issue a warning (cf. error) if a resource is consumed where the recorded
  #    abundance of the resource = 0 i.e. a resource that cannot be consumed in
  #    the null model. This may or may not require action by the user.
  if (!is.null(c.samples)) {
    spp.consumed <- stats::aggregate(consumers[, -1], by = list(c.samples), sum)
    spp.consumed <- spp.consumed[order(spp.consumed$Group.1, decreasing = FALSE), ]
    all.res <- cbind(r.samples, resources)
    all.res <- all.res[order(all.res$r.samples, decreasing = FALSE), ]
    if(sum(spp.consumed[, -1] > 0 & all.res[, -1] == 0) > 0) {warning(
      "One or more instances detected where a consumer interacted with a
       resource that has zero abundance in 'resources'")}
  } else {
    spp.consumed <- colSums(consumers[, -1])
    all.res <- resources
    if(sum(spp.consumed >0 & all.res == 0) > 0) {warning(
      "There is at least one instance where a resource was consumed but
      was zero in the abundance data (i.e. 'resources')")}
  }
  
  # 6. Correct types of data for data.type = "names" and "counts"
  if(data.type == "names" & {sum(consumers[, -1] == 1 | consumers[, -1] == 0) !=
      (nrow(consumers) * ncol(consumers[, -1]))}) {
    stop("Entries in the consumer data should equal 0 or 1")}
  if(data.type == "counts" & {sum(consumers[, -1] - round(consumers[, -1], 0)) > 0}) {
    stop("Entries in the consumer data should be integers")}
  # --------------------------------------
  
  # --------------------------------------
  # Set up a data frame that summarises the consumer data, containing:
  #   1. The original order of the data; 2. Consumer species; 3. Number of
  #      interactions per individual and total interactions (sum); 4. Sample
  #     codes (or "A" if samples not specified)
  if (!is.null(c.samples)){
    c.summary <- data.frame(ord = seq(1, nrow(consumers)), c.sample = c.samples,
                            species = consumers[, 1],
                            links = rowSums(consumers[, 2:ncol(consumers)] != 0),
                            total = rowSums(consumers[, 2:ncol(consumers)]))
  } else {
    c.summary <- data.frame(ord = seq(1, nrow(consumers)),
                            c.sample = rep("A", nrow(consumers)),
                            species = consumers[, 1],
                            links = rowSums(consumers[, 2:ncol(consumers)] != 0),
                            total = rowSums(consumers[, 2:ncol(consumers)]))
  }
  # --------------------------------------
  
  # --------------------------------------
  # Create data frame 'rd' by adding sample codes, where specified, to the
  #   resource abundance data. Otherwise a null sample value 'A' is added
  if (!is.null(c.samples)){
    rd <- cbind.data.frame(r.samples, resources)
    colnames(rd)[1] <- "r.sample"
  } else {
    rd <- cbind.data.frame(rep("A", nrow(resources)), resources)
    colnames(rd)[1] <- "r.sample"
  }
  # --------------------------------------
  
  # --------------------------------------
  #   1. Add in the summary data for individual consumers ('c.summary')
  #   2. Multiply the resource abundances by the r.weights matrix (if specified)
  rd2 <- merge(c.summary, rd, by.x = "c.sample", by.y = "r.sample")
  rd2 <- rd2[order(rd2$ord), ]
  if (!is.null(r.weights)) {
    abund.wghts <- merge(c.summary, r.weights, by.x = "species",
                         by.y = colnames(r.weights[1]))
    abund.wghts <- abund.wghts[order(abund.wghts$ord), ]
    
    # Error handling - check:
    #  1. The columns of consumer names match in the resource abundance 'rd2' and
    #      abundance weights 'abund.wghts' data frames
    #  2. The resource names and order match in 'rd2' and 'abund.wghts' before
    #      multiplying them together.
    if (!identical(abund.wghts$species, rd2$species))
    {stop("Consumer names in 'r.weights' and 'resources' do not match")}
    if (!identical(colnames(abund.wghts[, 6:ncol(abund.wghts)]),
                   colnames(rd2[6:ncol(rd2)])))
    {stop("Resource names in 'r.weights' and 'resources' do not match")}
    
    rd2[, 6:ncol(rd2)] <- rd2[, 6:ncol(rd2)] * abund.wghts[, 6:ncol(abund.wghts)]
  }
  
  # Error handling: check that the number of potential resources => the maximum
  #   number of interactions of individual consumers once any forbidden links
  #   are applied (within each sub-division of the data, where present)
  max.n.res <- data.frame(Subdiv = paste(rd2[, 1], rd2[, 3], sep = ""),
                          Links = rd2$links,
                          Nonzero = rowSums(rd2[, -c(1:5)] != 0))
  max.n.res <- stats::aggregate(max.n.res[, 2:3], by = list(max.n.res$Subdiv), max)
  if(sum(max.n.res$Links > max.n.res$Nonzero) > 0) {stop(
    "Some consumers interact with more resources than are available (>0 abundance)
     in 'resources', so the network cannot be modelled correctly")}
  # --------------------------------------
  
  # ======================================
  # Outer loop of the null model
  for (h in 1:sims){
    # Create a data frame to contain the simulated data
    res.df <- rd2[, 3:ncol(rd2)]
    res.df[, 4:ncol(res.df)] <- 0
    
    # -------------------------
    # 1. Qualitative (0/1 data)
    if (data.type == "names") {
      for (i in 1:nrow(res.df)) {
        res.choice <- sample(colnames(rd2[, 6:ncol(rd2)]), size = rd2$links[i],
                             prob = rd2[i, 6:ncol(rd2)], replace = FALSE)
        res.df[i, res.choice] <- 1
      }
    }
    
    # -------------------------
    # 2. Count data
    if (data.type == "counts") {
      if (is.null(maintain.d) | FALSE) { # Not maintaining the degree
        for (i in 1:nrow(res.df)) {
          res.choice <- sample(colnames(rd2[, 6:ncol(rd2)]), size = rd2$total[i],
                               prob = rd2[i, 6:ncol(rd2)], replace = TRUE)
          pc <- data.frame(table(res.choice))
          rownames(pc) <- pc$res.choice
          res.df[i, rownames(pc)] <- pc[rownames(pc), 2]
        }
      }
      if (!is.null(maintain.d) && TRUE) { # Maintaining the degree - 2 stage process
        for (i in 1:nrow(res.df)) {
          # Pt 1 - select the resource taxa
          choice <- sample(colnames(rd2[, 6:ncol(rd2)]), size = rd2$links[i],
                           prob = rd2[i, 6:ncol(rd2)], replace = FALSE)
          choice <- data.frame(Taxon = choice, pt1 = 1)
          # Pt 2 - allocate additional interations so the modelled total = observed total
          choice.2 <- sample(choice$Taxon, size = rd2$total[i] - rd2$links[i],
                             replace = TRUE)
          choice.2 <- data.frame(table(choice.2))
          choice <- merge(choice, choice.2, by.x = "Taxon", by.y = "choice.2")
          choice$count <- rowSums(choice[, -1])
          rownames(choice) <- choice$Taxon
          res.df[i, rownames(choice)] <- choice[rownames(choice), 4]
        }
      }
    }
    
    # -------------------------
    # 3. Measured quantities
    if (data.type == "quantities") {
      if (is.null(maintain.d) | FALSE) { # Not maintaining the degree
        for (i in 1:nrow(res.df)) {
          resource.props <- as.matrix(rd2[i, 6:ncol(rd2)] /
                                        sum(rd2[i, 6:ncol(rd2)]))
          res.choice <- gtools::rdirichlet(1, resource.props)
          res.choice <- res.choice * rd2[i, "total"]
          res.df[i, 4:ncol(res.df)] <- res.choice
        }
      }
      if (!is.null(maintain.d) && TRUE) { # Maintaining the degree - 2 stage process
        for (i in 1:nrow(res.df)) {
          # Pt 1 - select the resource taxa
          choice <- sample(colnames(rd2[, 6:ncol(rd2)]), size = rd2$links[i],
                           prob = rd2[i, 6:ncol(rd2)], replace = FALSE)
          choice <- data.frame(Taxon = choice, pt1 = 1)
          # Pt 2 - allocate interations among the selected taxa
          rp <- as.matrix(rd2[i, as.character(choice$Taxon)])
          colnames(rp) <- choice$Taxon
          resource.props <- as.matrix(rp / sum(rp))
          repeat{
            res.choice <- gtools::rdirichlet(1, resource.props)
            if(min(res.choice) > 0) break
          }
          colnames(res.choice) <- colnames(rp)
          res.choice <- res.choice * rd2[i, "total"]
          res.df[i, colnames(res.choice)] <- res.choice[, colnames(res.choice)]
        }
      }
    }
    # --------------------------------------
    
    # --------------------------------------
    # Prepare output from the iteration of the model (iterations indexed by
    #   the value of h)
    if (summary.type == "sum") {
      res.agg <- stats::aggregate(res.df[, 4:ncol(res.df)],
                                  by = list(res.df$species), sum)
      colnames(res.agg)[1] <- "Consumer"
      res.agg <- cbind(Iteration = as.matrix(rep(paste("It.", h, sep=""),
                                                 length(res.agg[, 1]))), res.agg)
      ifelse(h == 1, res.compiled <- res.agg, {
        res.compiled <- rbind(res.compiled, res.agg)})
    }
    
    if (summary.type == "mean") {
      res.agg <- stats::aggregate(res.df[, 4:ncol(res.df)],
                                  by = list(res.df$species), mean)
      colnames(res.agg)[1] <- "Consumer"
      res.agg <- cbind(Iteration = as.matrix(rep(paste("It.", h, sep=""),
                                                 length(res.agg[, 1]))), res.agg)
      ifelse(h == 1, {res.compiled <- res.agg}, {
        res.compiled <- rbind(res.compiled, res.agg)})
    }
    
    
    ### NEW CODE ###
    ### Will need a flag adding to the nullnetr object - not suitable for 
    ### test_interactions
    
    if (summary.type == "none") {
      res.agg <- res.df[, c(1, 4:ncol(res.df))]
      colnames(res.agg)[1] <- "Consumer"
      res.agg <- cbind(Iteration = as.matrix(rep(paste("It.", h, sep=""),
                                                 length(res.agg[, 1]))), res.agg)
      ifelse(h == 1, {res.compiled <- res.agg}, {
        res.compiled <- rbind(res.compiled, res.agg)})
    }
    
    ### NEW CODE ###
    
    
    # --------------------------------------
    
    # Progress count
    if (prog.count == TRUE){Sys.sleep(0.02); print(h); utils::flush.console()}
  }
  
  # --------------------------------------
  # Compile output objects:
  # 1. Observed consumer-resource interactions
  # 2. Null model outputs
  
  if (summary.type == "sum") {
    obs.interactions <- stats::aggregate(consumers[, 2:ncol(consumers)],
                                         list(consumers[, 1]), sum)
    colnames(obs.interactions)[1] <- "Consumer"
    null.interactions <- stats::aggregate(res.compiled[, 3:ncol(res.compiled)],
                                          list(res.compiled$Consumer), mean)
    colnames(null.interactions)[1] <- "Consumer"
    if(!identical(colnames(obs.interactions), colnames(null.interactions)) |
       !identical(obs.interactions$Consumer, null.interactions$Consumer))
    {stop("Different observed and expected taxon names")}
  }
  
  if (summary.type == "mean") {
    obs.interactions <- stats::aggregate(consumers[, 2:ncol(consumers)],
                                         list(consumers[, 1]), mean)
    colnames(obs.interactions)[1] <- "Consumer"
    null.interactions <- stats::aggregate(res.compiled[, 3:ncol(res.compiled)],
                                          list(res.compiled$Consumer), mean)
    colnames(null.interactions)[1] <- "Consumer"
    if(!identical(colnames(obs.interactions), colnames(null.interactions) )|
       !identical(obs.interactions$Consumer, null.interactions$Consumer))
    {stop("Different observed and expected names")}
  }
  
  
  ### NEW CODE ###
  ### The same as for summary.type == "sum" for now - just to avoid error
  
  
  if (summary.type == "none") {
    obs.interactions <- stats::aggregate(consumers[, 2:ncol(consumers)],
                                         list(consumers[, 1]), sum)
    colnames(obs.interactions)[1] <- "Consumer"
    null.interactions <- stats::aggregate(res.compiled[, 3:ncol(res.compiled)],
                                          list(res.compiled$Consumer), mean)
    colnames(null.interactions)[1] <- "Consumer"
    if(!identical(colnames(obs.interactions), colnames(null.interactions)) |
       !identical(obs.interactions$Consumer, null.interactions$Consumer))
    {stop("Different observed and expected taxon names")}
  }
  
  ### NEW CODE ###
  
  
  # --------------------------------------
  
  # --------------------------------------
  # Create output list
  null.net.res <- list(rand.data = res.compiled, obs.interactions = obs.interactions,
                       n.iterations = sims)
  class(null.net.res) <- "nullnet"
  return(null.net.res)
  # --------------------------------------
  
}

