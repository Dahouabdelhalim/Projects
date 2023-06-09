suppressPackageStartupMessages({
   require(R6) 
   require(tidyverse)
   require(RVenn)
   require(ggExtra)
   require(cowplot)
   #require(venneuler)
   
})

# R6 Class RBPbaseDataHandler
#' Class providing object with methods for RBPbase data handling
#' @docType class 
#' @importFrom R6 R6Class 
#' @export #' 
RBPbaseDataHandler <- R6::R6Class("RBPbaseDataHandler", 
                                  list(
                                     # list of organism table data 
                                     COMPILED_TABLE    = NULL, 
                                     # list of studies
                                     RIC_STUDIES       = NULL, 
                                     ANNO_STUDIES      = NULL, 
                                     ALLANNO           = NULL,
                                     DEBUG             = FALSE,
                                     DEBUG_LEVEL       = 0,
                                     annotation_prefix = "CANNO_",
                                     study_prefix      = "CSTUDY_",
                                     
                                     initialize = function(COMPILED_TABLE, 
                                                           RIC_STUDIES, 
                                                           ANNO_STUDIES,
                                                           DEBUG = FALSE,
                                                           DEBUG_LEVEL = 0) {
                                        stopifnot(is.list(COMPILED_TABLE))
                                        stopifnot(is_tibble(RIC_STUDIES))
                                        stopifnot(is_tibble(ANNO_STUDIES))
                                        stopifnot(is.logical(DEBUG))
                                        stopifnot(is.numeric(DEBUG_LEVEL))
                                        
                                        self$COMPILED_TABLE <- COMPILED_TABLE
                                        self$RIC_STUDIES    <- RIC_STUDIES
                                        self$ANNO_STUDIES   <- ANNO_STUDIES
                                        self$DEBUG          <- DEBUG
                                        self$DEBUG_LEVEL    <- DEBUG_LEVEL
                                        self$compute_allanno()
                                     },
                                     
                                     print = function(...) {
                                        cat("RBPbase Handler")
                                        invisible(self)
                                     },
                                     
                                     #' @param organism organism short code 
                                     #' @param selection_ids only specific IDs will be displayed
                                     #' @param descriptive_names rename ids to descriptive naems
                                     #' @param count_hits_for_studies logical, if a column with sums should be returned 
                                     #' @param numeric convert logical columns into numeric 
                                     #' @param annotation if custom annotations should be added 
                                     #' @param studies if custom studies should be added
                                     #' @param filter logic filter to filter rows
                                     #' @return tibble with compiled table for organism
                                     get_compiled_table_for_organism = function(organism, 
                                                                                selection_ids = NULL,
                                                                                descriptive_names = F,
                                                                                count_hits_for_studies = NULL,
                                                                                numeric = F,
                                                                                annotation = T,
                                                                                studies = T,
                                                                                filter = NULL) {
                                        stopifnot(is.null(selection_ids) || is.character(selection_ids))
                                        stopifnot(is.logical(descriptive_names))
                                        stopifnot(is.null(count_hits_for_studies) ||  is.character(count_hits_for_studies))
                                        stopifnot(is.logical(numeric))
                                        self$assert_organism_in_rbpbase(organism)
                                        stopifnot(is.logical(annotation))
                                        stopifnot(is.logical(studies))
                                        stopifnot(length(organism) == 1) 
                                        stopifnot(is.null(filter) || is.logical(filter))
                                        
                                        if(self$DEBUG && self$DEBUG_LEVEL >= 2) {
                                           cat(paste("get_compiled_table_for_organism: arg organism", organism, "\\n"))
                                           cat(paste("get_compiled_table_for_organism: arg selection_ids", paste(selection_ids, collapse = ", "), "\\n"))
                                           cat(paste("get_compiled_table_for_organism: arg descriptive_names", descriptive_names, "\\n"))
                                           cat(paste("get_compiled_table_for_organism: arg count_hits_for_studies", paste(as.character(count_hits_for_studies), collapse = ", "), "\\n"))
                                           cat(paste("get_compiled_table_for_organism: arg numeric", numeric, "\\n"))
                                           cat(paste("get_compiled_table_for_organism: arg annotation", annotation, "\\n"))
                                           cat(paste("get_compiled_table_for_organism: arg studies", studies, "\\n"))
                                        }
                                        
                                        x <- self$COMPILED_TABLE[[organism]]
                                        
                                        if(!is.null(selection_ids) && length(selection_ids) > 0) {
                                           self$assert_ids_in_rbpbase_for_organism(selection_ids, organism)  
                                           x <- x %>% dplyr::select(UNIQUE, any_of(selection_ids))
                                        }
                                        
                                        
                                        if(is.character(count_hits_for_studies)) {
                                           self$assert_ids_in_rbpbase(count_hits_for_studies)
                                           x <- x %>% mutate(number_hits_compare = x %>% 
                                                                dplyr::select(any_of(count_hits_for_studies)) %>% rowSums)
                                        }
                                        
                                        if(!annotation) {
                                           x <- x %>% dplyr::select(-starts_with(self$annotation_prefix))
                                        }
                                        
                                        if(!studies) {
                                           x <- x %>% dplyr::select(-starts_with(self$study_prefix))
                                        }          
                                        
                                        
                                        if(is.logical(filter)) {
                                           if(length(filter) != nrow(x)) {
                                              stop(paste0("Length of filter `", length(filter), "` needs to be length of COMPILED_TABLE `", nrow(x), "`")) 
                                           }
                                           
                                           x <- x %>% dplyr::filter(filter)
                                        }
                                        
                                        
                                        if(numeric) {
                                           x <- x %>% mutate_if(is.logical, ~(as.numeric(.)))
                                        }
                                        
                                        
                                        
                                        # change names to descriptive names
                                        if(descriptive_names) {
                                           x <- x %>% self$annotate_column_names()
                                        }
                                        
                                        x
                                     },
                                     
                                     # check if organism exists in rbpbase
                                     assert_organism_in_rbpbase = function(organism) {
                                        if(!all(organism %in% names(self$COMPILED_TABLE))) {
                                           x <- organism[!organism %in% names(self$COMPILED_TABLE)]
                                           stop(paste0("Organism(s) `", paste(x, collapse = ", ") , "` do(es) not exist in RBPbase"))  
                                        }
                                        
                                        invisible(self)
                                     },
                                     
                                     # check if id exists in rbpbase
                                     assert_ids_in_rbpbase = function(ids) {
                                        stopifnot(is.character(ids))
                                        
                                        if(!all(ids %in% self$ALLANNO$RBPBASEID)) {
                                           notin <- ids[ !ids %in% self$ALLANNO$RBPBASEID ]
                                           stop(paste0("Ids ", paste(notin, collapse = ", "), " does not exist in RBPbase"))  
                                        }
                                        
                                        invisible(self)
                                     },
                                     
                                     is_ids_in_rbpbase_for_organism = function(ids, organism) {
                                        stopifnot(is.character(ids))
                                        stopifnot(is.character(organism))
                                        
                                        x <- colnames(self$COMPILED_TABLE[[organism]])
                                        
                                        return(all(ids %in% x))
                                     },
                                     
                                     # check if id exists in rbpbase for organism  (includes any_* etc)
                                     assert_ids_in_rbpbase_for_organism = function(ids, organism) {
                                        stopifnot(is.character(ids))
                                        stopifnot(is.character(organism))
                                        
                                        x <- colnames(self$COMPILED_TABLE[[organism]])
                                        
                                        if(!self$is_ids_in_rbpbase_for_organism(ids, organism)) {
                                           notin <- ids[ !ids %in% x ]
                                           stop(paste0("Ids ", paste(notin, collapse = ", "), " does not exist in RBPbase"))  
                                        }
                                        
                                        invisible(self)
                                     },
                                     
                                     
                                     get_organisms_for_ids = function(ids) {
                                        self$assert_ids_in_rbpbase(selected_studies_ids)
                                        
                                        self$ALLANNO %>% 
                                           filter(RBPBASEID %in% selected_studies_ids) %>%
                                           pull(Organism) %>% 
                                           unique 
                                     },
                                     
                                     get_study_ids_for_organisms = function(organisms) {
                                        self$assert_organism_in_rbpbase(organisms)
                                        
                                        self$RIC_STUDIES %>% 
                                           filter(Organism %in% organisms) %>%
                                           pull(RBPBASEID) %>% 
                                           unique 
                                     },
                                     
                                     annotate_column_names = function(x) {
                                        colnames(x) <- self$get_descriptiveid_for_rbpbaseids(colnames(x))
                                        x
                                     },
                                     
                                     get_descriptiveid_for_rbpbaseids = function(x) {
                                        data.frame(base = x, stringsAsFactors = F) %>% 
                                           left_join(y = self$ALLANNO, by = c("base" = "RBPBASEID")) %>%
                                           mutate(DescriptiveID = ifelse(is.na(DescriptiveID), base, paste0(DescriptiveID))) %>%
                                           pull(DescriptiveID)
                                     },
                                     
                                     ##' @name add_annotation
                                     ##'        
                                     ##' @param name descriptive short identifier of study
                                     ##' @param values values for the COMPILED_TABLE column
                                     ##' @return RBPBASEID 
                                     add_annotation = function(name, values, organism, add_prefix = T,
                                                               description = "", curationInfo = "") {
                                        stopifnot(is.character(name))
                                        stopifnot(is.logical(values) || is.character(values) || is.factor(values)) 
                                        
                                        if(add_prefix) {
                                           id <- paste0(self$annotation_prefix, name)
                                        } else {
                                           id <- name 
                                        }
                                        
                                        self$add_study_to_compiled_table(id = id, 
                                                                         values = values,
                                                                         organism = organism)
                                        
                                        self$ANNO_STUDIES <- self$ANNO_STUDIES %>% 
                                           filter(RBPBASEID != id) %>%
                                           add_row(RBPBASEID = id,
                                                   DescriptiveID = name,
                                                   Organism = organism,
                                                   Format = class(values),
                                                   Description = description,
                                                   CurationInfo = curationInfo 
                                           )
                                                   
                                        
                                        self$compute_allanno()
                                        
                                        id
                                     },
                                     
                                     
                                     ##' @name add_study
                                     ##' 
                                     ##' @param name descriptive short identifier of study
                                     ##' @param values indicates which genes in the COMPILED_TABLE
                                     ##'                     are part of the custom study
                                     ##' @return RBPBASEID 
                                     add_study = function(name, values, organism) {
                                        stopifnot(is.character(name))
                                        stopifnot(is.logical(values) || is.character(values) || is.factor(values)) 
                                        
                                        id = paste0(self$study_prefix, name)
                                        
                                        
                                        self$add_study_to_compiled_table(id = id, 
                                                                         values = values,
                                                                         organism = organism)
                                        
                                        self$RIC_STUDIES <- self$RIC_STUDIES %>% 
                                           filter(RBPBASEID != id) %>%
                                           add_row(RBPBASEID = id,
                                                   DescriptiveID = name,
                                                   Organism = organism)
                                        self$compute_allanno()
                                        
                                        id
                                     },
                                     
                                     ##' @name add_study_to_compiled_table
                                     ##' 
                                     ##' @param id RBPBASEID
                                     ##' @param values indicates which genes in the COMPILED_TABLE
                                     ##'                     are part of the custom study
                                     ##' @param organism short organism descriptor e.g. Hs, Mm, ..
                                     ##' @return RBPBASEID 
                                     add_study_to_compiled_table = function(id, values, organism) {
                                        self$assert_organism_in_rbpbase(organism)
                                        stopifnot(is.logical(values) || is.character(values) || is.factor(values)) 
                                        
                                        if(length(values) != nrow(self$COMPILED_TABLE[[organism]])) {
                                           stop("logic vector has to have the same size as the numbers of rows in the compiled table")
                                        }
                                        
                                        self$COMPILED_TABLE[[organism]] <-  self$COMPILED_TABLE[[organism]] %>%
                                           mutate(!!id := values) 
                                        
                                        invisible(self)
                                     },
                                     
                                     remove_all_annotations = function(organism) {
                                        self$assert_organism_in_rbpbase(organism)
                                        
                                        
                                        self$COMPILED_TABLE[[organism]] <- self$COMPILED_TABLE[[organism]] %>%
                                           dplyr::select(-starts_with(self$annotation_prefix))
                                     },
                                     
                                     
                                     remove_all_studies = function(organism) {
                                        self$assert_organism_in_rbpbase(organism)
                                        
                                        self$COMPILED_TABLE[[organism]] <- self$COMPILED_TABLE[[organism]] %>%
                                           dplyr::select(-starts_with(self$study_prefix))
                                     },
                                     
                                     
                                     compute_allanno = function() {
                                        self$ALLANNO <- bind_rows(self$RIC_STUDIES, self$ANNO_STUDIES)
                                     }
                                  )
)



# R6 Class RBPbaseAnalysis
#' Class providing object with methods for analysis or RBPbase data
#' 
#' @docType class 
#' @importFrom R6 R6Class 
#' @import RVenn
#' @import ggExtra
#' @importFrom cowplot save_plot
#' @export #' 
#' 
RBPbaseAnalysis <- R6Class("RBPbaseAnalysis", 
                           list(
                              RBPbaseDataHandler          = NULL, 
                              selected_studies_ids        = NULL,
                              comparison_studies_ids      = NULL,
                              excluded_comparison_studies = NULL,
                              comparison_studies_organism = NULL,
                              selected_studies_organism   = NULL,
                              knownRBPs_id                = NULL,
                              known_RIC_desc_id           = NULL,
                              enzyme_desc_id              = NULL,
                              enzyme_id                   = NULL,
                              hasRBD_id                   = NULL,
                              metabolic_enzyme_id         = NULL,
                              knownRBPs_study_ids         = NULL,
                              anyRIC_id                   = NULL,
                              DEBUG                       = FALSE,
                              DEBUG_LEVEL                 = 0,
                              
                              initialize = function(RBPbaseDataHandler, 
                                                    selected_studies_ids, 
                                                    comparison_studies_organism,
                                                    excluded_comparison_studies = NULL,
                                                    DEBUG                       = FALSE,
                                                    DEBUG_LEVEL                 = 0,
                                                    enzyme_id                   = NULL, 
                                                    hasRBD_id                   = NULL,
                                                    metabolic_enzyme_id         = NULL,
                                                    knownRBPs_study_ids         = NULL,
                                                    anyRIC_id                   = NULL,
                                                    knownRBPs_id                = "knownRBPs",
                                                    known_RIC_desc_id           = "known_vs_RIC",
                                                    enzyme_desc_id              = "enzyme",
                                                    update                      = TRUE) {
                                 # checks
                                 stopifnot("RBPbaseDataHandler" %in% class(RBPbaseDataHandler))
                                 stopifnot(is.character(selected_studies_ids))
                                 stopifnot(is.character(comparison_studies_organism))
                                 stopifnot(is.null(excluded_comparison_studies) || is.character(excluded_comparison_studies))
                                 stopifnot(is.null(knownRBPs_id)                || is.character(knownRBPs_id))
                                 stopifnot(is.null(known_RIC_desc_id)           || is.character(known_RIC_desc_id))
                                 stopifnot(is.null(enzyme_desc_id)              || is.character(enzyme_desc_id))
                                 stopifnot(is.null(enzyme_id)                   || is.character(enzyme_id))
                                 stopifnot(is.null(knownRBPs_study_ids)         || is.character(knownRBPs_study_ids))
                                 stopifnot(is.null(hasRBD_id)                   || is.character(hasRBD_id))
                                 stopifnot(is.null(metabolic_enzyme_id)         || is.character(metabolic_enzyme_id))
                                 stopifnot(is.null(anyRIC_id)                   || is.character(anyRIC_id))
                                 stopifnot(is.logical(DEBUG))
                                 stopifnot(is.numeric(DEBUG_LEVEL))
                                 stopifnot(is.logical(update))
                                 
                                 
                                 if(DEBUG){
                                    cat("RBPbaseAnalysis initialize\\n") 
                                 }
                                 
                                 self$RBPbaseDataHandler          <- RBPbaseDataHandler
                                 self$selected_studies_ids        <- selected_studies_ids
                                 self$comparison_studies_organism <- comparison_studies_organism
                                 self$excluded_comparison_studies <- excluded_comparison_studies
                                 self$selected_studies_organism   <- RBPbaseDataHandler$get_organisms_for_ids(selected_studies_ids) 
                                 self$knownRBPs_study_ids         <- knownRBPs_study_ids
                                 self$knownRBPs_id                <- knownRBPs_id
                                 self$known_RIC_desc_id           <- known_RIC_desc_id
                                 self$enzyme_id                   <- enzyme_id
                                 self$enzyme_desc_id              <- enzyme_desc_id
                                 self$metabolic_enzyme_id         <- metabolic_enzyme_id
                                 self$DEBUG                       <- DEBUG
                                 self$DEBUG_LEVEL                 <- DEBUG_LEVEL
                                 
                                 if(!length(self$selected_studies_organism) == 1) {
                                    stop("Studies to be analysed cannot be from more than one organism.\\n")  
                                 }
                                 
                                 # set comparison studies
                                 x <- self$get_studies_ids_for_analysis() 
                                 self$comparison_studies_ids <- x[!x %in% self$selected_studies_ids]
                                 
                                 # setter 
                                 if(!is.null(anyRIC_id)) {
                                    self$set_anyRIC_id(anyRIC_id, update = FALSE)
                                 } else {
                                    # set default for any RIC id
                                    self$set_anyRIC_id(paste0("any_", self$selected_studies_organism), update = FALSE)
                                 }
                                 if(!is.null(enzyme_id)) {
                                    self$set_enzyme_id(enzyme_id, update = FALSE)
                                 }
                                 if(!is.null(hasRBD_id)) {
                                    self$set_hasRBD_id(hasRBD_id, update = FALSE)           
                                 }
                                 if(!is.null(metabolic_enzyme_id)) {
                                    self$set_metabolic_enzyme_id(metabolic_enzyme_id, update = FALSE)           
                                 }
                                 if(!is.null(knownRBPs_study_ids)) {
                                    self$set_knownRBPs_study_ids(knownRBPs_study_ids, update = FALSE)  
                                 }
                                 
                                 
                                 if(update) {
                                    if(self$DEBUG) {
                                       cat("initialize: update annotations\\n")
                                    }
                                    self$update_annotations()
                                 }
                                 
                                 invisible(self)
                              },
                              
                              print = function(...) {
                                 a <- function(...) { paste(..., collapse = ", ") }
                                 
                                 cat("RBPbase Analysis Handler\\n")
                                 cat(paste0("selected studies: ", a(self$selected_studies_ids), "\\n"))
                                 cat(paste0("comparison studies organism: ", a(self$comparison_studies_organism), "\\n"))
                                 cat(paste0("excluded comparison studies: ", a(self$excluded_comparison_studies), "\\n"))
                                 cat(paste0("selected_studies_organism: ", a(self$selected_studies_organism), "\\n"))
                                 cat(paste0("knownRBPs_id: ", a(self$knownRBPs_id ), "\\n"))
                                 cat(paste0("known_RIC_desc_id: ", a(self$known_RIC_desc_id ), "\\n"))
                                 cat(paste0("enzyme_desc_id: ", a(self$enzyme_desc_id ), "\\n"))
                                 cat(paste0("enzyme_id: ", a(self$enzyme_id ), "\\n"))
                                 cat(paste0("hasRBD_id: ", a(self$hasRBD_id ), "\\n"))
                                 cat(paste0("metabolic_enzyme_id: ", a(self$metabolic_enzyme_id ), "\\n"))
                                 cat(paste0("knownRBPs_study_ids: ", a(self$knownRBPs_study_ids ), "\\n"))
                                 cat(paste0("anyRIC_id: ", a(self$anyRIC_id ), "\\n"))
                                 cat(paste0("DEBUG: ", a(self$DEBUG ), "\\n"))
                                 
                                 invisible(self)
                              },
                              
                              
                              get_studies_ids_for_analysis = function() {
                                 self$RBPbaseDataHandler$get_study_ids_for_organisms(
                                    c(self$selected_studies_organism,
                                      self$comparison_studies_organism)) %>%
                                    self$filter_excluded_study_ids()
                              },
                              
                              filter_excluded_study_ids = function(x) {
                                 stopifnot(is.character(x))
                                 
                                 x[!x %in% self$excluded_comparison_studies]
                              },
                              
                              
                              get_main_table = function(...) {
                                 if(self$DEBUG && self$DEBUG_LEVEL >= 2) {
                                    args <- list(...)
                                    cat(paste("get_main_table: arg ...:", 
                                              paste(names(args), 
                                                    args,
                                                    collapse = ", "), 
                                              "\\n"))
                                 }
                                 
                                 self$RBPbaseDataHandler$get_compiled_table_for_organism(
                                    organism = self$selected_studies_organism,
                                    ...
                                 )
                              },
                              
                              
                              get_studies_only_ids = function(selected_studies      = T,
                                                              comparison_studies    = T,
                                                              main_organism_studies = F,
                                                              additional_ids        = NULL) {
                                 stopifnot(is.logical(selected_studies))
                                 stopifnot(is.logical(comparison_studies))
                                 stopifnot(is.logical(main_organism_studies))
                                 stopifnot(is.null(additional_ids) || is.character(additional_ids))
                                 
                                 if(!(selected_studies || comparison_studies || main_organism_studies)) 
                                    stop("one of selected_studies, comparison_studies or main_organism_studies must be TRUE.")   
                                 
                                 if(self$DEBUG && self$DEBUG_LEVEL >= 2) {
                                    cat(paste("get_studies_only_ids: arg selected_studies:", selected_studies, "\\n"))
                                    cat(paste("get_studies_only_ids: arg comparison_studies:", comparison_studies, "\\n"))
                                    cat(paste("get_studies_only_ids: arg main_organism_studies:", main_organism_studies, "\\n"))
                                    cat(paste("get_studies_only_ids: arg additional_ids:", paste(additional_ids, collapse = ", "), "\\n"))
                                 }
                                 
                                 
                                 # create the selection ids
                                 selection_ids <- c()
                                 if(selected_studies) 
                                    selection_ids <- c(selection_ids, self$selected_studies_ids)
                                 if(comparison_studies)  
                                    selection_ids <- c(selection_ids, self$comparison_studies_ids)
                                 if(main_organism_studies) {
                                    x <- self$RBPbaseDataHandler$get_study_ids_for_organisms(self$selected_studies_organism)
                                    x <- x[!x %in% self$selected_studies_ids]
                                    selection_ids <- c(selection_ids, x)
                                 }
                                 
                                 if(!is.null(additional_ids)) {
                                    selection_ids <- c(selection_ids, additional_ids) 
                                 }
                                 
                                 # exclude studies to be excluded
                                 selection_ids <- self$filter_excluded_study_ids(selection_ids)
                                 
                                 selection_ids
                              },
                              
                              
                              get_studies_only_table = function(selected_studies      = T,
                                                                comparison_studies    = T,
                                                                main_organism_studies = F,
                                                                count_comparison_hits = F,
                                                                additional_ids        = NULL,
                                                                ...) {
                                 stopifnot(is.logical(count_comparison_hits))
                                 stopifnot(is.logical(selected_studies))
                                 stopifnot(is.logical(comparison_studies))
                                 stopifnot(is.logical(main_organism_studies))
                                 stopifnot(is.null(additional_ids) || is.character(additional_ids))
                                 
                                 if(self$DEBUG && self$DEBUG_LEVEL >= 2) {
                                    cat(paste("get_studies_only_table: arg selected_studies:", selected_studies , "\\n"))
                                    cat(paste("get_studies_only_table: arg comparison_studies:", comparison_studies , "\\n"))
                                    cat(paste("get_studies_only_table: arg main_organism_studies:", main_organism_studies , "\\n"))
                                    cat(paste("get_studies_only_table: arg count_comparison_hits:", count_comparison_hits , "\\n"))
                                    cat(paste("get_studies_only_table: arg additional_ids:", additional_ids , "\\n"))
                                    args <- list(...)
                                    
                                    if(length(args) > 0) {
                                       cat(paste("get_studies_only_table: arg ...:", paste(names(args), args, collapse = ", "), "\\n"))
                                    }
                                 }
                                 
                                 selection_ids = self$get_studies_only_ids(selected_studies      = selected_studies,
                                                                           comparison_studies    = comparison_studies,
                                                                           main_organism_studies = main_organism_studies,
                                                                           additional_ids        = additional_ids)
                                 
                                 # ifelse throws error, not long version, reason unknown
                                 if(count_comparison_hits) {
                                    count_hits_for_studies <- self$comparison_studies_ids 
                                 } else {
                                    count_hits_for_studies <- NULL
                                 }
                                 
                                 
                                 if(self$DEBUG) {
                                    cat(paste("get_studies_only_table: count_hits_for_studies:", paste(count_hits_for_studies, collapse = ", "), "\\n"))
                                 }
                                 
                                 self$get_main_table(selection_ids          = selection_ids,
                                                     count_hits_for_studies = count_hits_for_studies,
                                                     ...)
                              },
                              
                              add_annotation = function(...) {
                                 self$RBPbaseDataHandler$add_annotation(..., organism = self$selected_studies_organism)
                              },
                              
                              add_study = function(...) {
                                 self$RBPbaseDataHandler$add_study(..., organism = self$selected_studies_organism)
                              },
                              
                              remove_all_annotations = function(...) {
                                 self$RBPbaseDataHandler$remove_all_annotations(..., organism = self$selected_studies_organism)
                              },
                              
                              remove_all_studies = function(...) {
                                 self$RBPbaseDataHandler$remove_all_studies(..., organism = self$selected_studies_organism)
                              },
                              
                              
                              add_knownRBPs_annotation = function(ids = NULL) {
                                 stopifnot(is.null(ids) || is.character(ids))
                                 
                                 
                                 if(is.null(ids)) {
                                    ids <- self$knownRBPs_study_ids
                                    
                                    if(is.null(ids) || length(ids) == 0) {
                                       stop(paste("add_knownRBPs_annotation can only be executed",
                                                  "when there are either ids provided or",
                                                  "`knownRBPs_study_ids` stored"))  
                                    }
                                 }
                                 
                                 if(self$DEBUG) {
                                    cat(paste("add_knownRBPs_annotation: arg ids:", paste0(ids, collapse = ", "), "\\n"))
                                 }
                                 
                                 
                                 self$add_annotation(
                                    name = self$knownRBPs_id,
                                    values = 
                                       self$get_main_table() %>% 
                                       dplyr::select(any_of(ids)) %>%
                                       mutate(property = rowSums(.) > 0) %>%
                                       pull(property),
                                    add_prefix = F)
                              },
                              
                              add_enzyme_annotation = function(enzyme_id = NULL, metabolic_enzyme_id = NULL) {
                                 stopifnot(is.null(enzyme_id) || is.character(enzyme_id))
                                 stopifnot(is.null(enzyme_id) || is.character(metabolic_enzyme_id))
                                 
                                 if(is.null(enzyme_id)) {
                                    enzyme_id <- self$enzyme_id 
                                 } 
                                 
                                 if(is.null(metabolic_enzyme_id)) {
                                    metabolic_enzyme_id <- self$metabolic_enzyme_id 
                                 }
                                 
                                 self$RBPbaseDataHandler$assert_ids_in_rbpbase_for_organism(enzyme_id, 
                                                                                            self$selected_studies_organism)
                                 self$RBPbaseDataHandler$assert_ids_in_rbpbase_for_organism(metabolic_enzyme_id, 
                                                                                            self$selected_studies_organism)
                                 
                                 
                                 self$add_annotation(
                                    name = self$enzyme_desc_id,
                                    values = 
                                       self$get_main_table() %>% 
                                       mutate(property = self$transform_enzyme(Enzyme      = !!sym(enzyme_id), 
                                                                               Metabolic.Enzyme = !!sym(metabolic_enzyme_id))) %>%
                                       pull(property),
                                    add_prefix = F)
                              },
                              
                              add_known_RIC_desc_annotation = function(knownRBPs_id = NULL, anyRIC_id = NULL) {
                                 stopifnot(is.null(knownRBPs_id) || is.character(knownRBPs_id))
                                 stopifnot(is.null(anyRIC_id)    || is.character(anyRIC_id))
                                 
                                 if(is.null(knownRBPs_id)) {
                                    knownRBPs_id <- self$knownRBPs_id
                                 }
                                 
                                 if(is.null(anyRIC_id)) {
                                    anyRIC_id <- self$anyRIC_id
                                 }
                                 
                                 if(self$DEBUG) {
                                    cat(paste0("add_known_RIC_desc_annotation: knownRBPs_id: ", knownRBPs_id, "\\n"))
                                    cat(paste0("add_known_RIC_desc_annotation: anyRIC_id: ", anyRIC_id, "\\n"))
                                 }
                                 
                                 self$RBPbaseDataHandler$assert_ids_in_rbpbase_for_organism(knownRBPs_id, 
                                                                                            self$selected_studies_organism)
                                 self$RBPbaseDataHandler$assert_ids_in_rbpbase_for_organism(anyRIC_id, 
                                                                                            self$selected_studies_organism)
                                 
                                 self$add_annotation(
                                    name = self$known_RIC_desc_id,
                                    values = 
                                       self$get_main_table() %>% 
                                       mutate(property = self$transform_known(eval(!!sym(knownRBPs_id)), 
                                                                              eval(!!sym(anyRIC_id)))) %>%
                                       pull(property),
                                    add_prefix = F 
                                 )
                              },
                              
                              
                              has_knownRBPs_study_ids = function() {
                                 !is.null(self$knownRBPs_study_ids) &&
                                    self$RBPbaseDataHandler$is_ids_in_rbpbase_for_organism(self$knownRBPs_study_ids, 
                                                                                           self$selected_studies_organism)
                              },
                              
                              has_enzyme_desc_annotation = function() {
                                 !is.null(self$enzyme_desc_id) &&
                                    self$RBPbaseDataHandler$is_ids_in_rbpbase_for_organism(self$enzyme_desc_id, 
                                                                                           self$selected_studies_organism)
                              },
                              
                              has_enzyme_annotation = function() {
                                 !is.null(self$enzyme_id) &&
                                    self$RBPbaseDataHandler$is_ids_in_rbpbase_for_organism(self$enzyme_id, 
                                                                                           self$selected_studies_organism)
                              },
                              
                              has_metabolic_enzyme_annotation = function() {
                                 !is.null(self$metabolic_enzyme_id) &&
                                    self$RBPbaseDataHandler$is_ids_in_rbpbase_for_organism(self$metabolic_enzyme_id, 
                                                                                           self$selected_studies_organism)
                              },
                              
                              has_hasRBD_annotation = function() {
                                 x <- (!is.null(self$hasRBD_id))
                                 y <- self$RBPbaseDataHandler$is_ids_in_rbpbase_for_organism(self$hasRBD_id, 
                                                                                             self$selected_studies_organism)
                                 return(x && y)
                              },
                              
                              has_anyRIC_annotation = function() {
                                 !is.null(self$anyRIC_id) &&
                                    self$RBPbaseDataHandler$is_ids_in_rbpbase_for_organism(self$anyRIC_id, 
                                                                                           self$selected_studies_organism)
                              },
                              
                              
                              get_knownRBPs_id = function() {
                                 self$knownRBPs_id
                              },
                              
                              get_enzyme_desc_id = function() {
                                 self$enzyme_desc_id 
                              },
                              
                              get_enzyme_id = function() {
                                 self$enzyme_id
                              },
                              
                              get_hasRBD_id = function() {
                                 self$hasRBD_id
                              },
                              
                              get_knownRBPs_study_ids = function() {
                                 self$knownRBPs_study_ids
                              },
                              
                              get_anyRIC_id = function() {
                                 self$anyRIC_id
                              },

                              
                              get_known_RIC_desc_id = function() {
                                 self$known_RIC_desc_id
                              },
                              
                              set_anyRIC_id = function(id, update = TRUE) {
                                 stopifnot(is.character(id) && length(id) == 1)
                                 self$RBPbaseDataHandler$assert_ids_in_rbpbase_for_organism(id, 
                                                                                            self$selected_studies_organism) 
                                 self$anyRIC_id = id
                                 
                                 if(update){
                                    self$update_annotations() 
                                 }
                              },
                              
                              set_enzyme_id = function(id, update = TRUE) {
                                 stopifnot(is.character(id) && length(id) == 1)
                                 self$RBPbaseDataHandler$assert_ids_in_rbpbase_for_organism(id, 
                                                                                            self$selected_studies_organism)
                                 
                                 self$enzyme_id = id
                                 
                                 if(update){
                                    self$update_annotations() 
                                 }
                              },
                              
                              set_hasRBD_id = function(hasRBD_id, update = TRUE) {
                                 stopifnot(is.character(hasRBD_id) && length(hasRBD_id) == 1)
                                 self$RBPbaseDataHandler$assert_ids_in_rbpbase_for_organism(hasRBD_id, 
                                                                                            self$selected_studies_organism)
                                 
                                 self$hasRBD_id = hasRBD_id
                                 
                                 if(update){
                                    self$update_annotations() 
                                 }
                              },
                              
                              set_metabolic_enzyme_id = function(metabolic_enzyme_id, update = TRUE) {
                                 stopifnot(is.character(metabolic_enzyme_id) && length(metabolic_enzyme_id) == 1)
                                 self$RBPbaseDataHandler$assert_ids_in_rbpbase_for_organism(metabolic_enzyme_id, 
                                                                                            self$selected_studies_organism)
                                 
                                 self$metabolic_enzyme_id = metabolic_enzyme_id
                                 if(update){
                                    self$update_annotations() 
                                 }
                              },
                              
                              set_knownRBPs_study_ids = function(knownRBPs_study_ids, update = TRUE) {
                                 stopifnot(is.character(knownRBPs_study_ids))
                                 self$RBPbaseDataHandler$assert_ids_in_rbpbase_for_organism(knownRBPs_study_ids, 
                                                                                            self$selected_studies_organism)
                                 
                                 self$knownRBPs_study_ids = knownRBPs_study_ids
                                 
                                 if(update){
                                    self$update_annotations() 
                                 }
                              },
                              
                              
                              update_annotations = function() {
                                 if(self$DEBUG){
                                    cat("update_annotations:\\n") 
                                 }
                                 
                                 if(self$has_knownRBPs_study_ids()) {
                                    if(self$DEBUG){
                                       cat("update_annotations: adding known annotation\\n") 
                                    }
                                    self$add_knownRBPs_annotation(ids = self$knownRBPs_study_ids)
                                 }
                                 
                                 if(self$has_enzyme_annotation() && self$has_metabolic_enzyme_annotation()) {
                                    if(self$DEBUG){
                                       cat("update_annotations: adding enzyme description annotation\\n") 
                                    }
                                    
                                    self$add_enzyme_annotation(enzyme_id = self$enzyme_id, metabolic_enzyme_id = self$metabolic_enzyme_id) 
                                 }
                                 
                                 if(self$has_knownRBPs_study_ids() && self$has_anyRIC_annotation()) {
                                    if(self$DEBUG){
                                       cat("RBPbaseAnalysis update_annotations: adding known vs RIC description annotation\\n") 
                                    }
                                    
                                    self$add_known_RIC_desc_annotation(knownRBPs_id = self$knownRBPs_id,
                                                                       anyRIC_id    = self$anyRIC_id)
                                    
                                 }
                              },
                              
                              
                              get_venn_data = function(selection_ids, names = NULL, filter = NULL) {
                                 stopifnot(is.character(selection_ids)) #&& length(selection_ids) <= 3)
                                 
                                 x <- self$get_main_table(selection_ids     = selection_ids, 
                                                          descriptive_names = is.null(names),
                                                          filter            = filter)
                                 
                                 
                                 
                                 if(!is.null(names)) {
                                    colnames(x) <- tibble(coln = colnames(x)) %>% 
                                       left_join(y = tibble(ids = selection_ids,
                                                            newn = names),
                                                 by = c("coln" = "ids")) %>%
                                       mutate(outn = ifelse(is.na(newn), coln, newn)) %>%
                                       pull(outn)
                                 }
                                 
                                 x %>% 
                                    mutate_if(is.logical, ~(ifelse(., UNIQUE, NA))) %>% 
                                    dplyr::select(-UNIQUE) %>% 
                                    as.list %>% 
                                    lapply(function(x) x[!is.na(x)]) %>%
                                    Venn()
                              },
                              
                              
                              plot_venn = function(selection_ids, names = NULL, ...) {
                                 stopifnot(is.character(selection_ids)&& length(selection_ids) <= 3)
                                 self$get_venn_data(selection_ids = selection_ids, 
                                                    names  = names) %>%
                                    ggvenn(...) +
                                    removeGrid()  + 
                                    theme(
                                       axis.line        = element_blank(),
                                       axis.text.x      = element_blank(),
                                       axis.text.y      = element_blank(),
                                       axis.ticks       = element_blank(),
                                       axis.title.x     = element_blank(),
                                       axis.title.y     = element_blank(),
                                       panel.background = element_blank(),
                                       panel.border     = element_blank(),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       plot.background  = element_blank()
                                    )  
                              },
                              
                              save_venn = function(selection_ids, 
                                                   names = NULL, 
                                                   file_root, 
                                                   formats = c(".png", ".pdf"), 
                                                   return_plot = T,
                                                   ...) {
                                 stopifnot(is.character(file_root) && length(file_root) == 1)
                                 stopifnot(is.character(formats) && length(file_root) > 0)
                                 
                                 p <- self$plot_venn(selection_ids = selection_ids, 
                                                     names  = names,
                                                     ...)
                                 
                                 for(i in formats) {
                                    cowplot::save_plot(p, 
                                                       filename = paste0(file_root, i), 
                                                       base_aspect_ratio = 1.5, 
                                                       base_height = ifelse(length(selection_ids) > 2, 8, 6))
                                 }
                                 
                                 if(return_plot) {
                                    return(p)
                                 } else {
                                    return(invisible(self)) 
                                 }
                              },
                              
                              
                              plot_euler = function(selection_ids, 
                                                    names                  = NULL,
                                                    fill                   = NULL,
                                                    filter                 = NULL,
                                                    quantities             = TRUE,
                                                    font                   = 1,
                                                    cex                    = 1,
                                                    alpha                  = 0.9,
                                                    descriptive_names      = F,
                                                    count_hits_for_studies = NULL,
                                                    numeric                = F,
                                                    annotation             = T,
                                                    studies                = T,
                                                    ...) {
                                 
                                 stopifnot(is.character(selection_ids)) #&& length(selection_ids) <= 3)
                                 if(!is.null(names)) {
                                    if (descriptive_names) {
                                       stop("descriptive_names has to be FALSE, when names argument is set")
                                    }
                                    if (length(names) != length(selection_ids)) {
                                       stop("Lengths of names has to be same length as selection_ids")  
                                    }
                                 }
                                 
                                 x <- self$get_main_table(selection_ids        = selection_ids, 
                                                          filter                 = filter,
                                                          descriptive_names      = descriptive_names,
                                                          count_hits_for_studies = count_hits_for_studies,
                                                          numeric                = numeric,
                                                          annotation             = annotation,
                                                          studies                = studies) %>% 
                                    dplyr::select(-UNIQUE)
                                 
                                 if(!is.null(names)) {
                                    x <- x %>% rename_at(vars(selection_ids), ~ names)
                                 }
                                 
                                 plot(euler(x),
                                      quantities = quantities,
                                      font = font,
                                      cex = cex,
                                      alpha = alpha, 
                                      ...)
                              },
                              
                              save_euler = function(file_root, 
                                                    formats = c(".png", ".pdf"), 
                                                    return_plot = T,
                                                    ...) {
                                 stopifnot(is.character(file_root) && length(file_root) == 1)
                                 stopifnot(is.character(formats) && length(file_root) > 0)
                                 
                                 p <- self$plot_euler(...)
                                 
                                 for(i in formats) {
                                    cowplot::save_plot(p, 
                                                       filename = paste0(file_root, i))
                                 }
                                 
                                 if(return_plot) {
                                    return(p)
                                 } else {
                                    return(invisible(self)) 
                                 }
                              },
                              
                              plot_upset_for_studies = function(selected_studies      = T,
                                                                comparison_studies    = T,
                                                                main_organism_studies = T,
                                                                n_intersections       = 50,
                                                                order_by              = "freq", 
                                                                ...) {
                                 ids <- self$get_studies_only_ids(selected_studies      = selected_studies,
                                                                  comparison_studies    = comparison_studies,
                                                                  main_organism_studies = main_organism_studies) 
                                 
                                 self$plot_upset(ids             = ids,
                                                 n_intersections = n_intersections,
                                                 order_by        = order_by,
                                                 ...)
                              },
                              
                              plot_upset = function(ids, 
                                                    n_intersections = 50,
                                                    order_by        = "freq",
                                                    ...) {
                                 self$get_main_table(selection_ids = ids,
                                                     ...) %>%
                                    gather(key   = "studies", 
                                           value = "found", 
                                           -UNIQUE) %>%
                                    filter(found) %>% 
                                    dplyr::select(studies, UNIQUE)%>%
                                    group_by(UNIQUE) %>%
                                    summarize(studies = list(studies)) %>%
                                    ggplot(aes(x = studies)) +
                                    geom_bar() +
                                    scale_x_upset(order_by        = order_by, 
                                                  n_intersections = n_intersections) 
                              },
                              
                              barplot = function(id, study_ids = NULL) {
                                 stopifnot(is.character(id))
                                 stopifnot(length(id) == 1)
                                 stopifnot(is.null(study_ids) || is.character(study_ids))
                                 
                                 if(is.null(study_ids)) {
                                    study_ids <- self$selected_studies_ids
                                 }
                                 
                                 
                                 self$get_main_table(selection_ids = c(id, study_ids)) %>% 
                                    gather(any_of(study_ids), 
                                           key   = "study", 
                                           value = "found") %>% 
                                    mutate(study =  self$RBPbaseDataHandler$get_descriptiveid_for_rbpbaseids(study)) %>% 
                                    filter(found)
                              },     
                              
                              plot_barplot = function(id, study_ids = NULL, col = NULL) {
                                 x <- self$barplot(id = id,
                                                   study_ids = study_ids)
                                 
                                 x <- x %>% 
                                    ggplot(aes(x = study, fill = !!sym(id))) + 
                                    geom_bar(position = "stack", width=0.5) +
                                    xlab("")  +  
                                    theme_minimal() + 
                                    theme(text        = element_text(size   = 20), 
                                          axis.text.x = element_text(size   = 20, 
                                                                     face   = "bold"), 
                                          plot.title  = element_text(size   = 25, 
                                                                     face   = "bold", 
                                                                     colour = "black", 
                                                                     vjust  = -1),
                                          legend.title = element_blank()
                                    )   + 
                                    removeGrid() 
                                 
                                 if(!is.null(fill)) 
                                    x <- x + scale_fill_manual(values = col)
                                 
                                 x + geom_text(stat      = 'count', 
                                               aes(label = ..count..),
                                               position  = position_stack(vjust = 0.5), 
                                               size      = 5.5) +
                                    ylab("") + 
                                    theme(plot.title  = element_text(hjust = 0.5),
                                          axis.text.y = element_blank(),
                                          axis.text.x = element_text(size=13), 
                                          axis.ticks  = element_blank(),
                                          axis.line   = element_blank())
                              },
                              
                              get_descriptiveid_for_rbpbaseids = function(...) {
                                 self$RBPbaseDataHandler$get_descriptiveid_for_rbpbaseids(...) 
                              },
                              
                              ## Functions for transformation
                              # Character string transformation for visualisation
                              # Known RBP or known RIC
                              transform_known = function(knownRBP, anyRIC) {
                                 stopifnot(is.logical(knownRBP))
                                 stopifnot(is.logical(anyRIC))
                                 
                                 unlist(
                                    lapply((as.numeric(knownRBP) * 10 + as.numeric(anyRIC)),
                                           function(x) {
                                              if(x == 11) { 
                                                 return("known & RIC")
                                              } else if (x == 10) {
                                                 return("known") 
                                              } else if (x == 1) {
                                                 return("RIC") 
                                              } else {
                                                 return("unknown")
                                              }
                                           })
                                 ) %>% factor(levels=c("known & RIC", "known", "RIC", "unknown"))
                              },
                              
                              ## Enzymes
                              transform_enzyme = function(Enzyme, Metabolic.Enzyme) {
                                 stopifnot(is.logical(Enzyme))
                                 stopifnot(is.logical(Metabolic.Enzyme))
                                 
                                 unlist(
                                    lapply((as.numeric(Enzyme) * 10 + as.numeric(Metabolic.Enzyme)),
                                           function(x) {
                                              if(x == 11) { 
                                                 return("metabolic enzyme")
                                              } else if (x == 10) {
                                                 return("enzyme") 
                                              } else if (x == 1) {
                                                 return("metabolic enzyme") 
                                              } else {
                                                 return("other")
                                              }
                                           })
                                 ) %>% factor(levels=c("other", "enzyme", "metabolic enzyme"))
                              }
                           )
)



new_RBPbaseAnalysis <- function(COMPILED_TABLE,
                                RIC_STUDIES,
                                ANNO_STUDIES,
                                selected_studies_ids,
                                comparison_studies_organism,
                                excluded_comparison_studies = NULL,
                                ...,
                                DEBUG = FALSE) {
   
   RBPbaseAnalysis$new(RBPbaseDataHandler = RBPbaseDataHandler$new(COMPILED_TABLE, 
                                                                   RIC_STUDIES,
                                                                   ANNO_STUDIES,
                                                                   DEBUG = DEBUG),
                       selected_studies_ids        = selected_studies_ids,
                       comparison_studies_organism = comparison_studies_organism,
                       excluded_comparison_studies = excluded_comparison_studies,
                       DEBUG                       = DEBUG,
                       ...)
}



