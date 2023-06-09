
library(RPostgreSQL)
library(graph)
library(igraph)
library(tidyverse)
library(ggplot2)
library(stringi)
library(compiler)
library(gsubfn)

#########
# @author: Scott A. Malec, PhD University of Pittsburgh School of Medicine DBMI
# @date: 2021/08/09
# @lastrevised: 2021/08/09
#########

#########
# DB setup
#########
drv <- dbDriver('PostgreSQL')
con <-
  dbConnect(
    drv,
    dbname = "smdb",
    host = "localhost",
    port = 5432,
    user = "smalec",
    password = "mandarin"
  )

##############
# This script requires the entirety of SemMedDB to be imported into Postgres.
#

#########
# Knowledge characteristics
#########
# * queries informing the PRISMA-style diagram for the relevant predications used from SemMedDB

##
# How much total knowledge [(top) box0 {zero}]
totalPredications <-dbGetQuery(con, paste(" select count(*) from predication;", sep = ""))
print(paste("totalPredications: ", as.integer(totalPredications)))

##
# How many causal predications? [ box1 = totalCausalPredications]
#                               [ box1_ex = totalPredications - totalCausalPredications ]
#predicates <- " ('CAUSES' ,'TREATS', 'PREVENTS', 'PREDISPOSES') "
predicates <- " ('CAUSES', 'PREDISPOSES', 'PREVENTS', 'AFFECTS', 'DISRUPTS', 'AUGMENTS') "

print(paste("Predicates used: ", predicates))

totalCausalPredications <- dbGetQuery(con, paste(" WITH ZX AS (SELECT cp.subject_cui AS covar, COUNT(*) AS scnt, cp.subject_name
                        FROM predication cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND cp.object_cui IN ('C0011581', 'C1269683', 'C0871546', 'C0041696', 'C0588008', 'C0588006', 'C0221480', 'C0154409', 'C0024517')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY subject_cui, subject_name),
                        ZY AS (SELECT cp.subject_cui AS covar , COUNT(*) AS scnt, cp.subject_name
                        FROM predication cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND cp.object_cui IN ('C0002395', 'C0494463', 'C0276496')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY cp.subject_cui, cp.subject_name)
                        SELECT COUNT(*) --DISTINCT ZX.covar, zx.scnt as theta1, zy.scnt AS theta2, zx.subject_name 
                                                 FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar; --ORDER BY theta2 desc, theta1 desc;
                                                 ", sep = ""))

print(paste("totalCausalPredications: ", totalCausalPredications))

box1_ex = totalPredications - totalCausalPredications

print(paste("totalPredications - totalCausalPredications: ", box1_ex))
  
#########
# get confounders
#########
getConfounders <- function(predicates) {
  confounders <-dbGetQuery(con, paste(" WITH ZX AS (SELECT cp.subject_cui AS covar, COUNT(*) AS scnt, cp.subject_name
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(c.pyear AS NUMERIC) >= 2010
                        AND cp.object_cui IN ('C0011581', 'C1269683', 'C0871546', 'C0041696', 'C0588008', 'C0588006', 'C0221480', 'C0154409', 'C0024517')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY subject_cui, subject_name),
                        ZY AS (SELECT cp.subject_cui AS covar , COUNT(*) AS scnt, cp.subject_name
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(pyear AS NUMERIC) >= 2010
                        AND cp.object_cui IN ('C0002395', 'C0494463', 'C0276496')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY cp.subject_cui, cp.subject_name)
                        SELECT DISTINCT ZX.covar, zx.scnt as theta1, zy.scnt AS theta2, zx.subject_name FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta2 desc, theta1 desc;", sep = ""))
  write.table(file = "covariates/rawConfounders_for_depression2AD.txt", x = confounders, row.names = FALSE, col.names = FALSE, quote = FALSE)  
  return(confounders)
}


#




#########
# DB setup
#########



#########
# DB setup
#########


#########
# DB setup
#########


#########
# DB setup
#########


#########
# DB setup
#########





