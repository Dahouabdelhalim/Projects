# Thomas Schwarzl
# schwarzl@embl.de
# This script compiles all figures and supplementary info

suppressPackageStartupMessages({
   require(rmarkdown)
   require(knitr)
   require(BiocStyle)
})


rmarkdown::render(
   input = "index.Rmd",
   output_format = NULL, 
   output_file = "index.html",
   output_dir = NULL,
   output_options = NULL,
   output_yaml = "index.yaml",
   intermediates_dir = NULL,
   knit_root_dir = NULL,
   runtime = "auto",
   clean = TRUE,
   params = NULL,
   knit_meta = NULL,
   envir = parent.frame(),
   run_pandoc = TRUE,
   quiet = TRUE,
   encoding = "UTF-8"
)

rmarkdown::render(
   input = "F1_Overlap.Rmd",
   output_format = NULL, 
   output_file = "F1_Overlap.html",
   output_dir = NULL,
   output_options = NULL,
   output_yaml = "F1_Overlap.yaml",
   intermediates_dir = NULL,
   knit_root_dir = NULL,
   runtime = "auto",
   clean = TRUE,
   params = NULL,
   knit_meta = NULL,
   envir = parent.frame(),
   run_pandoc = TRUE,
   quiet = TRUE,
   encoding = "UTF-8"
)

rmarkdown::render(
   input = "F2_Location.Rmd",
   output_format = NULL,
   output_file = "F2_Location.html",
   output_dir = NULL,
   output_options = NULL,
   output_yaml = "F2_Location.yaml",
   intermediates_dir = NULL,
   knit_root_dir = NULL,
   runtime = "auto",
   clean = TRUE,
   params = NULL,
   knit_meta = NULL,
   envir = parent.frame(),
   run_pandoc = TRUE,
   quiet = TRUE,
   encoding = "UTF-8"
)

rmarkdown::render(
   input = "F2_Location_alt.Rmd",
   output_format = NULL,
   output_file = "F2_Location_alt.html",
   output_dir = NULL,
   output_options = NULL,
   output_yaml = "F2_Location_alt.yaml",
   intermediates_dir = NULL,
   knit_root_dir = NULL,
   runtime = "auto",
   clean = TRUE,
   params = NULL,
   knit_meta = NULL,
   envir = parent.frame(),
   run_pandoc = TRUE,
   quiet = TRUE,
   encoding = "UTF-8"
)

rmarkdown::render(
   input = "F3_DiseaseNetworks.Rmd",
   output_format = NULL, 
   output_file = "F3_DiseaseNetworks.html",
   output_dir = NULL,
   output_options = NULL,
   output_yaml = "F3_DiseaseNetworks.yaml",
   intermediates_dir = NULL,
   knit_root_dir = NULL,
   runtime = "auto",
   clean = TRUE,
   params = NULL,
   knit_meta = NULL,
   envir = parent.frame(),
   run_pandoc = TRUE,
   quiet = TRUE,
   encoding = "UTF-8"
)


rmarkdown::render(
   input = "F4_TherapeuticArea.Rmd",
   output_format = NULL, 
   output_file = "F4_TherapeuticArea.html",
   output_dir = NULL,
   output_options = NULL,
   output_yaml = "F4_TherapeuticArea.yaml",
   intermediates_dir = NULL,
   knit_root_dir = NULL,
   runtime = "auto",
   clean = TRUE,
   params = NULL,
   knit_meta = NULL,
   envir = parent.frame(),
   run_pandoc = TRUE,
   quiet = TRUE,
   encoding = "UTF-8"
)


rmarkdown::render(
   input = "F5_List.Rmd",
   output_format = NULL, 
   output_file = "F5_List.html",
   output_dir = NULL,
   output_options = NULL,
   output_yaml = "F5_List.yaml",
   intermediates_dir = NULL,
   knit_root_dir = NULL,
   runtime = "auto",
   clean = TRUE,
   params = NULL,
   knit_meta = NULL,
   envir = parent.frame(),
   run_pandoc = TRUE,
   quiet = TRUE,
   encoding = "UTF-8"
)





