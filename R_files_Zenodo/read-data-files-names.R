## Copyright (c) 2018 Frédéric Vanwindekens, Dóra Drexler
## Centre wallon de Recherches agronomiques (CRA-W)
## Hungarian Research Institute of Organic Agriculture (ÖMKi),
## Projet DiverIMPACTS (https://www.diverimpacts.net/)
## License cc-by

data.full <- read.csv(
    "./data-raw/survey_438634_R_data_file-names.csv",
    quote = "'\\"",
    na.strings = c("", "\\"\\""),
    stringsAsFactors = FALSE,
    skip = 5)
