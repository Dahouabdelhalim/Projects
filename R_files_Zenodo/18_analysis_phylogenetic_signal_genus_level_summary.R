library(tidyverse)

# pagels lambda
li <- paste("output/", list.files("output", pattern = "genus_level_pagels"), sep = "")
lambda <- bind_rows(lapply(li, "read_csv"))

mean(lambda$lambda.lambda)
mean(lambda$lambda.P <= 0.05)

# Blombergs K

li <- paste("output/", list.files("output", pattern = "genus_level_blomberg"), sep = "")
lambda <- bind_rows(lapply(li, "read_csv"))

mean(lambda$BBK.K)
mean(lambda$BBK.P <= 0.05)
