## FIGURE S4 ##

# OENOTHEIN B, file: oenotheinB_herbivory.csv
friedman.test(oenotheinB_herbivory$abundance, oenotheinB_herbivory$treatment, oenotheinB_herbivory$species)

# QUERCETIN, file: quercetin_herbivory.csv
friedman.test(quercetin_herbivory$abundance, quercetin_herbivory$treatment, quercetin_herbivory$species)

# OENOTHEIN A, file: oenotheinA_herbivory.csv
friedman.test(oenotheinA_herbivory$abundance, oenotheinA_herbivory$treatment, oenotheinA_herbivory$species)

# KAEMPFEROL, file: kaempferol_herbivory.csv
friedman.test(kaempferol_herbivory$abundance, kaempferol_herbivory$treatment, kaempferol_herbivory$species)

# MYRICETIN, file: myricetin_herbivory.csv
friedman.test(myricetin_herbivory$abundance, myricetin_herbivory$treatment, myricetin_herbivory$species)

# NARINGENIN, file: naringenin_herbivory.csv
friedman.test(naringenin_herbivory$abundance, naringenin_herbivory$treatment, naringenin_herbivory$species)

# TOTAL, file: total_herbivory.csv
friedman.test(total_herbivory$abundance, total_herbivory$treatment, total_herbivory$species)
