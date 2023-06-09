library(tidyverse)

# Import file with DOI strings in column DOI
df <- read_csv("data/nwo_dois.csv")

# Clean DOIs
df_dois <- df %>%
  #copy column DOI for cleaning (adapt original column name as needed)
  select(DOI_corrected) %>%
  mutate(DOI_clean = DOI)

#General doi cleaning
df_dois <- df_dois %>%
  #detect presence of actual doi string, to exclude e.g 'no DOI' and non-doi identifiers
  filter(str_detect(DOI_clean, '10\\\\.')) %>%
  # split strings on whitespace
  mutate(DOI_clean = str_split(DOI_clean, '\\\\s')) %>%
  # get each substring into a row
  unnest_longer(DOI_clean) %>%
  #filter on presence of doi string '10.' to remove split strings that are not dois
  filter(str_detect(DOI_clean, '10\\\\.')) %>%
  #remove everything before first '10.' to get rid of url prefixes
  mutate(DOI_clean = str_replace(DOI_clean, '^.*?10\\\\.', '10\\\\.')) %>%
  # remove trailing punctuation marks
  mutate(DOI_clean = str_remove(DOI_clean, '[[:punct:]]*$'))

#replace URL encoding for deduplication and downstream matching
df_dois <- df_dois %>%
  mutate(DOI_clean = str_replace_all(DOI_clean, '%2F', '\\\\/')) %>%
  mutate(DOI_clean = str_replace_all(DOI_clean, '%28', '(')) %>%
  mutate(DOI_clean = str_replace_all(DOI_clean, '%29', ')'))

# convert all dois to lower case for deduplication and downstream matching
df_dois <- df_dois %>%
  mutate(DOI_clean = str_to_lower(DOI_clean)) %>%
  distinct()

# match back to orginal df
df_join <- df %>%
  left_join(df_dois, by = "DOI")

#unique list of cleaned dois
df_unique <- df_dois %>%
  select(DOI_clean) %>%
  distinct()

#write to file
write_csv(df_join, "data/nwo_dois_cleaned.csv")
write_csv(df_unique, "data/nwo_dois_cleaned_unique.csv")