## Save out two tables of literature curated targets. First is only those curated
## by Pavlab (Eric's paper as well as on-going curation). Second is all 
## interactions aggregated in Eric's paper (includes those from TRRUST etc) as 
## well as on-going curation. Note that on-going requires googlesheets access.
## https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009484
## -----------------------------------------------------------------------------

library(tidyverse)
library(googlesheets4)
source("R/setup-01_config.R")


# Tables downloaded from supplement of Chu 2021 (all includes external dbs)
lt_chu2021 <- read.delim(chu2021_path_records, stringsAsFactors = FALSE, skip = 1)
lt_all <- read.delim(chu2021_path_all, stringsAsFactors = FALSE, skip = 1)

# Protein coding genes for filtering
pc_hg <- read.delim(ref_path_hg, stringsAsFactors = FALSE)
pc_mm <- read.delim(ref_path_mm, stringsAsFactors = FALSE)

# On-going curation on googlesheets
current <- read_sheet(ss = gsheets_curated, sheet = "Master_Curation", trim_ws = TRUE, col_types = "c", range = "A:M")


# Format and join targets curated in Chu2021 and aggregated from external dbs. 
# Note that external interactions are missing details like pubmed ID or species.
# ------------------------------------------------------------------------------


keep_cols <- c("DTRI_ID", 
               "TF_Symbol", 
               "Target_Symbol", 
               "TF_Species", 
               "Target_Species",
               "Cell_Type",
               "Experiment_Type",
               "PubMed_ID",
               "Databases")


lt_chu2021 <- lt_chu2021 %>% 
  dplyr::rename(
    TF_Symbol = TF_Symbol_Human,
    Target_Symbol = Target_Symbol_Human) %>% 
  dplyr::select(any_of(keep_cols)) %>% 
  mutate(
    TF_Species = str_to_title(TF_Species),
    Target_Species = str_to_title(Target_Species),
    Databases = "Chu2021"
  )


lt_all <- lt_all %>%
  dplyr::rename(
    TF_Symbol = TF_Symbol_Human,
    Target_Symbol = Target_Symbol_Human) %>%
  dplyr::select(any_of(keep_cols)) %>% 
  mutate(Databases = str_replace(Databases, "Current", "Chu2021"))


lt_all <-
  left_join(lt_all, lt_chu2021, by = "DTRI_ID", suffix = c("", ".y")) %>%
  select(-c(ends_with(".y"), "DTRI_ID")) %>%
  mutate(
    TF_Symbol = ifelse(
      TF_Species == "Mouse" & !is.na(TF_Species),
      str_to_title(TF_Symbol),
      TF_Symbol
    ),
    Target_Symbol = ifelse(
      Target_Species == "Mouse" & !is.na(Target_Species),
      str_to_title(Target_Symbol),
      Target_Symbol
    )
  )


# For current/updated resource, only keep/add relevant cols, and only keep
# protein coding genes. Add species info, which has been encoded by the casing
# of the symbol name
# ------------------------------------------------------------------------------


current <- current %>% 
  dplyr::rename(TF_Symbol = TF_Gene_Name, 
                Target_Symbol = Target_Gene_Name) %>% 
  mutate(PubMed_ID = as.integer(PubMed_ID)) %>% 
  dplyr::select(any_of(keep_cols)) %>% 
  filter(TF_Symbol %in% pc_hg$Symbol | TF_Symbol %in% pc_mm$Symbol) %>% 
  filter(Target_Symbol %in% pc_hg$Symbol | Target_Symbol %in% pc_mm$Symbol) %>% 
  mutate(
    Databases = "Chu2021",
    TF_Species = ifelse(TF_Symbol %in% pc_hg$Symbol, "Human", "Mouse"),
    Target_Species = ifelse(Target_Symbol %in% pc_hg$Symbol, "Human", "Mouse"))


# Format and join Chu2021 with on going curation (Pavlab only curated targets)


lt_pavlab <-
  rbind(lt_chu2021[, setdiff(keep_cols, "DTRI_ID")], current) %>%
  mutate(
    TF_Symbol = ifelse(
      TF_Species == "Mouse",
      str_to_title(TF_Symbol),
      str_to_upper(TF_Symbol)
    ),
    Target_Symbol = ifelse(
      Target_Species == "Mouse",
      str_to_title(Target_Symbol),
      str_to_upper(Target_Symbol)
    )
  )


# Join all interactions with on-going (includes external curated targets)


all_final <- rbind(lt_all, current)


# Remove TCF4 from non Pavlab sources as heavily contaminated with TCF7L2.
# Keep TCF4 itself, autoreg shown in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5243153/
# Also remove TCF4 targets from Pavlab that are actually TCF7L2.
# Also format the perturbation experiment type (absent from external)
# ------------------------------------------------------------------------------


rm_target <- c("lef1", "cd36", "sox9", "glce", "vegfa")


all_final <- all_final %>% 
  mutate(Rm = !str_detect(Databases, "Chu2021") & TF_Symbol == "TCF4") %>% 
  mutate(
    Rm = ifelse(str_to_lower(TF_Symbol) == "tcf4" & str_to_lower(Target_Symbol) == "tcf4", FALSE, Rm),
    Experiment_Type = case_when(
      str_detect(str_to_lower(Experiment_Type), ".*perturbation,*") ~ "Perturbation",
      str_detect(str_to_lower(Experiment_Type), ".*binding*") ~ "Binding",
      str_detect(str_to_lower(Experiment_Type), ".*reporter*") ~ "Reporter",
      TRUE ~ "NA"
    )
  ) %>% 
  filter(!Rm) %>%
  filter(!(str_to_lower(TF_Symbol) == "tcf4" & str_to_lower(Target_Symbol) %in% rm_target)) %>% 
  select(-Rm)


# Summarize how many non TRRUST/current unique targets are in batch 1
# ------------------------------------------------------------------------------


batch1_tfs <- c("pax6", "ascl1", "neurod1", "mecp2", "hes1", "runx1", "mef2c", "tcf4")

final_batch1 <- all_final %>% 
  filter(str_to_lower(TF_Symbol) %in% batch1_tfs) %>% 
  mutate(Other = !str_detect(Databases, "Chu2021|TRRUST")) %>% 
  group_by(TF_Symbol) %>% 
  distinct(Target_Symbol, .keep_all = TRUE) %>% 
  ungroup()


prop_orth <- sum(final_batch1$Other)/nrow(final_batch1)  # 0.047


# Save
# ------------------------------------------------------------------------------


write.table(lt_pavlab,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            file = curated_path_pavlab)


write.table(all_final,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            file = curated_path_all)
