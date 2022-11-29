## Save out two tables of literature curated targets. First is only those curated
## by Pavlab (Eric's paper as well as on-going curation). Second is all 
## interactions aggregated in Eric's paper (includes those from TRRUST etc) as 
## well as on-going curation. Note that on-going is based of googlesheets access
## https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009484
## -----------------------------------------------------------------------------

library(tidyverse)
library(googlesheets4)

date <- "July2022"  # last accessed July 4 2022
gsheets_id <- "1ngjKoRGaOgF-8BlxUPK7o7XRg7wimTxYYQkSokSVYUM" # July
# gsheets_id <- "1PB2P-9Xk2zV0RSZnkY5BdnV6E4KkDpEKvo68Sw_Rnx8" # latest
out_all <- paste0("~/Data/Metadata/Curated_targets_all_", date, ".tsv")
out_pavlab <- paste0("~/Data/Metadata/Curated_targets_pavlab_", date, ".tsv")

# tables downloaded from supplement of eric's paper (all includes external dbs)
lt_eric <- read.delim("~/Data/Metadata/Eric_resource_records_DTRI.tsv", stringsAsFactors = FALSE, skip = 1)
lt_all <- read.delim("~/Data/Metadata/Eric_resource_external_DTRI.tsv", stringsAsFactors = FALSE, skip = 1)

# protein coding genes for filtering
pc_hg <- read.delim("~/Data/Metadata/refseq_select_hg38.tsv", stringsAsFactors = FALSE)
pc_mm <- read.delim("~/Data/Metadata/refseq_select_mm10.tsv", stringsAsFactors = FALSE)

# on-going curation on googlesheets
current <- read_sheet(ss = gsheets_id, sheet = "Master_Curation", trim_ws = TRUE, col_types = "c", range = "A:M")

# Format and join Eric's curated + external. Note that external resource 
# interactions are missing details like pubmed ID or species

keep_cols <- c("DTRI_ID", 
               "TF_Symbol", 
               "Target_Symbol", 
               "TF_Species", 
               "Target_Species",
               "Cell_Type",
               "Experiment_Type",
               "PubMed_ID",
               "Databases")


lt_eric <- lt_eric %>% 
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
  left_join(lt_all, lt_eric, by = "DTRI_ID", suffix = c("", ".y")) %>%
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


# Format and join Eric's resource with on going curation (Pavlab only curated)


lt_pavlab <-
  rbind(lt_eric[, setdiff(keep_cols, "DTRI_ID")], current) %>%
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


# Remove TCF4 from non Pavlab sources as sampled examples of TCF7L2 contamination.
# Keep TCF4 itself, autoreg shown in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5243153/
# Also remove TCF4 targets from Pavlab that are actually TCF7L2

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


tfs <- c("pax6", "ascl1", "neurod1", "mecp2", "hes1", "runx1", "mef2c", "tcf4")

final_batch1 <- all_final %>% 
  filter(str_to_lower(TF_Symbol) %in% tfs) %>% 
  mutate(Other = !str_detect(Databases, "Chu2021|TRRUST")) %>% 
  group_by(TF_Symbol) %>% 
  distinct(Target_Symbol, .keep_all = TRUE) %>% 
  ungroup()


# sum(final_batch1$Other)/nrow(final_batch1)  # 0.047


# Save


write.table(lt_pavlab,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            file = out_pavlab)


write.table(all_final,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            file = out_all)
