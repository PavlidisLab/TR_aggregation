## Save out two tables of literature curated targets. First is only those curated
## by Pavlab (Eric's paper as well as on-going curation). Second is all 
## interactions aggregated in Eric's paper (includes those from TRRUST etc) as 
## well as on-going curation. Note that on-going is based of googlesheets access
## https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009484
## -----------------------------------------------------------------------------

library(tidyverse)
library(googlesheets4)

date <- "July2022"  # last accessed July 4 2022
gsheets_id <- "1ngjKoRGaOgF-8BlxUPK7o7XRg7wimTxYYQkSokSVYUM"
out_all <- paste0("~/Data/Metadata/Curated_targets_all_", date, ".tsv")
out_pavlab <- paste0("~/Data/Metadata/Curated_targets_pavlab_", date, ".tsv")

# tables downloaded from supplement of eric's paper
lt_eric <- read.delim("~/Data/Metadata/Eric_resource_records_DTRI.tsv", stringsAsFactors = FALSE, skip = 1)
lt_all <- read.delim("~/Data/Metadata/Eric_resource_external_DTRI.tsv", stringsAsFactors = FALSE, skip = 1)

# protein coding genes for filtering
pc_hg <- read.delim("~/Data/Metadata/refseq_select_hg38.tsv", stringsAsFactors = FALSE)
pc_mm <- read.delim("~/Data/Metadata/refseq_select_mm10.tsv", stringsAsFactors = FALSE)

# on-going curation on googlesheets
current <- read_sheet(ss = gsheets_id, sheet = "Master_Curation", trim_ws = TRUE, col_types = "c")


# Format updated/current curation on googlesheets

current <- current[, 1:13]

current <- current %>% 
  dplyr::rename(TF_Symbol = TF_Gene_Name, 
                Target_Symbol = Target_Gene_Name) %>% 
  mutate(PubMed_ID = as.integer(PubMed_ID))


# Format the external/all table and join with Eric's resource to get relevant
# info when available. Note that external resource interactions are missing 
# details like pubmed ID or species


lt_all <- lt_all %>%
  dplyr::rename(TF_Symbol = TF_Symbol_Human, 
                Target_Symbol = Target_Symbol_Human) %>% 
  dplyr::select(DTRI_ID, TF_Symbol, Target_Symbol, Databases) %>% 
  left_join(., lt_eric, by = "DTRI_ID") %>% 
  mutate(TF_Species = str_to_title(TF_Species),
         Target_Species = str_to_title(Target_Species)) %>% 
  dplyr::select(-c(DTRI_ID))


# For current/updated resource, only keep/add relevant cols, and only keep
# protein coding genes. Add species info, which has been encoded by the casing
# of the symbol name

current <- current[, c("TF_Symbol", "Target_Symbol", "PubMed_ID", "Cell_Type", "Experiment_Type")]
current$Databases <- "Current"

current <- current %>% 
  filter(TF_Symbol %in% pc_hg$Symbol | TF_Symbol %in% pc_mm$Symbol) %>% 
  filter(Target_Symbol %in% pc_hg$Symbol | Target_Symbol %in% pc_mm$Symbol)

current$Target_Species <- current$TF_Species <- NA

for (i in 1:nrow(current)) {
  
  if (current$TF_Symbol[i] %in% pc_mm$Symbol) {
    current$TF_Species[i] <- "Mouse"
  } else if (current$TF_Symbol[i] %in% pc_hg$Symbol) {
    current$TF_Species[i] <- "Human"
  }
 
  if (current$Target_Symbol[i] %in% pc_mm$Symbol) {
    current$Target_Species[i] <- "Mouse"
  } else if (current$Target_Symbol[i] %in% pc_hg$Symbol) {
    current$Target_Species[i] <- "Human"
  }

}


# Format and join Eric resource with on going curation (Pavlab only curated)


lt_eric <- lt_eric %>%
  rename(TF_Symbol = TF_Symbol_Human,
         Target_Symbol = Target_Symbol_Human) %>%
  dplyr::select(one_of(intersect(colnames(current), colnames(.))))
  

pavlab_final <- rbind(lt_eric, dplyr::select(current, -Databases)) %>%
  mutate(
    TF_Species = str_to_title(TF_Species),
    Target_Species = str_to_title(Target_Species),
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

all_final <-
  rbind(lt_all[, intersect(colnames(current), colnames(lt_all))], current)



# Summarize how many non TRRUST/current unique targets are in batch 1

tfs <- c("pax6", "ascl1", "neurod1", "mecp2", "hes1", "runx1", "mef2c", "tcf4")

final_batch1 <- all_final %>% 
  filter(str_to_lower(TF_Symbol) %in% tfs) %>% 
  mutate(Other = !str_detect(Databases, "Current|TRRUST")) %>% 
  group_by(TF_Symbol) %>% 
  distinct(Target_Symbol, .keep_all = TRUE) %>% 
  ungroup()


# sum(final_batch1$Other)/nrow(final_batch1)  # 0.054


# Save


write.table(pavlab_final,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            file = out_pavlab)


write.table(all_final,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            file = out_all)

