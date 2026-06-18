suppressPackageStartupMessages({library(dplyr); library(tidyr); library(tibble)})
setwd("C:/Users/aneta/Documents/GitHub/avian_mixed_paternity/mixed_paternity_final")
data_overall_study <- readRDS("data/data_overall_study.rds")
cat("Counts with both nest variables:\n")
print(data_overall_study %>% summarise(
  n_total = n(),
  n_nest_box = sum(!is.na(nest_box_model)),
  n_nest_type = sum(!is.na(nest_type_model)),
  n_both = sum(!is.na(nest_box_model) & !is.na(nest_type_model))
))
cat("\nCross-tab:\n")
print(data_overall_study %>% filter(!is.na(nest_box_model), !is.na(nest_type_model)) %>% count(nest_type_model, nest_box_model) %>% group_by(nest_type_model) %>% mutate(row_percent = 100*n/sum(n)) %>% ungroup())
cat("\nNest box by type percent:\n")
print(data_overall_study %>% filter(!is.na(nest_box_model), !is.na(nest_type_model)) %>% count(nest_box_model, nest_type_model) %>% group_by(nest_box_model) %>% mutate(percent = 100*n/sum(n)) %>% ungroup())
