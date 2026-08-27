library(tidyverse)
library(arrow)
library(ggplot2)



diann_output <- "path/to/file.parquet"

# inspect all columns
read_parquet(diann_output) %>% names()

############################################################################################################################################################
# visualize peptide coverage of a protein

#compute protein coverage 



# BSA 
arrow::read_parquet(diann_output) %>% 
  separate(Run, into = c("Run_ID", "Run_number", "sample"), sep = "_", remove = F) %>% 
  unite("test", sep = "_" , c(Run_ID, Run_number, sample), remove = F) %>% 
  filter(
    grepl("BSA", Protein.Group)) %>% 

    mutate(Stripped.Sequence = fct_reorder(Stripped.Sequence, pept_start), 
         Proteotypic = factor(ifelse(Proteotypic == 1, "proteotypic", "not proteotypic"))) %>% 
  dplyr::rename("annotation" = Proteotypic) %>% 
  ggplot(aes(pept_start, Precursor.Quantity, text = Stripped.Sequence, text0 = Precursor.Quantity, text1 = PEP, text2 = pept_start, text3 = pept_stop, fill = annotation)) + 
  geom_rect(aes(xmin= pept_start, xmax = pept_stop, ymin = 0, ymax = Precursor.Quantity) , color = "#333333")+
  scale_fill_manual(values = c("proteotypic" = "#00FF00",
                               "not proteotypic"="#009900", 
                               "signal" = "#CC9900",
                               "linker" = "black", 
                               "chain A" = "#FFFFCC", 
                               "chain B" = "#FFCC33")) +
  facet_grid(prot_ID ~ test) + 
  theme_bw(base_size = 20)

