#Load libraries
library(biomaRt)
library(dplyr)

#1. Extract the data from the table with the TFs and the table with the dN/dS values.
human_dNdS <- read.table("HumanDnDsW.txt", header = TRUE, sep = "\t")
tf_data <- read.table("HumanTFs_DBD.txt", header = TRUE, sep = "\t")
tf_data <- tf_data[, 2:4] 


#2. Start Ensembl Mart.
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://oct2014.archive.ensembl.org")
listDatasets(ensembl)
attributes <- listAttributes(ensembl)
listFilters(ensembl)


#3. Create a list with all the TFs (that were in HumanTFs_DBD.txt).
tf_ids <- tf_data$Ensembl_ID 


#4. Identify genes that were not already in the dNdS table and create a new dataframe with these new genes.
new_rows <- data.frame(
  ensembl_gene_id = setdiff(tf_ids, human_dNdS$GeneID), 
  Homo_sapiens.dN = NA, 
  Homo_sapiens.dS = NA,
  Homo_sapiens.dw = NA, 
  Homo_sapiens.w = NA
)


#5. Put the new data frame below with the names of the TFs that were not in human_dNdS.
colnames(human_dNdS)[1] <- "ensembl_gene_id" #Rename the GeneID column in human_dNdS to match the new rows.
human_dNdS <- rbind(human_dNdS, new_rows)
human_dNdS <- human_dNdS %>%
  group_by(ensembl_gene_id) %>%
  summarise(across(everything(), ~ {
    # Remove NA values and collapse
    non_na_values <- na.omit(.)
    if (length(non_na_values) == 0) {
      return(NA)  # If all values are NA, return NA
    } else {
      return(non_na_values[1])  # Choose the first non-NA value
    }
  }))


#6. Download from Ensembl.
name_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                   filters = 'ensembl_gene_id',
                   values = human_dNdS$ensembl_gene_id,
                   mart = ensembl)

human_dNdS = full_join(name_info, human_dNdS, by = "ensembl_gene_id")


#7. Define a list of species to compare the dNdS with.
species_list <- list(
  #Human = c("hsapiens_paralog_ds", "hsapiens_paralog_dn"),
  Chimpanzee = c("ptroglodytes_homolog_ds", "ptroglodytes_homolog_dn"),
  #Rat = c("rnorvegicus_homolog_ds", "rnorvegicus_homolog_dn"),
  Mouse = c("mmusculus_homolog_ds", "mmusculus_homolog_dn")
  #Macaque = c("mmulatta_homolog_ds", "mmulatta_homolog_dn"),
  #Marmoset = c("cjacchus_homolog_ds", "cjacchus_homolog_dn")
)


#8. Crear un data frame vacío para guardar los resultados
result_df <- data.frame(ensembl_gene_id = human_dNdS$ensembl_gene_id) #es un dataframe vacío con todos los ids de la tabla human_dNdS


#9. Extraer de Ensembl el dN y dS de todas las especies de la lista
for (species in names(species_list)) {
  species_info <- getBM(
    attributes = c('ensembl_gene_id', species_list[[species]]),
    filters = 'ensembl_gene_id',
    values = human_dNdS$ensembl_gene_id,
    mart = ensembl
  )
  # Unir species_info en el dataframe result_df mediante su ensembl_gene_id
  result_df <- merge(result_df, species_info, by = 'ensembl_gene_id')
}


#10. Remplazar los 0 con NA en result_df
result_df[result_df == 0] <- NA
result_df <- result_df %>%
  group_by(ensembl_gene_id) %>%
  summarise(across(everything(), ~ {
    # Remove NA values and collapse
    non_na_values <- na.omit(.)
    if (length(non_na_values) == 0) {
      return(NA)  # If all values are NA, return NA
    } else {
      return(non_na_values[1])  # Choose the first non-NA value
    }
  }))


#11. Join the result_df data frame (which has the dNdS of the homologs and paralogs) with the one we already had with the human dNdS.
dNdS = full_join(human_dNdS, result_df, by = "ensembl_gene_id")


#12. Export it to an excel.
library(writexl)
write_xlsx(dNdS, path = "dNdS_chimp&mouse.xlsx")

library(readxl)
df <- read_excel("dNdS_chimp&mouse_w.xlsx", sheet = "Sheet1")



###################################### GRAPHICS ############################################

#Load libraries
library(ggplot2)
library(ggstatsplot)
library(hrbrthemes)
library(viridis)
library(dplyr)

#Getting OLIG2 values
mm_result <- df %>%
  filter(hgnc_symbol == "OLIG2") %>%
  pull(mmusculus_homolog_w)

pt_result <- df %>%
  filter(hgnc_symbol == "OLIG2") %>%
  pull(ptroglodytes_homolog_w)


#1. All genes plot.
t_test <- t.test(df$ptroglodytes_homolog_w, mu = pt_result)
p_value <- t_test$p.value

ggplot(data = df, aes(x = "", y = ptroglodytes_homolog_w)) +
  geom_violin(fill = '#cad197') +
  geom_boxplot(width = 0.1, color = "gray", alpha = 0.2) +
  geom_point(aes(x = "", y = pt_result), 
             color = "red", size = 2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position = "none",
    plot.title = element_text(size=11)
  ) + ylim(c(0,1))
#repeat for mus musculus.


#2. TFs plot.
df_filtered <- df %>%
  filter(ensembl_gene_id %in% tf_ids)

write_xlsx(df_filtered, path = "df_filtered.xlsx")

t_test <- t.test(df_filtered$ptroglodytes_homolog_w, mu = pt_result)
p_value <- t_test$p.value

ggplot(data = df_filtered, aes(x = "", y = ptroglodytes_homolog_w)) +
  geom_violin(fill = "#cad197") +
  geom_boxplot(width = 0.1, color = "gray", alpha = 0.2) +
  geom_point(aes(x = "", y = pt_result), 
             color = "red", size = 2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position = "none",
    plot.title = element_text(size=11)
  ) + ylim(c(0,1))
#repeat for mus musculus.


#3. bHLH plot.
tf_data <- read.table("HumanTFs_DBD.txt", header = TRUE, sep = "\t")
tf_data <- tf_data[, 2:4] 

filt <- tf_data %>%
  filter(DBD == "bHLH") %>%
  pull(Ensembl_ID)

df_filtered <- df %>%
  filter(ensembl_gene_id %in% filt)

t_test <- t.test(df_filtered$ptroglodytes_homolog_w, mu = pt_result)
p_value <- t_test$p.value

ggplot(data = df_filtered, aes(x = "", y = ptroglodytes_homolog_w)) +
  geom_violin(fill = "#cad197") +
  geom_boxplot(width = 0.1, color = "gray", alpha = 0.2) +
  geom_point(aes(x = "", y = pt_result), 
             color = "red", size = 2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position = "none",
    plot.title = element_text(size=11)
  ) + ylim(c(0,1))
#repeat for mus musculus.
