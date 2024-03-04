# Se calcula la media, la variación absoluta y la relativa según las réplicas para cada péptido

# Se cargan los paquetes a utilizar
library(data.table)
library(dplyr)

# Se establece la carpeta donde están los datos de las señales crudas
data_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/Array Roche E1224/00_inputs/02_serums/raw_data"
setwd(data_folder)

# Se establece la carpeta donde están los datos del diseño del experimento
data_folder2 <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/design"

# Se establece la carpeta donde se guardará el resultado
output_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/raw_data_mean_replicas_w_variation/no_formula"


t_cruzi_map_data_file <- sprintf("%s/01_individual_serum_design_v3.tsv", data_folder2)
t_cruzi_map <- fread(t_cruzi_map_data_file, header = T, sep = "\t", na.strings = NULL)

# Se determinan los orígenes y suborígenes a usar
origins_to_use <- c("analyzed_proteome_cutoff_4_SD_2_pep",
                    "peaks_from_cardiacos_cutoff_4_SD_2_pep",
                    "peaks_from_serodiscordantes_cutoff_8_SD_2_pep",
                    "peaks_from_BLAST_against_231_pident_95_length_similarity_80",
                    "peaks_from_other_pathogens_cutoff_4_SD_2_pep")

suborigins_to_use <- c("","herpes_strain17", "herpes_HG52", "HCMV") # Como control

# Se filtran los orígenes y suborígenes a usar dentro de la totalidad de los datos
t_cruzi_map_filtered_by_origin <- t_cruzi_map[origin %in% origins_to_use]
t_cruzi_map_filtered_by_origin <- t_cruzi_map_filtered_by_origin[suborigin %in% suborigins_to_use]

# Se seleccionan los péptidos a utilizar en el análisis
peptides_to_use <- unique(t_cruzi_map_filtered_by_origin$truncated_peptide)

# Se obtiene la lista de archivos de muestras a analizar
samples<-list.files(data_folder)

# Se itera sobre cada muestra
for(sample in samples){
  
  # Se cargan los datos crudos del archivo
  raw_data <- fread(sample, header = T, sep = "\t", na.strings = NULL)
  
  # Se renombra la columna "type" como "time" para mayor claridad
  raw_data <- raw_data %>% 
    rename(time=type)
  
  # Se filtran los datos para seleccionar solo los péptidos a utilizar en el análisis
  raw_data <- raw_data[sequence %in% peptides_to_use]   
  
  # Se calcula la media de las señales de las réplicas, y la variación absoluta y relativa
  raw_data_mean <- raw_data %>%
  raw_data_mean <- raw_data %>%
    group_by(sequence,source,time) %>%
    summarize(replica_1_signal = mean(signal[replica == "1"]),
              replica_2_signal = mean(signal[replica == "2"]),
              mean_signal = mean(signal, na.rm = TRUE)) %>%
    mutate(relative_variation = (max(replica_1_signal, replica_2_signal) / min(replica_1_signal, replica_2_signal)),
           absolute_variation =  abs(replica_1_signal - replica_2_signal)) %>%
    ungroup()
  
  # Se guarda el resultado en un archivo de salida para cada muestra
  raw_data_mean_output_file <- sprintf("%s/%s.tsv", output_folder,sample)
  write.table(raw_data_mean, file = raw_data_mean_output_file, sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)
}