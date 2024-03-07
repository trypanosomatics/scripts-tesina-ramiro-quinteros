# A cada data frame obtenido en el script 11, con los marcadores para cada muestra, le agrego la cantidad de pacientes en los que 
# esa proteína cambia a lo largo del conjunto de tiempos

#Se cargan los paquetes a utilizar
library(data.table)
library(dplyr)

#Se establece la carpeta donde se encuentran los datos de las proteína a usar obtenidas en el script 11
data_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/changing_markers"
setwd(data_folder)

#Se establece la carpeta donde se guardarán los resultados
output_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/changing_markers_final_w_number_of_patients_change"

# Se cargan los datos de la cantidad de pacientes en los que cada proteína cambia a lo largo del tiempo para el tiempo TMZ-D65
# (obtenidos en el scipt anterior)
markers_w_number_of_patients_D65_path <- sprintf("%s/markers_w_number_of_patients_D65_vs_TMZ.tsv", output_folder)
markers_w_number_of_patients_D65 <- fread(markers_w_number_of_patients_D65_path, header = T, sep = "\t", na.strings = NULL)

# Se cargan los datos de la cantidad de pacientes en los que cada proteína cambia a lo largo del tiempo para el tiempo TMZ-D65
# (obtenidos en el scipt anterior)
markers_w_number_of_patients_12M_path <- sprintf("%s/markers_w_number_of_patients_12M_vs_TMZ.tsv", output_folder)
markers_w_number_of_patients_12M <- fread(markers_w_number_of_patients_12M_path, header = T, sep = "\t", na.strings = NULL)


# Se obtiene la lista de archivos de proteínas a analizar al tiempo TMZ-D65
samples_D65_folder <- sprintf("%s/D65-TMZ", data_folder)
samples_D65 <-list.files(samples_D65_folder)

# Se obtiene la lista de archivos de proteínas a analizar al tiempo TMZ-12M
samples_12M_folder <- sprintf("%s/12M-TMZ", data_folder)
samples_12M <-list.files(samples_12M_folder)

# Para cada muestra en el tiempo TMZ-D65
for (sample_D65 in samples_D65) {
  sample_D65_path <- sprintf("%s/%s", samples_D65_folder, sample_D65)
  data_frame_sample_D65 <- fread(sample_D65_path, header = T, sep = "\t", na.strings = NULL)
  
  # Se fusionan los datos de la muestra con la información sobre la cantidad de pacientes en los que cada proteína cambia a lo largo del tiempo TMZ-D65
  final_markers_data_frame_D65 <- merge(data_frame_sample_D65,markers_w_number_of_patients_D65, by = "protein_identifier", all.x = TRUE)
  
  # Se seleccionan las columnas relevantes
  final_markers_data_frame_D65 <- final_markers_data_frame_D65 %>%
    select(protein_identifier, number_of_patients_w_change, number_of_patients_decrease, number_of_patients_increase, number_of_patients_w_no_change,
           mean_log2_fold_change, mean_log10_P_value, peptide_count, mean_signal_D65, mean_signal_TMZ, change_in_patients, 
           direction_of_change,sd_D65, sd_TMZ, sd_log2_Fold_change, sd_log10_P_value, patient_id)
  
  # Se ordenan los datos por la cantidad de pacientes en los que hay cambios
  final_markers_data_frame_D65 <- final_markers_data_frame_D65 %>%
    arrange(desc(number_of_patients_w_change))
  
  # Se obtiene el ID del paciente
  patient_D65 <- final_markers_data_frame_D65$patient_id[1]
  
  # Se guardan los datos para el conjunto de tiempos TMZ-D65
  final_markers_data_frame_D65_output_file <- sprintf("%s/D65-TMZ/%s_final_changing_markers_D65_vs_TMZ.tsv", output_folder, patient_D65)
  write.table(final_markers_data_frame_D65, file = final_markers_data_frame_D65_output_file, sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)
}

## El proceso anterior se repite para los tiempos TMZ y 12M, con sus respectivas operaciones

for (sample_12M in samples_12M) {
  sample_12M_path <- sprintf("%s/%s", samples_12M_folder, sample_12M)
  data_frame_sample_12M <- fread(sample_12M_path, header = T, sep = "\t", na.strings = NULL)
  
  final_markers_data_frame_12M <- merge(data_frame_sample_12M,markers_w_number_of_patients_12M, by = "protein_identifier", all.x = TRUE)
  
  final_markers_data_frame_12M <- final_markers_data_frame_12M %>%
    select(protein_identifier, number_of_patients_w_change, number_of_patients_decrease, number_of_patients_increase, number_of_patients_w_no_change,
           mean_log2_fold_change, mean_log10_P_value, peptide_count, mean_signal_12M, mean_signal_TMZ, change_in_patients, 
           direction_of_change, sd_12M, sd_TMZ, sd_log2_Fold_change, sd_log10_P_value, patient_id)
  
  final_markers_data_frame_12M <- final_markers_data_frame_12M %>%
    arrange(desc(number_of_patients_w_change))
  
  patient_12M <- final_markers_data_frame_12M$patient_id[1]
  
  final_markers_data_frame_12M_output_file <- sprintf("%s/12M-TMZ/%s_final_changing_markers_12M_vs_TMZ.tsv", output_folder, patient_12M)
  write.table(final_markers_data_frame_12M, file = final_markers_data_frame_12M_output_file, sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)
}