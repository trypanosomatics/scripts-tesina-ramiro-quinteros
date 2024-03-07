# Se crea un data frame por conjunto de tiempos que establece para cada proteína la cantidad de pacientes en las que cambia 
# a lo largo del tiempo (o no lo hace), y la dirección del mismo

#Se cargan los paquetes a utilizar
library(data.table)
library(dplyr)

#Se establece la carpeta donde se encuentran los datos de las proteína a usar obtenidas en el script anterior
data_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/all_antigens"
setwd(data_folder)

#Se establece la carpeta donde se guardarán los resultados
output_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/changing_markers_final_w_number_of_patients_change"

# Se obtiene la lista de archivos de proteínas a analizar al tiempo TMZ-D65
samples_D65_folder <- sprintf("%s/D65-TMZ", data_folder)
samples_D65 <-list.files(samples_D65_folder)

# Se obtiene la lista de archivos de proteínas a analizar al tiempo TMZ-12M
samples_12M_folder <- sprintf("%s/12M-TMZ", data_folder)
samples_12M <-list.files(samples_12M_folder)

# Se leen y combinan los datos de todas las proteínas para el tiempo TMZ-D65
combined_df_D65 <- data.frame()
for (sample_D65 in samples_D65) {
  sample_D65_path <- sprintf("%s/%s", samples_D65_folder, sample_D65)
  combined_df_D65_aux <- fread(sample_D65_path, header = T, sep = "\t", na.strings = NULL)
  combined_df_D65 <- rbind(combined_df_D65,combined_df_D65_aux)
}

# Se leen y combinan los datos de todas las proteínas para el tiempo TMZ-12M
combined_df_12M <- data.frame()
for (sample_12M in samples_12M) {
  sample_12M_path <- sprintf("%s/%s", samples_12M_folder, sample_12M)
  combined_df_12M_aux <- fread(sample_12M_path, header = T, sep = "\t", na.strings = NULL)
  combined_df_12M <- rbind(combined_df_12M,combined_df_12M_aux)
}


# Se calcula la cantidad de pacientes en los que cada proteína experimenta cambios a lo largo del tiempo o permanece sin cambios,
# así como la dirección del cambio para el tiempo TMZ-D65
markers_w_number_of_patients_D65 <- combined_df_D65 %>%
  group_by(protein_identifier) %>%
  summarise(
    number_of_patients_w_change = sum(change_in_patients == "Change", na.rm = TRUE),
    number_of_patients_w_no_change = sum(change_in_patients == "No_change", na.rm = TRUE),
    number_of_patients_increase = sum(direction_of_change == "Increase", na.rm = TRUE),
    number_of_patients_decrease = sum(direction_of_change == "Decrease", na.rm = TRUE)
  )

# Se calcula la cantidad de pacientes en los que cada proteína experimenta cambios a lo largo del tiempo o permanece sin cambios,
# así como la dirección del cambio para el tiempo TMZ-12M
markers_w_number_of_patients_12M <- combined_df_12M %>%
  group_by(protein_identifier) %>%
  summarise(
    number_of_patients_w_change = sum(change_in_patients == "Change", na.rm = TRUE),
    number_of_patients_w_no_change = sum(change_in_patients == "No_change", na.rm = TRUE),
    number_of_patients_increase = sum(direction_of_change == "Increase", na.rm = TRUE),
    number_of_patients_decrease = sum(direction_of_change == "Decrease", na.rm = TRUE)
  )

# Se ordenan los datos, la cantidad de pacientes en los que hay más cambios aparecen primeros
markers_w_number_of_patients_D65 <- markers_w_number_of_patients_D65 %>%
  arrange(desc(number_of_patients_w_change))

markers_w_number_of_patients_12M <- markers_w_number_of_patients_12M %>%
  arrange(desc(number_of_patients_w_change))

# Se guardan los datos para el conjunto de tiempos TMZ-D65
markers_w_number_of_patients_D65_output_file <- sprintf("%s/markers_w_number_of_patients_D65_vs_TMZ.tsv", output_folder)
write.table(markers_w_number_of_patients_D65, file = markers_w_number_of_patients_D65_output_file, sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)

# Se guardan los datos para el conjunto de tiempos TMZ-12M
markers_w_number_of_patients_12M_output_file <- sprintf("%s/markers_w_number_of_patients_12M_vs_TMZ.tsv", output_folder)
write.table(markers_w_number_of_patients_12M, file = markers_w_number_of_patients_12M_output_file, sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)



