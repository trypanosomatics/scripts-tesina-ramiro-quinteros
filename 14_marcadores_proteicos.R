#Se cargan los paquetes a utilizar
library(data.table)
library(dplyr)
library(tidyr)

#Se establece la carpeta donde se encuentran los datos de las señales suavizadas
data_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/data_w_rolling"
setwd(data_folder)

#Se establece la carpeta donde se encuentran los datos de las señales sin procesar
data_folder2 <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/raw_data_mean_replicas_filtered_w_formula"

#Se establece la carpeta donde se guardarán los marcadores de cambio
output_folder1 <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/changing_markers"

#Se establece la carpeta donde se guardarán todos los antígenos
output_folder2 <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/all_antigens"

# Se establece el sufijo de los archivos
suffix <- "cruzi_protein_data.tsv"
suffix2 <- "raw_data_filtered.tsv"

# Se especifica el paciente o listas de pacientes que se van a utilizar
# patients<- "2-159-304"

patients<- c("1-005-101","1-049-116","1-124-150","1-160-163","1-227-184","1-264-203","2-098-277","2-132-293",
             "2-190-322","2-212-332","2-240-353","2-261-364", "1-052-118","1-105-143","1-143-157",
             "1-218-183","1-248-193","2-011-253","2-103-280","2-125-286", "2-159-304","2-177-313","2-203-331",
             "2-242-355","1-044-115","1-063-121","1-064-124","1-110-146","1-130-153","1-150-159","1-153-161","1-225-187",
             "1-245-195","2-032-256","2-211-335","2-236-352")

# Bucle principal que itera sobre cada paciente
for (patient in patients) {
  
  # Se determinan los archivos donde se encuentran los datos suavizados
  smoothed_data_file_TMZ <-  sprintf("%s/%s_TMZ_%s", data_folder, patient, suffix)
  smoothed_data_file_D65 <-  sprintf("%s/%s_D65_%s", data_folder, patient, suffix)
  smoothed_data_file_12M <-  sprintf("%s/%s_12M_%s", data_folder, patient, suffix)
  
  #  Se cargan los datos suavizados
  smoothed_data_TMZ <- fread(smoothed_data_file_TMZ, header = T, sep = "\t", na.strings = NULL)
  smoothed_data_D65 <- fread(smoothed_data_file_D65, header = T, sep = "\t", na.strings = NULL)
  smoothed_data_12M <- fread(smoothed_data_file_12M, header = T, sep = "\t", na.strings = NULL)
  
  # Se unen los datos suavizados de los tres tiempos en un único conjunto de datos
  smoothed_data <- rbind(smoothed_data_TMZ,smoothed_data_D65,smoothed_data_12M)
  
  # Se determinan los archivos donde se encuentran los sin procesar
  raw_data_file_TMZ <-  sprintf("%s/%s_TMZ_%s", data_folder2, patient, suffix2)
  raw_data_file_D65 <-  sprintf("%s/%s_D65_%s", data_folder2, patient, suffix2)
  raw_data_file_12M <-  sprintf("%s/%s_12M_%s", data_folder2, patient, suffix2)
  
  #  Se cargan los datos sin procesar
  raw_data_TMZ <- fread(raw_data_file_TMZ, header = T, sep = "\t", na.strings = NULL)
  raw_data_D65 <- fread(raw_data_file_D65, header = T, sep = "\t", na.strings = NULL)
  raw_data_12M <- fread(raw_data_file_12M, header = T, sep = "\t", na.strings = NULL)
  
  # Se unen los datos suavizados sin procesar
  raw_data <- rbind(raw_data_TMZ,raw_data_D65,raw_data_12M)
  
  # Se ordenan los datos
  smoothed_data <- smoothed_data[order(source, time, sequence)]
  raw_data <- raw_data[order(source, time, sequence)]
  
  # Se filtran los datos suavizados para cada tiempo
  smoothed_data1 <- filter(smoothed_data, time == "TMZ" | time == "D65")
  smoothed_data2 <- filter(smoothed_data, time == "TMZ" | time == "12M")
  
  # Se filtran los datos sin procesar para cada tiempo
  raw_data1 <- filter(raw_data, time == "TMZ" | time == "D65")
  raw_data2 <- filter(raw_data, time == "TMZ" | time == "12M")
  
  # Se sacan todos los datos que no tengan 2 tiempos
  smoothed_data1 <- smoothed_data1 %>%
    group_by(sequence, source) %>%
    mutate(time_count = n_distinct(time)) %>%
    filter(time_count == 2) %>%
    ungroup() %>%
    select(-time_count)
  
  # Se seleccionan las columnas a utilizar, y se renombra una de ellas
  smoothed_data1 <- smoothed_data1 %>%
    select(protein, time, truncated_sequence, start, source, mean_smoothed_signal) %>%
    rename(sequence = truncated_sequence)
  
  
  sub_smoothed_data1 <- as.data.table(smoothed_data1)
  
  # Se sacan todos los datos que no tengan 2 tiempos
  raw_data1 <- raw_data1 %>%
    group_by(sequence, source) %>%
    mutate(time_count = n_distinct(time)) %>%
    filter(time_count == 2) %>%
    ungroup() %>%
    select(-time_count)
  
  sub_raw_data1 <- as.data.table(raw_data1)
  
  # Se fusionan los datos suavizados y sin procesar para el primer conjunto de tiempos
  all_data1 <- merge(sub_smoothed_data1,sub_raw_data1[,c("sequence","time","replica_1_signal","replica_2_signal")], by = c("sequence","time"), all.x = TRUE)
  
  
  sub_all_data1 <- as.data.table(all_data1)
  # unique_times1 <- unique(sub_all_data1$time)
  unique_times1 <- c("D65","TMZ")
  
  # Se realiza un pivot de los datos para crear un formato más adecuado para realizar las demás operaciones (de un formato largo a ancho)
  plot_data1 <- sub_all_data1 %>%
    pivot_wider(
      id_cols = c(sequence, protein, start, source),
      names_from = time,
      values_from = c(mean_smoothed_signal, replica_1_signal, replica_2_signal),
      names_glue = "{.value}_{time}"
    )

  # Se renombran las columnas
  plot_data1 <- plot_data1 %>%
    rename(protein_identifier = protein)
  
  # Se eliminan filas con valores NA
  plot_data1 <- na.omit(plot_data1)
  
  # Se realiza una prueba t de student sobre los datos y se calcula  el p-value
  plot_data1 <- plot_data1 %>%
    rowwise() %>%
    mutate(p_value = t.test(c(mean_smoothed_signal_D65, replica_1_signal_D65, replica_2_signal_D65), c(mean_smoothed_signal_TMZ, replica_1_signal_TMZ, replica_2_signal_TMZ), paired = F)$p.value) %>%
    ungroup()
  
  plot_data1 <- as.data.table(plot_data1)
  
  # Se calcula el logaritmo del fold change y el log del p-value
  plot_data1$log2_Fold_change <- log2(plot_data1$mean_smoothed_signal_D65) - log2(plot_data1$mean_smoothed_signal_TMZ)
  plot_data1$log10_P_value <- -log10(plot_data1$p_value)
  
  # Se eliminan filas duplicadas en el caso de que las haya
  plot_data1 <- plot_data1 %>%
    distinct()
  
  
  ## El proceso anterior se repite para los tiempos TMZ y 12M, con sus respectivas operaciones
  smoothed_data2 <- smoothed_data2 %>%
    group_by(sequence, source) %>%
    mutate(time_count = n_distinct(time)) %>%
    filter(time_count == 2) %>%
    ungroup() %>%
    select(-time_count)
  
  smoothed_data2 <- smoothed_data2 %>%
    select(protein, time, truncated_sequence, start, source, mean_smoothed_signal) %>%
    rename(sequence = truncated_sequence)
  
  sub_smoothed_data2 <- as.data.table(smoothed_data2)
  
  
  raw_data2 <- raw_data2 %>%
    group_by(sequence, source) %>%
    mutate(time_count = n_distinct(time)) %>%
    filter(time_count == 2) %>%
    ungroup() %>%
    select(-time_count)
  
  sub_raw_data2 <- as.data.table(raw_data2)
  
  all_data2 <- merge(sub_smoothed_data2,sub_raw_data2[,c("sequence","time","replica_1_signal","replica_2_signal")], by = c("sequence","time"), all.x = TRUE)
  
  sub_all_data2 <- as.data.table(all_data2)
  # unique_times2 <- unique(sub_all_data2$time)
  unique_times2 <- c("12M","TMZ")
  
  plot_data2 <- sub_all_data2 %>%
    pivot_wider(
      id_cols = c(sequence, protein, start, source),
      names_from = time,
      values_from = c(mean_smoothed_signal, replica_1_signal, replica_2_signal),
      names_glue = "{.value}_{time}"
    )
  
  plot_data2 <- plot_data2 %>%
    rename(protein_identifier = protein)
  
  plot_data2 <- na.omit(plot_data2)
  
  plot_data2 <- plot_data2 %>%
    rowwise() %>%
    mutate(p_value = t.test(c(mean_smoothed_signal_12M, replica_1_signal_12M, replica_2_signal_12M), c(mean_smoothed_signal_TMZ, replica_1_signal_TMZ, replica_2_signal_TMZ), paired = F)$p.value) %>%
    ungroup()
  
  plot_data2 <- as.data.table(plot_data2)
  
  plot_data2$log2_Fold_change <- log2(plot_data2$mean_smoothed_signal_12M) - log2(plot_data2$mean_smoothed_signal_TMZ)
  plot_data2$log10_P_value <- -log10(plot_data2$p_value)
  
  plot_data2 <- plot_data2 %>%
    distinct()

  # Se fusionan los datos de plot_data1 y plot_data2 basados en la secuencia, la protein ID, y el punto de inicio,
  # manteniendo solo las columnas de plot_data2 relacionadas con el tiempo 12M, para poder tener las señales de este tiempo y realizar
  # el cálculo de los picos antigénicos pertenecientes a algún tiempo
  markers1_w_id <- merge(plot_data1,plot_data2[ ,c("sequence", "protein_identifier", "start", "mean_smoothed_signal_12M")], by = c("sequence","protein_identifier", "start"), all.x = TRUE, all.y=FALSE)
  
  # Se filtran los datos para mantener solo las filas donde al menos una de las señales (TMZ, D65, o 12M) es mayor o igual a 10000
  # Esto son los picos antigénicos
  markers1_w_id <- markers1_w_id %>%
  markers1_w_id <- markers1_w_id %>%
    filter(mean_smoothed_signal_TMZ >= 10000 | mean_smoothed_signal_D65 >= 10000 | mean_smoothed_signal_12M >= 10000)
  
  # Se agrupan los datos por la protein ID y se calculan diversas estadísticas: la media y la desviación estándar de diversas
  # métricas: log2(fold change), log(p-value), señal a los diferentes tiempos. Se cuenta la cantidad de péptidos pertenecientes al pico
  markers1_w_id <- markers1_w_id %>%
    group_by(protein_identifier) %>%
    summarize(mean_log2_fold_change = round(mean(log2_Fold_change),2), mean_log10_P_value = round(mean(log10_P_value),2),
              mean_signal_D65 = round(mean(mean_smoothed_signal_D65),2), mean_signal_TMZ = round(mean(mean_smoothed_signal_TMZ),2), 
              peptide_count = n(), sd_D65 = round(sd(mean_smoothed_signal_D65),2), sd_TMZ = round(sd(mean_smoothed_signal_TMZ),2),
              sd_log2_Fold_change = round(sd(log2_Fold_change),2),sd_log10_P_value = round(sd(log10_P_value),2))
  
  # Se agregan nuevas columnas para indicar si hay un cambio significativo entre tiempos TMZ-D65 en los pacientes y la dirección del cambio
  markers1_w_id <- markers1_w_id %>%
    mutate(
      change_in_patients = case_when(
        (mean_log2_fold_change <= -1 | mean_log2_fold_change >= 1) & mean_log10_P_value >= 3 ~ "Change",
        mean_log2_fold_change > -1 & mean_log2_fold_change < 1 & mean_log10_P_value < 3 ~ "No_change",
        TRUE ~ "Indeterminate"
      ),
      direction_of_change = case_when(
        change_in_patients %in% c("No_change", "Indeterminate") ~ "-",
        mean_log2_fold_change <= -1 ~ "Decrease",
        mean_log2_fold_change >= 1 ~ "Increase"
      )
    )
  
  # Se agrega una columna para indicar el ID del paciente
  markers1_w_id <- markers1_w_id %>%
    mutate(patient_id = patient)
  
  # Se filtran los marcadores que muestran un cambio significativo en los pacientes
  changing_markers1_w_id <- markers1_w_id %>%
    filter((mean_log2_fold_change <= -1 | mean_log2_fold_change >= 1) & mean_log10_P_value >= 3)
  
  
  ## El proceso anterior se repite para los tiempos TMZ y 12M, con sus respectivas operaciones
  markers2_w_id <- merge(plot_data2,plot_data1[ ,c("sequence", "protein_identifier", "start", "mean_smoothed_signal_D65")], by = c("sequence","protein_identifier", "start"), all.x = TRUE, all.y=FALSE)
  
  markers2_w_id <- markers2_w_id %>%
    filter(mean_smoothed_signal_TMZ >= 10000 | mean_smoothed_signal_D65 >= 10000 | mean_smoothed_signal_12M >= 10000)
  
  markers2_w_id <- markers2_w_id %>%
    group_by(protein_identifier) %>%
    summarize(mean_log2_fold_change = round(mean(log2_Fold_change),2), mean_log10_P_value = round(mean(log10_P_value),2),
              mean_signal_12M = round(mean(mean_smoothed_signal_12M),2), mean_signal_TMZ = round(mean(mean_smoothed_signal_TMZ),2),
              peptide_count = n(), sd_12M = round(sd(mean_smoothed_signal_12M),2), sd_TMZ = round(sd(mean_smoothed_signal_TMZ),2), 
              sd_log2_Fold_change = round(sd(log2_Fold_change),2), sd_log10_P_value = round(sd(log10_P_value),2))
  
  markers2_w_id <- markers2_w_id %>%
    mutate(
      change_in_patients = case_when(
        (mean_log2_fold_change <= -1 | mean_log2_fold_change >= 1) & mean_log10_P_value >= 3 ~ "Change",
        mean_log2_fold_change > -1 & mean_log2_fold_change < 1 & mean_log10_P_value < 3 ~ "No_change",
        TRUE ~ "Indeterminate"
      ),
      direction_of_change = case_when(
        change_in_patients %in% c("No_change", "Indeterminate") ~ "-",
        mean_log2_fold_change <= -1 ~ "Decrease",
        mean_log2_fold_change >= 1 ~ "Increase"
      )
    )
  
  markers2_w_id <- markers2_w_id %>%
    mutate(patient_id = patient)
  
  changing_markers2_w_id <- markers2_w_id %>%
    filter((mean_log2_fold_change <= -1 | mean_log2_fold_change >= 1) & mean_log10_P_value >= 3)
  
  
  # Se guardan en un data frame los marcadores que muestran un cambio significativo en los pacientes para el primer conjunto de tiempos (D65-TMZ)
  changing_markers1_w_id_output_file <- sprintf("%s/D65-TMZ/changing_antigens_%s_%s_vs_%s.tsv", output_folder1, patient, unique_times1[1],unique_times1[2])
  write.table(changing_markers1_w_id, file = changing_markers1_w_id_output_file, sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)
  
  # Se guardan en un data frame todos los marcadores para el primer conjunto de tiempos (D65-TMZ)
  all_markers1_w_id_output_file <- sprintf("%s/D65-TMZ/all_antigens_%s_%s_vs_%s.tsv", output_folder2, patient, unique_times1[1],unique_times1[2])
  write.table(markers1_w_id, file = all_markers1_w_id_output_file, sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)
  
  # Se guardan en un data frame los marcadores que muestran un cambio significativo en los pacientes para el segundo conjunto de tiempos (12M-TMZ)
  changing_markers2_w_id_output_file <- sprintf("%s/12M-TMZ/changing_antigens_%s_%s_vs_%s.tsv", output_folder1, patient, unique_times2[1],unique_times2[2])
  write.table(changing_markers2_w_id, file = changing_markers2_w_id_output_file, sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)
  
  # Se guardan en un data frame todos los marcadores para el primer conjunto de tiempos (12M-TMZ)
  all_markers2_w_id_output_file <- sprintf("%s/12M-TMZ/all_antigens_%s_%s_vs_%s.tsv", output_folder2, patient, unique_times2[1],unique_times2[2])
  write.table(markers2_w_id, file = all_markers2_w_id_output_file, sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)
}
