# En estos gráficos se representó en el eje x el nivel de cambio en la señal (fold change) de las réplicas y en el eje y 
# la significancia de este cambio (dada por el p-valor de la prueba t de student entre ambas réplicas)
# Esto puede realizarse para los diferentes pacientes del ensayo clínico, para la comparación entre los tiempos TMZ-D65 o TMZ-12M.
# # Se colorea los puntos según la densidad de puntos en esa área. Se marcan en rojo los péptidos del compartimento variable.

# Se cargan los paquetes a utilizar
library(data.table)
library(dplyr)

# Se establece la carpeta donde están los datos de las señales filtradas
data_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/raw_data_mean_replicas_filtered_w_formula"
setwd(data_folder)

# Se establece la carpeta donde se guardará el resultado
output_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/volcano_plots"

# Se establece el archivo donde están definidas funciones que se utilizarán más adelante
source("C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R/Antigenos analizados/Function_density_plots.R")

# Se establece el sufijo de los archivos
suffix <- "raw_data_filtered.tsv"

# Se especifica el paciente o listas de pacientes que se van a utilizar
# patients<- "2-159-304"

patients<- c("1-005-101","1-049-116","1-124-150","1-160-163","1-227-184","1-264-203","2-098-277","2-132-293",
             "2-190-322","2-212-332","2-240-353","2-261-364", "1-052-118","1-105-143","1-143-157",
             "1-218-183","1-248-193","2-011-253","2-103-280","2-125-286", "2-159-304","2-177-313","2-203-331",
             "2-242-355","1-044-115","1-063-121","1-064-124","1-110-146","1-130-153","1-150-159","1-153-161","1-225-187",
             "1-245-195","2-032-256","2-211-335","2-236-352")


##############-
#### MAIN ####
##############-
# Se establecen las configuraciones de los gráficos
fixed_scale <- 0
# min_fixed_x_scale <- -6
# max_fixed_x_scale <- 5
# min_fixed_y_scale <- 0
# max_fixed_y_scale <- 5


# Bucle principal que itera sobre cada paciente
for (patient in patients) {
  
  # Se determinan los archivos donde se encuentran los datos filtrados
  raw_data_file_TMZ <-  sprintf("%s/%s_TMZ_%s", data_folder, patient, suffix)
  raw_data_file_D65 <-  sprintf("%s/%s_D65_%s", data_folder, patient, suffix)
  raw_data_file_12M <-  sprintf("%s/%s_12M_%s", data_folder, patient, suffix)
  
  # Se cargan los datos filtrados
  raw_data_TMZ <- fread(raw_data_file_TMZ, header = T, sep = "\t", na.strings = NULL)
  raw_data_D65 <- fread(raw_data_file_D65, header = T, sep = "\t", na.strings = NULL)
  raw_data_12M <- fread(raw_data_file_12M, header = T, sep = "\t", na.strings = NULL)
  
  # Se unen los datos de los tres tiempos en un único conjunto de datos
  raw_data <- rbind(raw_data_TMZ,raw_data_D65,raw_data_12M)
  
  # Se ordenan los datos
  raw_data <- raw_data[order(source, time, sequence)]
  
  # Se Filtran los datos para los tiempos necesarios
  raw_data1 <- filter(raw_data, time == "TMZ" | time == "D65")
  raw_data2 <- filter(raw_data, time == "TMZ" | time == "12M")
  
  
  ## Se obtiene el volcano plot comparando entre los tiempos "TMZ" y "D65"
  # Se sacan todos los datos que no tengan 2 tiempos
  raw_data1 <- raw_data1 %>%
    group_by(sequence, source) %>%
    mutate(time_count = n_distinct(time)) %>%
    filter(time_count == 2) %>%
    ungroup() %>%
    select(-time_count)
  
  sub_raw_data1 <- as.data.table(raw_data1)
  
  unique_times1 <- unique(sub_raw_data1$time)
  
  # Se cargan los datos necesarios para el gráfico en un data table
  plot_data1 <- data.table(sequence = unique(sub_raw_data1$sequence),
                           x1 = sub_raw_data1[time == unique_times1[1]]$replica_1_signal,
                           x2 = sub_raw_data1[time == unique_times1[1]]$replica_2_signal,
                           y1 = sub_raw_data1[time == unique_times1[2]]$replica_1_signal,
                           y2 = sub_raw_data1[time == unique_times1[2]]$replica_2_signal)
  # x=time, y=replica
  # x1 = D65, replica 1
  # x2 = D65, replica 2
  # y1 = TMZ, replica 1
  # y2 = TMZ, replica 2
  
  plot_data1 <- na.omit(plot_data1)
  
  # Se calcula el p-valor para cada péptido (mediante una prueba t de student)
  plot_data1 <- plot_data1 %>%
    rowwise() %>%
    mutate(p_value = t.test(c(x1, x2), c(y1, y2), paired = F)$p.value) %>%
    ungroup()
  
  # Se calcula el valor promedio de cada péptido
  plot_data1 <- plot_data1 %>%
    rowwise() %>%
    mutate(average_x = mean(c(x1, x2)),
           average_y = mean(c(y1, y2))) %>%
    ungroup()
  
  plot_data1 <- as.data.table(plot_data1)
  
  # Se calcula el log2(fold change) y el -log10(p-value) de cada péptido
  plot_data1$log2_Fold_change <- log2(plot_data1$average_x) - log2(plot_data1$average_y)
  plot_data1$log10_P_value <- -log10(plot_data1$p_value)
  
  # Se seleccionan las columnas necesarias
  plot_data1 <- plot_data1 %>%
    select(log2_Fold_change, log10_P_value)
  
  plot_data1 <- as.data.table(plot_data1)
  
  # Se define el título del gráfico 
  plot_title1 <- sprintf("Volcano plot | %s vs %s | %s" , unique_times1[1], unique_times1[2], patient)
  
  # Se genera el volcano plot comparando entre los tiempos "TMZ" y "D65", se establecen los umbrales para cada compartimento
  p <- plotManuallyColoredDensityScatterplot(plot_data = plot_data1, x_column_name = "log2_Fold_change", y_column_name = "log10_P_value", 
                                             bin_amount_per_axis = 1000, max_data_density_to_consider_per_cell = 100000, gradient_steps = 200,
                                             geom_point_size = 1,
                                             plot_title = plot_title1, x_label = "log2(fold change)", y_label = "-log10(p-value)", 
                                             axis_title_size = 16, axis_text_size = 16, axis_ticks_length = 0.2, panel_background_size = 1,
                                             fixed_scale = 0, 
                                             # fixed_scale_x_min = min_fixed_x_scale, fixed_scale_x_max = max_fixed_x_scale,
                                             # fixed_scale_y_min = min_fixed_y_scale, fixed_scale_y_max = max_fixed_y_scale,
                                             vt1 = -1, vt2 = 1, ht = 3)
  # Se guarda la imagen en formato png
  output_file <- sprintf("%s/Volcano_plot_%s_%s_vs_%s.png", output_folder, patient, unique_times1[1],unique_times1[2])
  png(file = output_file)
  print(p)
  dev.off()
  
  
  ## Se obtiene el volcano plot comparando entre los tiempos "TMZ" y "12M"
  ## El proceso anterior se repite para los tiempos "TMZ" y "12M", con sus respectivos gráficos
  
  raw_data2 <- raw_data2 %>%
    group_by(sequence, source) %>%
    mutate(time_count = n_distinct(time)) %>%
    filter(time_count == 2) %>%
    ungroup() %>%
    select(-time_count)
  
  sub_raw_data2 <- as.data.table(raw_data2)
  
  unique_times2 <- unique(sub_raw_data2$time)
  
  plot_data2 <- data.table(sequence = unique(sub_raw_data2$sequence),
                           x1 = sub_raw_data2[time == unique_times2[1]]$replica_1_signal,
                           x2 = sub_raw_data2[time == unique_times2[1]]$replica_2_signal,
                           y1 = sub_raw_data2[time == unique_times2[2]]$replica_1_signal,
                           y2 = sub_raw_data2[time == unique_times2[2]]$replica_2_signal)
  # x=time, y=replica
  # x1 = 12M, replica 1
  # x2 = 12M, replica 2
  # y1 = TMZ, replica 1
  # y2 = TMZ, replica 2
  
  plot_data2 <- na.omit(plot_data2)
  
  plot_data2 <- plot_data2 %>%
    rowwise() %>%
    mutate(p_value = t.test(c(x1, x2), c(y1, y2), paired = F)$p.value) %>%
    ungroup()
  
  plot_data2 <- plot_data2 %>%
    rowwise() %>%
    mutate(average_x = mean(c(x1, x2)),
           average_y = mean(c(y1, y2))) %>%
    ungroup()
  
  plot_data2 <- as.data.table(plot_data2)
  
  plot_data2$log2_Fold_change <- log2(plot_data2$average_x) - log2(plot_data2$average_y)
  plot_data2$log10_P_value <- -log10(plot_data2$p_value)
  
  plot_data2 <- plot_data2 %>%
    select(log2_Fold_change, log10_P_value)
  
  plot_data2 <- as.data.table(plot_data2)
  
  plot_title2 <- sprintf("Volcano plot | %s vs %s | %s" , unique_times2[1], unique_times2[2], patient)
  
  q <- plotManuallyColoredDensityScatterplot(plot_data = plot_data2, x_column_name = "log2_Fold_change", y_column_name = "log10_P_value", 
                                             bin_amount_per_axis = 1000, max_data_density_to_consider_per_cell = 100000, gradient_steps = 200,
                                             geom_point_size = 1,
                                             plot_title = plot_title2, x_label = "log2(fold change)", y_label = "-log10(p-value)", 
                                             axis_title_size = 16, axis_text_size = 16, axis_ticks_length = 0.2, panel_background_size = 1,
                                             fixed_scale = 0, 
                                             # fixed_scale_x_min = min_fixed_x_scale, fixed_scale_x_max = max_fixed_x_scale,
                                             # fixed_scale_y_min = min_fixed_y_scale, fixed_scale_y_max = max_fixed_y_scale,
                                             vt1 = -1, vt2 = 1, ht = 3)
  
  output_file <- sprintf("%s/Volcano_plot_%s_%s_vs_%s.png", output_folder, patient, unique_times2[1],unique_times2[2])
  png(file = output_file)
  print(q)
  dev.off()
  
}
