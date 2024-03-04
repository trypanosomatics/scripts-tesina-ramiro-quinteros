# Se realiza un control de calidad de las réplicas técnicas. Para esto se hace un gráfico recíproco comparando 
# las señales de fluorescencia provenientes de todos los péptidos de la réplica 1 (eje x) contra los de la réplica 2 (eje y). 
# Se colorea los puntos según la densidad de puntos en esa área.
# Se calculan los coeficientes de correlación de pearson para estas réplicas (lo ideal es 1).
# Esto puede realizarse para los diferentes pacientes del ensayo clínico, para los 3 tiempos de cada paciente.

# Para evaluar el filtrado de outliers se debe utilizar la carpeta (variable data_folder) donde se encuentran los datos
# filtrados, además de usar el sufijo (variable suffix) que corresponde para estos archivos. Además se pueden modificar 
# el nombre de la carpeta de salida (variable output_folder) o el nombre de salida (variable output_pdf).

# Se cargan los paquetes a utilizar
library(gridExtra)
library(data.table)
library(ggplot2)
library(dplyr)

# Se establece la carpeta donde están los datos de las señales crudas
data_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/Array Roche E1224/00_inputs/02_serums/raw_data"
setwd(data_folder)

# Se establece la carpeta donde están los datos del diseño del experimento
data_folder2 <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/design"

# Se establece la carpeta donde se guardará el resultado
output_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/raw_data_plots/compare_replicas_raw_data"

# Se establece el archivo donde están definidas las funciones plotManuallyColoredDensityScatterplot y addQuarterCircleToPlot,
# que se utilizarán más tarde
source("C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/Regiones ale/plot_functions_scatterplots.R")

# Se establece el sufijo de los archivos
suffix <- "raw_data"

# Se establece el nombre del tiempo que luego va a aparecer en el gráfico
name_TMZ <- "TMZ"
name_D65 <- "D65"
name_12M <- "12M"


t_cruzi_map_data_file <- sprintf("%s/01_individual_serum_design_v3.tsv", data_folder2)
t_cruzi_map <- fread(t_cruzi_map_data_file, header = T, sep = "\t", na.strings = NULL)

# Se determinan los orígenes y suborígenes a usar
#  origins_to_use <- unique(t_cruzi_map$origin) #para el análisis de todos los orígenes se utiliza este
#  suborigins_to_use <- unique(t_cruzi_map$suborigin) #para el análisis de todos los suborígenes se utiliza este

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
fixed_scale <- 1
min_fixed_scale <- 0
max_fixed_scale <- 60000
quantile_cutoff <- 0.99 #for top values pearson

# Bucle principal que itera sobre cada paciente
for (patient in patients) {
  
  # Se determinan los archivos donde se encuentran los datos sin procesar
  raw_data_file_TMZ <-  sprintf("%s/%s_TMZ_%s",data_folder, patient, suffix)
  raw_data_file_D65 <-  sprintf("%s/%s_D65_%s",data_folder, patient, suffix)
  raw_data_file_12M <-  sprintf("%s/%s_12M_%s",data_folder, patient, suffix)
  
  # Se cargan los datos sin procesar
  raw_data_TMZ <- fread(raw_data_file_TMZ, header = T, sep = "\t", na.strings = NULL)
  raw_data_D65 <- fread(raw_data_file_D65, header = T, sep = "\t", na.strings = NULL)
  raw_data_12M <- fread(raw_data_file_12M, header = T, sep = "\t", na.strings = NULL)
  
  # Se unen los datos de los tres tiempos en un único conjunto de datos
  raw_data <- rbind(raw_data_TMZ,raw_data_D65,raw_data_12M) %>%
    rename(time = type) # Se renombra la columna "type" como "time" para mayor claridad
  
  # Se filtran los datos de acuerdo a los péptidos seleccionados
  raw_data <- raw_data[sequence %in% peptides_to_use]   
  
  # Se ordenan los datos
  raw_data <- raw_data[order(source, time, replica, sequence)]
  
  # Se filtran los datos para cada tiempo
  raw_data_TMZ <- filter(raw_data, time == "TMZ")
  raw_data_D65 <- filter(raw_data, time == "D65")
  raw_data_12M <- filter(raw_data, time == "12M")
  
  
  # Se utilizan los datos solo con dos réplicas, por las dudas
  raw_data_TMZ <- raw_data_TMZ %>%
    group_by(sequence, source, time) %>%
    mutate(replica_count = n_distinct(replica)) %>%
    filter(replica_count == 2) %>%
    ungroup() %>%
    select(-replica_count)
  
  sub_raw_data_TMZ <- as.data.table(raw_data_TMZ)
  
  unique_replicas_TMZ <- unique(sub_raw_data_TMZ$replica)
  # Se cargan los datos necesarios para el gráfico en un data table
  if (length(unique_replicas_TMZ) == 2) {
    plot_data_TMZ <- data.table(x_TMZ = sub_raw_data_TMZ[replica == unique_replicas_TMZ[1]]$signal,
                                y_TMZ = sub_raw_data_TMZ[replica == unique_replicas_TMZ[2]]$signal)
    
    # Se calcula el coeficiente de correlación de Pearson entre las dos réplicas
    pearson_aux_TMZ <- cor(plot_data_TMZ$x_TMZ, plot_data_TMZ$y_TMZ, method = "pearson")
    
    # Se calcula el punto de corte para los valores superiores en el coeficiente de Pearson
    plot_data_TMZ[, distance := (x_TMZ^2 + y_TMZ^2)^(1/2)]
    quantile_cutoff_distance_TMZ <- quantile(plot_data_TMZ$distance, quantile_cutoff)
    quantile_cutoff_plot_data_TMZ <- plot_data_TMZ[distance >= quantile_cutoff_distance_TMZ]
    
    # Se calcula el coeficiente de correlación de Pearson para los datos sobre el punto de corte
    quantile_cutoff_pearson_aux_TMZ <- cor(quantile_cutoff_plot_data_TMZ$x, quantile_cutoff_plot_data_TMZ$y, method = "pearson")
    
    # Se define el título del gráfico y el de los ejes para el tiempo TMZ
    plot_title_TMZ <- sprintf("%s | %s | Pearson = %s / %s", name_TMZ, patient, round(pearson_aux_TMZ, 2), round(quantile_cutoff_pearson_aux_TMZ, 2))
    
    plot_x_label_TMZ <- paste("Replica", unique_replicas_TMZ[1], "signal", sep = " ")
    plot_y_label_TMZ <- paste("Replica", unique_replicas_TMZ[2], "signal", sep = " ")
    
    # Se genera el gráfico de dispersión con densidad para el tiempo TMZ
    p <- plotManuallyColoredDensityScatterplot(plot_data = plot_data_TMZ, x_column_name = "x_TMZ", y_column_name = "y_TMZ", 
                                               bin_amount_per_axis = 1000, max_data_density_to_consider_per_cell = 100000, gradient_steps = 200,
                                               geom_point_size = 1,
                                               plot_title = plot_title_TMZ, x_label = plot_x_label_TMZ, y_label = plot_y_label_TMZ, 
                                               axis_title_size = 16, axis_text_size = 16, axis_ticks_length = 0.2, panel_background_size = 1,
                                               fixed_scale = 1, fixed_scale_min = min_fixed_scale, fixed_scale_max = max_fixed_scale)
    
    # Se agrega el círculo del cuantil al gráfico
    a <- addQuarterCircleToPlot(p = p, radius = quantile_cutoff_distance_TMZ,
                                geom_line_color = "orange", geom_line_size = 1.5, geom_line_linetype = "22",
                                origin_x = 0, origin_y = 0)
  }
  
  
  ## El proceso anterior se repite para el tiempo D65, con sus respectivos gráficos y cálculos de Pearson
  raw_data_D65 <- raw_data_D65 %>%
    group_by(sequence, source, time) %>%
    mutate(replica_count = n_distinct(replica)) %>%
    filter(replica_count == 2) %>%
    ungroup() %>%
    select(-replica_count)
  
  sub_raw_data_D65 <- as.data.table(raw_data_D65)
  
  
  unique_replicas_D65 <- unique(sub_raw_data_D65$replica)
  if (length(unique_replicas_D65) == 2) {
    plot_data_D65 <- data.table(x_D65 = sub_raw_data_D65[replica == unique_replicas_D65[1]]$signal,
                                y_D65 = sub_raw_data_D65[replica == unique_replicas_D65[2]]$signal)
    
    pearson_aux_D65 <- cor(plot_data_D65$x_D65, plot_data_D65$y_D65, method = "pearson")
    
    plot_data_D65[, distance := (x_D65^2 + y_D65^2)^(1/2)]
    quantile_cutoff_distance_D65 <- quantile(plot_data_D65$distance, quantile_cutoff)
    quantile_cutoff_plot_data_D65 <- plot_data_D65[distance >= quantile_cutoff_distance_D65]
    
    quantile_cutoff_pearson_aux_D65 <- cor(quantile_cutoff_plot_data_D65$x, quantile_cutoff_plot_data_D65$y, method = "pearson")
    
    plot_title_D65 <- sprintf("%s | %s | Pearson = %s / %s", name_D65, patient, round(pearson_aux_D65, 2), round(quantile_cutoff_pearson_aux_D65, 2))
    
    plot_x_label_D65 <- paste("Replica", unique_replicas_D65[1], "signal", sep = " ")
    plot_y_label_D65 <- paste("Replica", unique_replicas_D65[2], "signal", sep = " ")
    
    q <- plotManuallyColoredDensityScatterplot(plot_data = plot_data_D65, x_column_name = "x_D65", y_column_name = "y_D65", 
                                               bin_amount_per_axis = 1000, max_data_density_to_consider_per_cell = 100000, gradient_steps = 200,
                                               geom_point_size = 1,
                                               plot_title = plot_title_D65, x_label = plot_x_label_D65, y_label = plot_y_label_D65, 
                                               axis_title_size = 16, axis_text_size = 16, axis_ticks_length = 0.2, panel_background_size = 1,
                                               fixed_scale = 1, fixed_scale_min = min_fixed_scale, fixed_scale_max = max_fixed_scale)
    
    b <- addQuarterCircleToPlot(p = q, radius = quantile_cutoff_distance_D65,
                                geom_line_color = "orange", geom_line_size = 1.5, geom_line_linetype = "22",
                                origin_x = 0, origin_y = 0)
  }
  
  
  ## El proceso anterior se repite para el tiempo 12M, con sus respectivos gráficos y cálculos de Pearson
  raw_data_12M <- raw_data_12M %>%
    group_by(sequence, source,  time) %>%
    mutate(replica_count = n_distinct(replica)) %>%
    filter(replica_count == 2) %>%
    ungroup() %>%
    select(-replica_count)
  
  sub_raw_data_12M <- as.data.table(raw_data_12M)
  
  
  unique_replicas_12M <- unique(sub_raw_data_12M$replica)
  if (length(unique_replicas_12M) == 2) {
    plot_data_12M <- data.table(x_12M = sub_raw_data_12M[replica == unique_replicas_12M[1]]$signal,
                                y_12M = sub_raw_data_12M[replica == unique_replicas_12M[2]]$signal)
    
    pearson_aux_12M <- cor(plot_data_12M$x_12M, plot_data_12M$y_12M, method = "pearson")
    
    plot_data_12M[, distance := (x_12M^2 + y_12M^2)^(1/2)]
    quantile_cutoff_distance_12M <- quantile(plot_data_12M$distance, quantile_cutoff)
    quantile_cutoff_plot_data_12M <- plot_data_12M[distance >= quantile_cutoff_distance_12M]
    
    quantile_cutoff_pearson_aux_12M <- cor(quantile_cutoff_plot_data_12M$x, quantile_cutoff_plot_data_12M$y, method = "pearson")
    
    plot_title_12M <- sprintf("%s | %s | Pearson = %s / %s", name_12M, patient, round(pearson_aux_12M, 2), round(quantile_cutoff_pearson_aux_12M, 2))
    # plot_title_12M <- sprintf("%s | Replica %s vs Replica %s | %s | Pearson = %s / %s", name_12M , unique_replicas_12M[1], unique_replicas_12M[2], patient, round(pearson_aux_12M, 2), round(quantile_cutoff_pearson_aux_12M, 2))
    
    plot_x_label_12M <- paste("Replica", unique_replicas_12M[1], "signal", sep = " ")
    plot_y_label_12M <- paste("Replica", unique_replicas_12M[2], "signal", sep = " ")
    
    r <- plotManuallyColoredDensityScatterplot(plot_data = plot_data_12M, x_column_name = "x_12M", y_column_name = "y_12M", 
                                               bin_amount_per_axis = 1000, max_data_density_to_consider_per_cell = 100000, gradient_steps = 200,
                                               geom_point_size = 1,
                                               plot_title = plot_title_12M, x_label = plot_x_label_12M, y_label = plot_y_label_12M, 
                                               axis_title_size = 16, axis_text_size = 16, axis_ticks_length = 0.2, panel_background_size = 1,
                                               fixed_scale = 1, fixed_scale_min = min_fixed_scale, fixed_scale_max = max_fixed_scale)
    
    c <- addQuarterCircleToPlot(p = r, radius = quantile_cutoff_distance_12M,
                                geom_line_color = "orange", geom_line_size = 1.5, geom_line_linetype = "22",
                                origin_x = 0, origin_y = 0)
  }
  
  # Se combinan los gráficos en una sola imagen
  combined_plots <- list(a,b,c)
  
  
  # Se guarda la imagen en formato PDF
  output_pdf <- sprintf("%s/Compare_replicas_all_times_%s.pdf", output_folder, patient)
  pdf(file = output_pdf, width = 20, height = 7)
  do.call(grid.arrange, c(combined_plots, ncol = 3))
  dev.off()
  
  # Se guarda la imagen en formato TIFF
  # output_tiff <- sprintf("%s/Compare_replicas_all_times_%s.tiff", output_folder, patient)
  # tiff(file = output_tiff, width = 720, height = 480)
  # do.call(grid.arrange, c(combined_plots, ncol = 3))
  # dev.off()
}