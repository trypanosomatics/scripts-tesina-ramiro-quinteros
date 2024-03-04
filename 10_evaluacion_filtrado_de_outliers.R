# Se realiza un control de calidad de las réplicas técnicas. Para esto se hace un gráfico recíproco comparando 
# las señales de fluorescencia provenientes de todos los péptidos de la réplica 1 (eje x) contra los de la réplica 2 (eje y). 
# Se colorea los puntos según la densidad de puntos en esa área. Se marcan en rojo los puntos eliminados en el proceso de filtrado de los outliers
# Se calculan los coeficientes de correlación de pearson para estas réplicas (lo ideal es 1).
# Esto puede realizarse para los diferentes pacientes del ensayo clínico, para los 3 tiempos de cada paciente.


# Se cargan los paquetes a utilizar
library(gridExtra)
library(data.table)
library(ggplot2)
library(dplyr)

# Se establece la carpeta donde están los datos de las señales crudas, con puntaje de variación
data_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/raw_data_mean_replicas_w_variation"

# Se establece la carpeta donde se guardará el resultado
output_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/compare_filtered_replicas_eliminated_peptides"

# Se establece el archivo donde están definidas las funciones plotManuallyColoredDensityScatterplot y addQuarterCircleToPlot,
# que se utilizarán más tarde
source("C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/Regiones ale/plot_functions_scatterplots.R")

# Se establece el sufijo de los archivos
suffix <- "raw_data.tsv"


# Se establece el nombre del tiempo que luego va a aparecer en el gráfico
name_TMZ <- "TMZ"
name_D65 <- "D65"
name_12M <- "12M"

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
  raw_data_file_TMZ <-  sprintf("%s/%s_TMZ_%s", data_folder, patient, suffix)
  raw_data_file_D65 <-  sprintf("%s/%s_D65_%s", data_folder, patient, suffix)
  raw_data_file_12M <-  sprintf("%s/%s_12M_%s", data_folder, patient, suffix)
  
  # Se cargan los datos sin procesar
  raw_data_TMZ <- fread(raw_data_file_TMZ, header = T, sep = "\t", na.strings = NULL)
  raw_data_D65 <- fread(raw_data_file_D65, header = T, sep = "\t", na.strings = NULL)
  raw_data_12M <- fread(raw_data_file_12M, header = T, sep = "\t", na.strings = NULL)
  
  # Se Seleccionan las columnas necesarias 
  raw_data_TMZ <- raw_data_TMZ %>%
    select(source, time, replica_1_signal, replica_2_signal, sequence, variation_score, relative_variation, absolute_variation)
  
  raw_data_D65 <- raw_data_D65 %>%
    select(source, time, replica_1_signal, replica_2_signal, sequence, variation_score, relative_variation, absolute_variation)
  
  raw_data_12M <- raw_data_12M %>%
    select(source, time, replica_1_signal, replica_2_signal, sequence, variation_score, relative_variation, absolute_variation)
  
  # Se consiguen los datos de los péptidos eliminados para cada tiempo
  eliminated_data_TMZ <- raw_data_TMZ %>%
    filter(!(variation_score < 0.16 | relative_variation < 1.5 | absolute_variation < 1000))
  
  eliminated_data_D65 <- raw_data_D65 %>%
    filter(!(variation_score < 0.16 | relative_variation < 1.5 | absolute_variation < 1000))
  
  eliminated_data_12M <- raw_data_12M %>%
    filter(!(variation_score < 0.16 | relative_variation < 1.5 | absolute_variation < 1000))
  
  eliminated_plot_data_TMZ <- as.data.table(eliminated_data_TMZ)
  eliminated_plot_data_D65 <- as.data.table(eliminated_data_D65)
  eliminated_plot_data_12M <- as.data.table(eliminated_data_12M)
  
  # Se consiguen los datos de los péptidos que se mantienen en el análisis para cada tiempo
  raw_data_TMZ <- raw_data_TMZ %>%
    filter(variation_score < 0.16 | relative_variation < 1.5 | absolute_variation < 1000)
  
  raw_data_D65 <- raw_data_D65 %>%
    filter(variation_score < 0.16 | relative_variation < 1.5 | absolute_variation < 1000)
  
  raw_data_12M <- raw_data_12M %>%
    filter(variation_score < 0.16 | relative_variation < 1.5 | absolute_variation < 1000)
  
  
  ## Se obtiene el gráfico (reciprocal plot) comparando entre ambas réplicas del tiempo TMZ, marcando en color
  ## rojo los puntos eliminados en el filtrado de outliers
  
  sub_raw_data_TMZ <- as.data.table(raw_data_TMZ)
  # Se cargan los datos necesarios para el gráfico en un data table
  plot_data_TMZ <- data.table(x = sub_raw_data_TMZ$replica_1_signal,
                              y = sub_raw_data_TMZ$replica_2_signal)
  
  # Se calcula el coeficiente de correlación de Pearson entre las dos réplicas
  pearson_aux_TMZ <- cor(plot_data_TMZ$x, plot_data_TMZ$y, method = "pearson")
  
  # Se calcula el punto de corte para los valores superiores en el coeficiente de Pearson
  plot_data_TMZ[, distance := (x^2 + y^2)^(1/2)]
  quantile_cutoff_distance_TMZ <- quantile(plot_data_TMZ$distance, quantile_cutoff)
  quantile_cutoff_plot_data_TMZ <- plot_data_TMZ[distance >= quantile_cutoff_distance_TMZ]
  
  # Se calcula el coeficiente de correlación de Pearson para los datos sobre el punto de corte
  quantile_cutoff_pearson_aux_TMZ <- cor(quantile_cutoff_plot_data_TMZ$x, quantile_cutoff_plot_data_TMZ$y, method = "pearson")
  
  # Se define el título del gráfico y el de los ejes para el tiempo TMZ
  plot_title_TMZ <- sprintf("%s | %s | Pearson = %s / %s", name_TMZ, patient, round(pearson_aux_TMZ, 2), round(quantile_cutoff_pearson_aux_TMZ, 2))
  
  plot_x_label_TMZ <- "Replica 1 signal"
  plot_y_label_TMZ <- "Replica 2 signal"
  
  # Se genera el gráfico de dispersión con densidad para el tiempo TMZ (señales no eliminadas)
  p <- plotManuallyColoredDensityScatterplot(plot_data = plot_data_TMZ, x_column_name = "x", y_column_name = "y", 
                                             bin_amount_per_axis = 1000, max_data_density_to_consider_per_cell = 100000, gradient_steps = 200,
                                             geom_point_size = 1,
                                             plot_title = plot_title_TMZ, x_label = plot_x_label_TMZ, y_label = plot_y_label_TMZ, 
                                             axis_title_size = 16, axis_text_size = 16, axis_ticks_length = 0.2, panel_background_size = 1,
                                             fixed_scale = 1, fixed_scale_min = min_fixed_scale, fixed_scale_max = max_fixed_scale)
  
  # Se agrega el círculo del cuantil al gráfico
  a <- addQuarterCircleToPlot(p = p, radius = quantile_cutoff_distance_TMZ,
                              geom_line_color = "orange", geom_line_size = 1.5, geom_line_linetype = "22",
                              origin_x = 0, origin_y = 0)
  
  # Se agregan las señales eliminadas, y son marcadas en color rojo
  a <- a + geom_point(data = eliminated_plot_data_TMZ, aes(x = replica_1_signal, y = replica_2_signal), color = "#b00b0b")
  
  
  
  ## Se obtiene el gráfico (reciprocal plot) comparando entre ambas réplicas del tiempo D65, marcando en color
  ## rojo los puntos eliminados en el filtrado de outliers
  ## El proceso anterior se repite para el tiempo "D65", con sus respectivos gráficos y cálculos de Pearson
  
  sub_raw_data_D65 <- as.data.table(raw_data_D65)
  
  plot_data_D65 <- data.table(x = sub_raw_data_D65$replica_1_signal,
                              y = sub_raw_data_D65$replica_2_signal)
  
  pearson_aux_D65 <- cor(plot_data_D65$x, plot_data_D65$y, method = "pearson")
  
  plot_data_D65[, distance := (x^2 + y^2)^(1/2)]
  quantile_cutoff_distance_D65 <- quantile(plot_data_D65$distance, quantile_cutoff)
  quantile_cutoff_plot_data_D65 <- plot_data_D65[distance >= quantile_cutoff_distance_D65]
  
  quantile_cutoff_pearson_aux_D65 <- cor(quantile_cutoff_plot_data_D65$x, quantile_cutoff_plot_data_D65$y, method = "pearson")
  
  plot_title_D65 <- sprintf("%s | %s | Pearson = %s / %s", name_D65, patient, round(pearson_aux_D65, 2), round(quantile_cutoff_pearson_aux_D65, 2))
  
  plot_x_label_D65 <- "Replica 1 signal"
  plot_y_label_D65 <- "Replica 2 signal"
  
  q <- plotManuallyColoredDensityScatterplot(plot_data = plot_data_D65, x_column_name = "x", y_column_name = "y", 
                                             bin_amount_per_axis = 1000, max_data_density_to_consider_per_cell = 100000, gradient_steps = 200,
                                             geom_point_size = 1,
                                             plot_title = plot_title_D65, x_label = plot_x_label_D65, y_label = plot_y_label_D65, 
                                             axis_title_size = 16, axis_text_size = 16, axis_ticks_length = 0.2, panel_background_size = 1,
                                             fixed_scale = 1, fixed_scale_min = min_fixed_scale, fixed_scale_max = max_fixed_scale)
  
  b <- addQuarterCircleToPlot(p = q, radius = quantile_cutoff_distance_D65,
                              geom_line_color = "orange", geom_line_size = 1.5, geom_line_linetype = "22",
                              origin_x = 0, origin_y = 0)
  
  b <- b + geom_point(data = eliminated_plot_data_D65, aes(x = replica_1_signal, y = replica_2_signal), color = "#b00b0b")
  
  
  ## Se obtiene el gráfico (reciprocal plot) comparando entre ambas réplicas del tiempo 12M, marcando en color
  ## rojo los puntos eliminados en el filtrado de outliers
  ## El proceso anterior se repite para el tiempo "12M", con sus respectivos gráficos y cálculos de Pearson
  
  sub_raw_data_12M <- as.data.table(raw_data_12M)
  
  plot_data_12M <- data.table(x = sub_raw_data_12M$replica_1_signal,
                              y = sub_raw_data_12M$replica_2_signal)
  
  pearson_aux_12M <- cor(plot_data_12M$x, plot_data_12M$y, method = "pearson")
  
  plot_data_12M[, distance := (x^2 + y^2)^(1/2)]
  quantile_cutoff_distance_12M <- quantile(plot_data_12M$distance, quantile_cutoff)
  quantile_cutoff_plot_data_12M <- plot_data_12M[distance >= quantile_cutoff_distance_12M]
  
  quantile_cutoff_pearson_aux_12M <- cor(quantile_cutoff_plot_data_12M$x, quantile_cutoff_plot_data_12M$y, method = "pearson")
  
  plot_title_12M <- sprintf("%s | %s | Pearson = %s / %s", name_12M, patient, round(pearson_aux_12M, 2), round(quantile_cutoff_pearson_aux_12M, 2))
  
  plot_x_label_12M <- "Replica 1 signal"
  plot_y_label_12M <- "Replica 2 signal"
  
  r <- plotManuallyColoredDensityScatterplot(plot_data = plot_data_12M, x_column_name = "x", y_column_name = "y", 
                                             bin_amount_per_axis = 1000, max_data_density_to_consider_per_cell = 100000, gradient_steps = 200,
                                             geom_point_size = 1,
                                             plot_title = plot_title_12M, x_label = plot_x_label_12M, y_label = plot_y_label_12M, 
                                             axis_title_size = 16, axis_text_size = 16, axis_ticks_length = 0.2, panel_background_size = 1,
                                             fixed_scale = 1, fixed_scale_min = min_fixed_scale, fixed_scale_max = max_fixed_scale)
  
  c <- addQuarterCircleToPlot(p = r, radius = quantile_cutoff_distance_12M,
                              geom_line_color = "orange", geom_line_size = 1.5, geom_line_linetype = "22",
                              origin_x = 0, origin_y = 0)
  
  c <- c + geom_point(data = eliminated_plot_data_12M, aes(x = replica_1_signal, y = replica_2_signal), color = "#b00b0b")
  
  
  # Se combinan los gráficos
  combined_plots <- list(a,b,c)
  
  
  # Se guarda la imagen en formato PDF
  output_pdf <- sprintf("%s/Compare_filtered_replicas_all_times_eliminated_peptides_%s.pdf", output_folder, patient)
  pdf(file = output_pdf, width = 20, height = 7)
  do.call(grid.arrange, c(combined_plots, ncol = 3))
  dev.off()
}