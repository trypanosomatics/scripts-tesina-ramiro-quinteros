# Se compararon las señales de fluorescencia provenientes de todos los péptidos del tiempo D65 o 12M (eje x) 
# contra los del tiempo TMZ (eje y).
# Se colorea los puntos según la densidad de puntos en esa área.
# Se calculan los coeficientes de correlación de pearson para estos tiempos
# Esto puede realizarse para los diferentes pacientes del ensayo clínico.

# Se cargan los paquetes a utilizar
library(gridExtra)
library(data.table)
library(ggplot2)
library(dplyr)

# Se establece la carpeta donde están los datos de las señales filtradas
data_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/raw_data_mean_replicas_filtered_w_formula"
setwd(data_folder)

# Se establece la carpeta donde se guardará el resultado
output_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/compare_filtered_times"

# Se establece el archivo donde están definidas las funciones plotManuallyColoredDensityScatterplot y addQuarterCircleToPlot,
# que se utilizarán más tarde
source("C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/Regiones ale/plot_functions_scatterplots.R")

# Se establece el sufijo de los archivos
suffix <- "raw_data_filtered.tsv"

# Se establece el nombre del tiempo que luego va a aparecer en el gráfico
name1 <- "TMZ vs D65"
name2 <- "TMZ vs 12M"

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

# Data frame para almacenar los coeficientes de pearson
pearson_data <- data.table(pearson = as.numeric(), pearson_above_cutoff = as.numeric())

# Bucle principal que itera sobre cada paciente
for (patient in patients) {
  
  # Se determinan los archivos donde se encuentran los datos filtrados
  raw_data_file_TMZ <-  sprintf("%s/%s_TMZ_%s",data_folder, patient, suffix)
  raw_data_file_D65 <-  sprintf("%s/%s_D65_%s",data_folder, patient, suffix)
  raw_data_file_12M <-  sprintf("%s/%s_12M_%s",data_folder, patient, suffix)
  
  # Se cargan los datos filtrados
  raw_data_TMZ <- fread(raw_data_file_TMZ, header = T, sep = "\t", na.strings = NULL)
  raw_data_D65 <- fread(raw_data_file_D65, header = T, sep = "\t", na.strings = NULL)
  raw_data_12M <- fread(raw_data_file_12M, header = T, sep = "\t", na.strings = NULL)
  
  # Se unen los datos de los tres tiempos en un único conjunto de datos
  raw_data <- rbind(raw_data_TMZ,raw_data_D65,raw_data_12M)
  
  # Se Seleccionan las columnas necesarias 
  raw_data <- raw_data %>%
    select(source, time, sequence, mean_signal)
  
  # Se Filtran los datos para los tiempos necesarios
  raw_data1 <- filter(raw_data, time == "TMZ" | time == "D65")
  raw_data2 <- filter(raw_data, time == "TMZ" | time == "12M")
  
  
  ## Se obtiene el gráfico (scatterplot) comparando entre los tiempos "TMZ" y "D65"

  # Se sacan todos los datos que no tengan 2 tiempos
  raw_data1 <- raw_data1 %>%
    group_by(sequence, source) %>%
    mutate(time_count = n_distinct(time)) %>%
    filter(time_count == 2) %>%
    ungroup() %>%
    select(-time_count)
  
  sub_raw_data1 <- as.data.table(raw_data1)                  
  
  unique_times1 <- c("D65","TMZ")
  
  # Se cargan los datos necesarios para el gráfico en un data table
  if (length(unique_times1) == 2) {                  
    plot_data <- data.table(x = sub_raw_data1[time == unique_times1[1]]$mean_signal,
                            y = sub_raw_data1[time == unique_times1[2]]$mean_signal)
    
    # Se calcula el coeficiente de correlación de Pearson entre los dos tiempos
    pearson_aux <- cor(plot_data$x, plot_data$y, method = "pearson")
    
    # Se calcula el punto de corte para los valores superiores en el coeficiente de Pearson
    plot_data[, distance := (x^2 + y^2)^(1/2)]
    quantile_cutoff_distance <- quantile(plot_data$distance, quantile_cutoff)
    quantile_cutoff_plot_data <- plot_data[distance >= quantile_cutoff_distance]
    
    # Se calcula el coeficiente de correlación de Pearson para los datos sobre el punto de corte
    quantile_cutoff_pearson_aux <- cor(quantile_cutoff_plot_data$x, quantile_cutoff_plot_data$y, method = "pearson")
    
    # Se define el título del gráfico y el de los ejes 
    plot_title <- sprintf("%s vs %s | %s | Pearson = %s / %s", unique_times1[1], unique_times1[2], patient, round(pearson_aux, 2), round(quantile_cutoff_pearson_aux, 2))
    
    plot_x_label <- paste(unique_times1[1], "signal", sep = " ")
    plot_y_label <- paste(unique_times1[2], "signal", sep = " ")
    
    # Se genera el gráfico de dispersión comparando entre los tiempos "TMZ" y "D65"
    p <- plotManuallyColoredDensityScatterplot(plot_data = plot_data, x_column_name = "x", y_column_name = "y", 
                                               bin_amount_per_axis = 1000, max_data_density_to_consider_per_cell = 100000, gradient_steps = 200,
                                               geom_point_size = 1,
                                               plot_title = plot_title, x_label = plot_x_label, y_label = plot_y_label, 
                                               axis_title_size = 16, axis_text_size = 16, axis_ticks_length = 0.2, panel_background_size = 1,
                                               fixed_scale = 1, fixed_scale_min = min_fixed_scale, fixed_scale_max = max_fixed_scale)
    
    # Se agrega el círculo del cuantil al gráfico
    a <- addQuarterCircleToPlot(p = p, radius = quantile_cutoff_distance,
                                geom_line_color = "orange", geom_line_size = 1.5, geom_line_linetype = "22",
                                origin_x = 0, origin_y = 0)
    
  }
  
  ## Se obtiene el gráfico (scatterplot) comparando entre los tiempos "TMZ" y "D65"
  ## El proceso anterior se repite para los tiempos "TMZ" y "12M", con sus respectivos gráficos y cálculos de Pearson
  
  raw_data2 <- raw_data2 %>%
    group_by(sequence, source) %>%
    mutate(time_count = n_distinct(time)) %>%
    filter(time_count == 2) %>%
    ungroup() %>%
    select(-time_count)
  
  sub_raw_data2 <- as.data.table(raw_data2)
  
  unique_times2 <- c("12M","TMZ")
  
  if (length(unique_times2) == 2) {
    plot_data2 <- data.table(x2 = sub_raw_data2[time == unique_times2[1]]$mean_signal,
                             y2 = sub_raw_data2[time == unique_times2[2]]$mean_signal)
    
    pearson_aux2 <- cor(plot_data2$x2, plot_data2$y2, method = "pearson")
    
    plot_data2[, distance := (x2^2 + y2^2)^(1/2)]
    quantile_cutoff_distance2 <- quantile(plot_data2$distance, quantile_cutoff)
    quantile_cutoff_plot_data2 <- plot_data2[distance >= quantile_cutoff_distance2]
    
    quantile_cutoff_pearson_aux2 <- cor(quantile_cutoff_plot_data2$x2, quantile_cutoff_plot_data2$y2, method = "pearson")
    
    plot_title2 <- sprintf("%s vs %s | %s | Pearson = %s / %s", unique_times2[1], unique_times2[2], patient, round(pearson_aux2, 2), round(quantile_cutoff_pearson_aux2, 2))
    
    plot_x_label2 <- paste(unique_times2[1], "signal", sep = " ")
    plot_y_label2 <- paste(unique_times2[2], "signal", sep = " ")
    
    q <- plotManuallyColoredDensityScatterplot(plot_data = plot_data2, x_column_name = "x2", y_column_name = "y2", 
                                               bin_amount_per_axis = 1000, max_data_density_to_consider_per_cell = 100000, gradient_steps = 200,
                                               geom_point_size = 1,
                                               plot_title = plot_title2, x_label = plot_x_label2, y_label = plot_y_label2, 
                                               axis_title_size = 16, axis_text_size = 16, axis_ticks_length = 0.2, panel_background_size = 1,
                                               fixed_scale = 1, fixed_scale_min = min_fixed_scale, fixed_scale_max = max_fixed_scale)
    
    b <- addQuarterCircleToPlot(p = q, radius = quantile_cutoff_distance2,
                                geom_line_color = "orange", geom_line_size = 1.5, geom_line_linetype = "22",
                                origin_x = 0, origin_y = 0)
    
    # Se combinan los gráficos
    combined_plots <- list(a, b)
    
    
    # Se guarda la imagen en formato PDF
    output_pdf <- sprintf("%s/Compare_times_with_circle_%s.pdf", output_folder, patient)
    pdf(file = output_pdf, width = 13, height = 6)
    do.call(grid.arrange, c(combined_plots, ncol = 2))
    dev.off()
    
    # Se guarda la imagen en formato TIFF
    # output_tiff <- sprintf("Compare_times_with_circle_%s_%s.tiff", patient, normalization)
    # tiff(file = output_tiff)
    # do.call(grid.arrange, c(combined_plots, ncol = 2))
    # dev.off()
  }
}