# Se realiza un control de calidad de las réplicas técnicas. Para esto se crea un MA plot para las réplicas
# (utilizado en microarrays de ADN para evaluar réplicas).
# Calcula la pendiente de la regresión de estos puntos (lo ideal es 0).
# En los MA plots, cada la señal se reemplaza por su log2 y luego la señal promedio de las dos réplicas (A) 
# se traza contra la diferencia de las dos señales (M). Se colorea los puntos según la densidad de puntos en esa área.

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

# Se establece la carpeta donde se guardará el resultado
output_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/raw_data_plots/ma_plots_raw_data"

# Se carga un archivo con funciones específicas para los gráficos
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
min_fixed_scale_x <- 4
max_fixed_scale_x <- 16
min_fixed_scale_y <- -10
max_fixed_scale_y <- 10

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
  
  # Se seleccionan las columnas relevantes para el análisis
  raw_data_TMZ <- raw_data_TMZ %>%
    select(source, time, replica_1_signal, replica_2_signal, sequence)
  
  raw_data_D65 <- raw_data_D65 %>%
    select(source, time, replica_1_signal, replica_2_signal, sequence)
  
  raw_data_12M <- raw_data_12M %>%
    select(source, time, replica_1_signal, replica_2_signal, sequence)
  
  
  ## Se generan los MA plots para el tiempo "TMZ"
  
  sub_raw_data_TMZ <- as.data.table(raw_data_TMZ)
  
  # Se cargan los datos necesarios para el gráfico en un data table
  plot_data_TMZ <- data.table(x_TMZ = sub_raw_data_TMZ$replica_1_signal,
                              y_TMZ = sub_raw_data_TMZ$replica_2_signal)
  
  # Se calculan las diferencias y medias logarítmicas
  plot_data_TMZ$M <- log2(plot_data_TMZ$x_TMZ) - log2(plot_data_TMZ$y_TMZ)
  plot_data_TMZ$A <- (log2(plot_data_TMZ$x_TMZ) + log2(plot_data_TMZ$y_TMZ)) / 2
  
  # Se ajusta una regresión lineal a los datos MA
  fit_TMZ <- lm(plot_data_TMZ$M ~ plot_data_TMZ$A)
  
  # Se crea el título del gráfico
  plot_title_TMZ <- sprintf("%s Replicas | %s | Regression slope: %s", name_TMZ, patient, round(fit_TMZ$coefficients[[2]], 2))
  
  # Se genera el gráfico MA con densidad de puntos coloreada
  p <- plotManuallyColoredDensityScatterplot(plot_data = plot_data_TMZ, x_column_name = "A", y_column_name = "M", 
                                             bin_amount_per_axis = 1000, max_data_density_to_consider_per_cell = 100000, gradient_steps = 200,
                                             geom_point_size = 1,
                                             plot_title = plot_title_TMZ, x_label = "A", y_label = "M", 
                                             axis_title_size = 16, axis_text_size = 16, axis_ticks_length = 0.2, panel_background_size = 1,
                                             fixed_scale = 0)
  # Se ajusta el rango de los ejes
  p <- p + coord_cartesian(xlim = c(min_fixed_scale_x, max_fixed_scale_x), ylim = c(min_fixed_scale_y, max_fixed_scale_y))
  
  
  ## Se generan los MA plots para el tiempo "D65". Se repite el proceso anterior
  
  sub_raw_data_D65 <- as.data.table(raw_data_D65)
  
  plot_data_D65 <- data.table(x_D65 = sub_raw_data_D65$replica_1_signal,
                              y_D65 = sub_raw_data_D65$replica_2_signal)
  
  plot_data_D65$M <- log2(plot_data_D65$x_D65) - log2(plot_data_D65$y_D65)
  plot_data_D65$A <- (log2(plot_data_D65$x_D65) + log2(plot_data_D65$y_D65)) / 2
  
  fit_D65 <- lm(plot_data_D65$M ~ plot_data_D65$A)
  
  plot_title_D65 <- sprintf("%s Replicas | %s | Regression slope: %s", name_D65, patient, round(fit_D65$coefficients[[2]], 2))
  
  
  q <- plotManuallyColoredDensityScatterplot(plot_data = plot_data_D65, x_column_name = "A", y_column_name = "M", 
                                             bin_amount_per_axis = 1000, max_data_density_to_consider_per_cell = 100000, gradient_steps = 200,
                                             geom_point_size = 1,
                                             plot_title = plot_title_D65, x_label = "A", y_label = "M", 
                                             axis_title_size = 16, axis_text_size = 16, axis_ticks_length = 0.2, panel_background_size = 1,
                                             fixed_scale = 0)
  
  q <- q + coord_cartesian(xlim = c(min_fixed_scale_x, max_fixed_scale_x), ylim = c(min_fixed_scale_y, max_fixed_scale_y))
  
  
  ## Se generan los MA plots para el tiempo "12M". Se repite el proceso anterior
  
  sub_raw_data_12M <- as.data.table(raw_data_12M)
  
  plot_data_12M <- data.table(x_12M = sub_raw_data_12M$replica_1_signal,
                              y_12M = sub_raw_data_12M$replica_2_signal)
  
  plot_data_12M$M <- log2(plot_data_12M$x_12M) - log2(plot_data_12M$y_12M)
  plot_data_12M$A <- (log2(plot_data_12M$x_12M) + log2(plot_data_12M$y_12M)) / 2
  
  fit_12M <- lm(plot_data_12M$M ~ plot_data_12M$A)
  
  plot_title_12M <- sprintf("%s Replicas | %s | Regression slope: %s", name_12M, patient, round(fit_12M$coefficients[[2]], 2))
  
  
  r <- plotManuallyColoredDensityScatterplot(plot_data = plot_data_12M, x_column_name = "A", y_column_name = "M", 
                                             bin_amount_per_axis = 1000, max_data_density_to_consider_per_cell = 100000, gradient_steps = 200,
                                             geom_point_size = 1,
                                             plot_title = plot_title_12M, x_label = "A", y_label = "M", 
                                             axis_title_size = 16, axis_text_size = 16, axis_ticks_length = 0.2, panel_background_size = 1,
                                             fixed_scale = 0)
  
  r <- r + coord_cartesian(xlim = c(min_fixed_scale_x, max_fixed_scale_x), ylim = c(min_fixed_scale_y, max_fixed_scale_y))
  
  
  # Se combinan los gráficos en una sola imagen
  combined_plots <- list(p, q, r)
  
  # Se guarda la imagen en formato TIFF
  # output_file <- sprintf("%s/MA_compare_replicas_%s.tiff", output_folder, patient)
  # tiff(file = output_file)
  # do.call(grid.arrange, c(combined_plots, ncol = 3))
  # dev.off()
  
  # Se guarda la imagen en formato PDF
  output_pdf <- sprintf("%s/MA_compare_replicas_%s.pdf", output_folder, patient)
  pdf(file = output_pdf, width = 20, height = 7)
  do.call(grid.arrange, c(combined_plots, ncol = 3))
  dev.off()
}

