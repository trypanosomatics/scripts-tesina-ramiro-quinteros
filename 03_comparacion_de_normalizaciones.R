# Permite realizar la comparación de señal de todos los péptidos, a través de gráficos recíprocos, comparando por un lado las señales
# crudas contra las normalizadas por cuantiles, y por otro las señales crudas contra las normalizadas por GRSN . Cada tiempo del ensayo 
# es marcado en un color diferente. Esto puede realizarse para los diferentes pacientes del ensayo clínico.

# Se cargan los paquetes a utilizar
library(gridExtra)
library(data.table)
library(ggplot2)
library(dplyr)

# Se establece la carpeta donde están los datos de las señales crudas
data_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/Array Roche E1224/00_inputs/02_serums/raw_data"

# Se establece la carpeta donde se encuentran los datos normalizados
data_folder2 <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R/Datos normalizados"

# Se establece la carpeta donde se guardará la comparación
output_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/normalization_compare"


# Se establecen los sufijos de los archivos para cada tipo de normalización
suffix <- "raw_data"
suffix2 <- "cruzi_normalized_data.tsv"
suffix3 <- "cruzi_GRSN_normalized_data.tsv"

# Normalizaciones utilizadas
normalization <- "Raw"
normalization2 <- "Quantile"
normalization3 <- "GRSN"

# Se especifica el paciente o listas de pacientes que se van a comparar
# patients<- "2-159-304"
patients<- c("1-005-101","1-049-116","1-124-150","1-160-163","1-227-184","1-264-203","2-098-277","2-132-293",
             "2-190-322","2-240-353","2-261-364", "1-052-118","1-105-143","1-143-157",
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

# Bucle principal para el procesamiento de cada paciente
for (patient in patients) {
  
  # Se determinan los archivos donde se encuentran los datos sin procesar
  raw_data_file_TMZ <-  sprintf("%s/%s_TMZ_%s",data_folder, patient, suffix)
  raw_data_file_D65 <-  sprintf("%s/%s_D65_%s",data_folder, patient, suffix)
  raw_data_file_12M <-  sprintf("%s/%s_12M_%s",data_folder, patient, suffix)
  
  # Se cargan los datos sin procesar
  raw_data_TMZ <- fread(raw_data_file_TMZ, header = T, sep = "\t", na.strings = NULL)
  raw_data_D65 <- fread(raw_data_file_D65, header = T, sep = "\t", na.strings = NULL)
  raw_data_12M <- fread(raw_data_file_12M, header = T, sep = "\t", na.strings = NULL)
  
  # Se combinan los datos sin procesar en una sola tabla
  raw_data <- rbind(raw_data_TMZ,raw_data_D65,raw_data_12M)
  
  # Se renombran las columnas para mayor claridad
  raw_data <- raw_data %>%
    rename(time=type)
  
  # Se calcula la media de la señal para cada secuencia, paciente y tiempo
  raw_data <- raw_data %>%
    group_by(sequence,source,time) %>%
    summarise(mean_signal= mean(signal))
  
  # Se agrega una columna indicando el tipo de normalización
  raw_data$normalization <- normalization
  
  # Se seleccionan las columnas relevantes para el análisis
  raw_data <- raw_data %>%
    select(source, time, sequence, mean_signal, normalization)
  
## Se repite el mismo proceso para las otras dos normalizaciones (Quantile y GRSN), modificando las rutas y etiquetas correspondientes.
  
  # Normalización por cuantiles
  quantile_data_file_TMZ <-  sprintf("%s/todos/%s_TMZ_%s",data_folder2, patient, suffix2)
  quantile_data_file_D65 <-  sprintf("%s/todos/%s_D65_%s",data_folder2, patient, suffix2)
  quantile_data_file_12M <-  sprintf("%s/todos/%s_12M_%s",data_folder2, patient, suffix2)
  
  quantile_data_TMZ <- fread(quantile_data_file_TMZ, header = T, sep = "\t", na.strings = NULL)
  quantile_data_D65 <- fread(quantile_data_file_D65, header = T, sep = "\t", na.strings = NULL)
  quantile_data_12M <- fread(quantile_data_file_12M, header = T, sep = "\t", na.strings = NULL)
  
  quantile_data <- rbind(quantile_data_TMZ,quantile_data_D65,quantile_data_12M)
  
  quantile_data <- quantile_data %>%
    group_by(sequence,source,time) %>%
    summarise(mean_signal= mean(signal))
  
  quantile_data$normalization <- normalization2
  
  quantile_data <- quantile_data %>%
    select(source, time, sequence, mean_signal, normalization)
  
  # Normalización por GRSN
  GRSN_data_file_TMZ <-  sprintf("%s/GRSN/%s_TMZ_%s",data_folder2, patient, suffix3)
  GRSN_data_file_D65 <-  sprintf("%s/GRSN/%s_D65_%s",data_folder2, patient, suffix3)
  GRSN_data_file_12M <-  sprintf("%s/GRSN/%s_12M_%s",data_folder2, patient, suffix3)
  
  GRSN_data_TMZ <- fread(GRSN_data_file_TMZ, header = T, sep = "\t", na.strings = NULL)
  GRSN_data_D65 <- fread(GRSN_data_file_D65, header = T, sep = "\t", na.strings = NULL)
  GRSN_data_12M <- fread(GRSN_data_file_12M, header = T, sep = "\t", na.strings = NULL)
  
  GRSN_data <- rbind(GRSN_data_TMZ,GRSN_data_D65,GRSN_data_12M)
  
  GRSN_data <- GRSN_data %>%
    group_by(sequence,source,time) %>%
    summarise(mean_signal= mean(signal))
  
  GRSN_data$normalization <- normalization3
  
  GRSN_data <- GRSN_data %>%
    select(source, time, sequence, mean_signal, normalization)
  
 
  # Finalmente, se combinan los datos de las tres normalizaciones en una sola tabla 
  all_data <- rbind(raw_data, quantile_data, GRSN_data)
  
  # Se filtran los datos para comparar únicamente las normalizaciones "Raw" y "Quantile"
  all_data1 <- filter(all_data, normalization == "Raw" | normalization == "Quantile")
  # Se repite el proceso para comparar las normalizaciones "Raw" y "GRSN"
  all_data2 <- filter(all_data, normalization == "Raw" | normalization == "GRSN")
  
  # Por las dudas se eliminan los grupos que no tienen ambas normalizaciones para cada secuencia, paciente y tiempo.
  all_data1 <- all_data1 %>%
    group_by(sequence, source, time) %>%
    mutate(normalization_count = n_distinct(normalization)) %>%
    filter(normalization_count == 2) %>%
    ungroup() %>%
    select(-normalization_count)
  
  sub_all_data1 <- as.data.table(all_data1)
  
  unique_normalizations1 <- unique(sub_all_data1$normalization)
  
  ## Se generan los gráficos comparativos para las normalizaciones "Raw" y "Quantile", donde se utilizan las señales promedio de cada
  ## normalización, señalando cada tiempo con un color distinto
  
  # Se cargan los datos necesarios para el gráfico en un data table
  if (length(unique_normalizations1) == 2) {
    plot_data1 <- data.table(x = sub_all_data1[normalization == unique_normalizations1[1]]$mean_signal,
                             y = sub_all_data1[normalization == unique_normalizations1[2]]$mean_signal,
                             x1 = sub_all_data1[normalization == unique_normalizations1[1]]$time,
                             x2 = sub_all_data1[normalization == unique_normalizations1[2]]$time)
    
    # Se define el título del gráfico y el de los ejes
    plot_title1 <- sprintf("%s vs %s Normalization | %s ", normalization, normalization2, patient)
    
    plot_x_label1 <- paste(normalization, "signal", sep = " ")
    plot_y_label1 <- paste(normalization2, "signal", sep = " ")
    
    # Se realiza el gráfico
    p <- ggplot(plot_data1, aes(x = x, y = y, color = x1)) +
      geom_point(size = 1) +
      labs(
        title = plot_title1,
        x = plot_x_label1,
        y =  plot_y_label1,
        color = "Time"
      ) +
      xlim(min_fixed_scale, max_fixed_scale) +
      ylim(min_fixed_scale, max_fixed_scale) +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.ticks.length = unit(0.2, "cm"),
        panel.background = element_rect(fill = "white", size = 1)
      )
  }  
  
  
  ## Se realiza el mismpro procedimientos pero para los datos crudos("raw") y la normalización "GRSN"
  all_data2 <- all_data2 %>%
    group_by(sequence, source, time) %>%
    mutate(normalization_count = n_distinct(normalization)) %>%
    filter(normalization_count == 2) %>%
    ungroup() %>%
    select(-normalization_count)
  
  sub_all_data2 <- as.data.table(all_data2)
  
  unique_normalizations2 <- unique(sub_all_data2$normalization)
  
  ## Se generan los gráficos comparativos para las normalizaciones "Raw" y "GRSN", donde se utilizan las señales promedio de cada
  ## normalización, señalando cada tiempo con un color distinto
  if (length(unique_normalizations2) == 2) {
    plot_data2 <- data.table(x = sub_all_data2[normalization == unique_normalizations2[1]]$mean_signal,
                             y = sub_all_data2[normalization == unique_normalizations2[2]]$mean_signal,
                             x1 = sub_all_data1[normalization == unique_normalizations2[1]]$time,
                             x2 = sub_all_data1[normalization == unique_normalizations2[2]]$time)
    
    plot_title2 <- sprintf("%s vs %s Normalization | %s ", normalization, normalization3, patient)
    
    plot_x_label2 <- paste(normalization, "signal", sep = " ")
    plot_y_label2 <- paste(normalization3, "signal", sep = " ")
    
    q <- ggplot(plot_data2, aes(x = x, y = y, color = x1)) +
      geom_point(size = 1) +
      labs(
        title = plot_title2,
        x = plot_x_label2,
        y =  plot_y_label2,
        color = "Time"
      ) +
      xlim(min_fixed_scale, max_fixed_scale) +
      ylim(min_fixed_scale, max_fixed_scale) +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.ticks.length = unit(0.2, "cm"),
        panel.background = element_rect(fill = "white", size = 1)
      )
    
  }
  
  # Se combinan los gráficos en una sola imagen
  combined_plots <- list(p,q)
  
  
  # Se guarda la imagen en formato PDF
  output_pdf <- sprintf("%s/Compare_normalization_color_time_%s.pdf", output_folder, patient)
  pdf(file = output_pdf, width = 13, height = 6)
  do.call(grid.arrange, c(combined_plots, ncol = 2))
  dev.off()
  
  # Se guarda la imagen en formato PNG
  output_png <- sprintf("%s/Compare_normalization_color_time_%s.png", output_folder, patient)
  png(file = output_png, width = 1000, height = 500)
  do.call(grid.arrange, c(combined_plots, ncol = 2))
  dev.off()
  
  # Se guarda la imagen en formato TIFF
  # output_tiff <- sprintf("Compare_times_with_circle_%s_%s.tiff", patient, normalization)
  # tiff(file = output_tiff)
  # do.call(grid.arrange, c(combined_plots, ncol = 2))
  # dev.off()
}
