# Permite realizar la comparación de señal de los péptidos (pico antigénico) de determinada proteína, comparando las señales 
# sin procesar, y las normalizadas por GRSN o por cuantiles. Esto puede realizarse para los diferentes pacientes del ensayo clínico.

# Se cargan los paquetes a utilizar
library(data.table)
library(ggplot2)
library(dplyr)

# Se establece la carpeta donde están los datos de las señales crudas
data_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/Array Roche E1224/00_inputs/02_serums/raw_data"

# Se establece la carpeta donde se encuentran los datos normalizados
data_folder2 <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R/Datos normalizados"

# Se establece la carpeta donde están los datos del diseño del experimento
data_folder3 <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/design"

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

# Proteína a analizar
protein <- "TcCLB.507953.200"

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

# Se carga el diseño del experimento, se seleccionan las columnas necesarias, y se renombran varias para mayor claridad
t_cruzi_map_data_file <- sprintf("%s/01_individual_serum_design_v3.tsv", data_folder3)
t_cruzi_map <- fread(t_cruzi_map_data_file, header = T, sep = "\t", na.strings = NULL)

t_cruzi_map_final <- t_cruzi_map %>%
  select(protein, peptide, truncated_peptide, start, origin) %>%
  rename(sequence = truncated_peptide, protein_identifier = protein)


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
  
  ## Se repite el mismo proceso para las otras dos normalizaciones (Quantile y GRSN), modificando las rutas correspondientes.
  
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
    summarise(mean_signal = mean(signal))
  
  GRSN_data$normalization <- normalization3
  
  GRSN_data <- GRSN_data %>%
    select(source, time, sequence, mean_signal, normalization)
  
  
  # Finalmente, se combinan los datos de las tres normalizaciones en una sola tabla 
  all_data <- rbind(raw_data, quantile_data, GRSN_data)
  
  # Se filtran los datos para mantener solo aquellos que tienen al menos una señal promedio mayor a 10000 (pico antigénico)
  all_data <- all_data %>%
    group_by(time, sequence) %>%
    filter(any(mean_signal > 10000))
  
  # Se fusionan los datos con el diseño del experimento
  all_data_w_id <- merge(all_data,t_cruzi_map_final[,c("sequence","protein_identifier", "start")], by = "sequence", all.x = TRUE)
  
  # Se calcula la media de la señal de pico para cada proteína
  all_data_w_id <- all_data_w_id %>%
    group_by(protein_identifier, time, normalization) %>%
    summarise(mean_peak_signal = round(mean(mean_signal),2))
  
  # Se filtran los datos para mantener solo aquellos relacionados con una proteína específica
  all_data_w_id <- all_data_w_id %>%
    filter(protein_identifier == protein)
  
  # Se establece el orden de los niveles del factor de tiempo
  all_data_w_id$time <- factor(all_data_w_id$time, levels = c("TMZ", "D65", "12M"))
  
  # Se establece el título del gráfico                                        
  plot_title <- sprintf(" %s | %s ", protein , patient)
  
  # Se crea el gráfico para comparar las normalizaciones de la proteína                                        
  p <- ggplot(all_data_w_id, aes(x = time, y = mean_peak_signal, color = normalization)) +
    geom_point(size = 1) +
    geom_line(aes(group = normalization), size = 0.7) +
    labs(
      title = plot_title,
      x = "Time",
      y = "Mean peak signal",
      color = "Normalization"
    ) +
    theme_minimal() 
  
  # Se guarda la imagen en formato PNG
  output_png <- sprintf("%s/Compare_peak_normalization_%s.png", output_folder, patient)
  ggsave(file = output_png, p , width = 10, height = 6, units = "in", dpi = 600)
  
}