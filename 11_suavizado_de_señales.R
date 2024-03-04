# Se realiza el suavizado de las señales promedio de cada réplica mediante el método de “rolling mean” o “media móvil”, 
# que consiste en calcular cada señal como el promedio de un conjunto de valores adyacentes en una ventana deslizante que 
# se mueve a lo largo de la serie. Este suavizado se realiza para todas las muestras de todos los pacientes, para sus señales
# luego del filtrado de outliers.
# Debido a que un péptido puede tener diferentes entornos, es decir, puede estar presente en distintas proteínas con péptidos
# adyacentes diferentes en cada una, este péptido puede tener más de una señal que lo representa. Debido a este problema de 
# un péptido con señales diferentes, las señales obtenidas para un péptido a partir de este punto se van a tener un cuenta en 
# el contexto de cada proteína a la que pertenece.

#Se cargan los paquetes a utilizar
library(data.table)
library(zoo) #Para rollmean y rollmedian

# Se establece la carpeta donde se encuentran los datos de procesamiento
data_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing"
setwd(data_folder)

################-
#### CONFIG ####
################-

# Archivo que contiene el diseño del experimento
design_data_file <- sprintf("%s/design/01_individual_serum_design_v3.tsv", data_folder)

# Orígenes a considerar en el diseño
design_origins_to_parse <- c()

# design_origins_to_parse <- c("analyzed_proteome_cutoff_4_SD_2_pep",
# "peaks_from_cardiacos_cutoff_4_SD_2_pep",
# "peaks_from_serodiscordantes_cutoff_8_SD_2_pep",
# "peaks_from_BLAST_against_231_pident_95_length_similarity_80",
# "peaks_from_other_pathogens_cutoff_4_SD_2_pep") #si esta vacio usa todos


#Se da el path a donde están los archivos con las señales filtradas
filtered_data_folder <- sprintf("%s/raw_data_mean_replicas_filtered_w_formula", data_folder)
filtered_data_suffix <- "_raw_data_filtered.tsv"

#Se obtiene la lista de sueros de pacientes
layout_data_file <- sprintf("%s/design/parsed_layout_266_2.tsv", data_folder)
times_to_parse <- c() #si esta vacio usa los tres tiempos, TMZ, D65 y 12M

#Estas son las opciones para el rolling mean y el rolling median (smootheado de signals en una proteína mediante la funcion smoothVector)
# Determino el tamaño de la ventana y opción en el borde
smoothing_median_window_size <- 3
smoothing_mean_window_size <- 3
smooth_borders_option <- "zeros"

# Sufijo para los archivos de salida
output_suffix <- "_cruzi_protein_data.tsv"

# Decimales para la señal de salida
output_signal_decimals <- 2


#############################-
#### AUXILIARY FUNCTIONS ####
#############################-
# Función para suavizar un vector, los péptidos en este caso
smoothVector <- function(vector, median_window_size = 3, mean_window_size = 0, borders = "zeros") {
  # borders puede ser "repeat" o "zeros"
  
  if (median_window_size > 0) {
    # Llena los bordes para la mediana
    if (borders == "repeat") {
      #Rellena los lados con el primer y último número para tener la misma cantidad de datos después del suavizado
      prefix <- rep(vector[1], floor((median_window_size - 1) / 2))
      suffix <- rep(vector[length(vector)], ceiling((median_window_size - 1) / 2))        
    } else if (borders == "zeros") {
      #Rellena los lados con 0 para tener la misma cantidad de datos después del suavizado
      ### Esto aplana un poco los bordes.
      prefix <- rep(0, floor((median_window_size - 1) / 2))
      suffix <- rep(0, ceiling((median_window_size - 1) / 2))
    } else {
      writeLines("WARNING: Incorrect border option.")
    }
    vector_aux <- c(prefix, vector, suffix)
    
    # Calcular la mediana móvil
    smoothed_vector <- round(rollmedian(vector_aux, median_window_size), 3)    
  } else {
    smoothed_vector <- vector #this is because the name change
  }
  
  if (mean_window_size > 0) {
    # Llena los bordes para la media
    if (borders == "repeat") {
      #Rellena los lados con el primer y último número para tener la misma cantidad de datos después del suavizado
      prefix <- rep(smoothed_vector[1], floor((mean_window_size - 1) / 2))
      suffix <- rep(smoothed_vector[length(smoothed_vector)], ceiling((mean_window_size - 1) / 2))        
    } else if (borders == "zeros") {
      #Rellena los lados con 0 para tener la misma cantidad de datos después del suavizado
      ### Esto aplana un poco los bordes.
      prefix <- rep(0, floor((mean_window_size - 1) / 2))
      suffix <- rep(0, ceiling((mean_window_size - 1) / 2))
    } else {
      writeLines("WARNING: Incorrect border option.")
    }
    vector_aux <- c(prefix, smoothed_vector, suffix)
    
    # Calcula la media móvil
    smoothed_vector <- round(rollmean(vector_aux, mean_window_size), 3)        
  }
  
  smoothed_vector
}

##############-
#### MAIN ####
##############-
#Se consigue la lista de sueros
layout_data <- fread(layout_data_file, header = T, sep = "\t", na.strings = NULL)

#Se filtra los tiempos si hace falta
if (length(times_to_parse) > 0) {
  layout_data <- layout_data[time %in% times_to_parse]
}

#Se obtiene la lista de sueros de pacientes a utilizar
unique_sources <- unique(layout_data$source)

###Se lee el diseño del experimento y se filtra según orígenes a utilizar
design_data <- fread(design_data_file, header = T, sep = "\t", na.strings = NULL)

if (length(design_origins_to_parse) > 0) {
  design_data <- design_data[origin %in% design_origins_to_parse]
}

# Obtengo la secuencia de los péptidos a utilizar
sequences_in_design <- unique(design_data$truncated_peptide)

# Se quedan solo las columnas relevantes para este proceso
design_data <- design_data[, .(protein, start, sequence = peptide, truncated_sequence = truncated_peptide)]

# Se eliminan datos duplicados
design_data <- dplyr::distinct(design_data)

### Por cada suero se cargan los datos filtrados y se reconstruyen las proteinas
for (source_for in unique_sources) {
  # source_for <- unique_sources[1]
  
  writeLines(source_for)
  
  filtered_data_file <- sprintf("%s/%s%s", filtered_data_folder, source_for, filtered_data_suffix)
  filtered_data <- fread(filtered_data_file, header = T, sep = "\t", na.strings = NULL)
  setnames(filtered_data, "sequence", "truncated_sequence")
  filtered_data <- filtered_data[, .(truncated_sequence, source, time, mean_signal)]
  
  # Se saca cualquier peptido que no está en el diseño a utilizar
  filtered_data <- filtered_data[truncated_sequence %in% sequences_in_design]
  
  # Se Agregan los datos de proteina a las signals
  #Se usa "truncated_sequence" ya que es la secuencia que realmente se testeó en el microarray
  filtered_data <- merge(filtered_data,
                         design_data,
                         by = "truncated_sequence",
                         allow.cartesian = T)
  
  # Se odenan los datos
  filtered_data <- filtered_data[order(source, time, protein, start)]
  
  # Se suaviza la señal
  smoothed_filtered_data_aux <- filtered_data[, .(mean_smoothed_signal = smoothVector(mean_signal, smoothing_median_window_size, smoothing_mean_window_size, smooth_borders_option)),
                                              by = .(source, time, protein)]
  
  filtered_data$mean_smoothed_signal <- smoothed_filtered_data_aux$mean_smoothed_signal
  
  # Se saca los que tienen señal 0, de esta forma elimino los orígenes que solo tienen 1 peptido y para los que no sirve rolling mean
  #  (Vienen de alanin-scan y esos orígenes)
  filtered_data <- dplyr::filter(filtered_data, mean_smoothed_signal != 0)
  
  # Se redondea las señales a 2 decimales
  filtered_data <- filtered_data[, .(mean_smoothed_signal = round(mean_smoothed_signal, output_signal_decimals)),
                                 by = .(source, time, protein, start)]
  
  # Se agregan los datos de la secuencia
  filtered_data <- merge(filtered_data,
                         design_data[, .(protein, start, sequence, truncated_sequence)],
                         by = c("protein", "start"))
  
  # Se ordenan las columnas y se guarda el archivo
  filtered_data <- filtered_data[order(source, time, protein, start)]
  setcolorder(filtered_data, c("source", "time", "protein", "start",
                               "mean_smoothed_signal", "sequence",
                               "truncated_sequence"))
  
  output_file <- sprintf("%s/data_w_rolling_mean/%s%s", data_folder, source_for, output_suffix)
  write.table(filtered_data, file = output_file, col.names = T, row.names = F, sep = "\t", quote = T)
}
