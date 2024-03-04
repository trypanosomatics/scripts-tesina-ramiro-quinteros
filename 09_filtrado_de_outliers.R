# Se filtran los outliers, se descartan aquellos que poseen un valor mayor a 0.16 de puntaje de variación, excepto aquellos
# con variación absoluta < 1000, o variación relativa < 1.5, ya que los mismo no son considerados outliers.

# Se cargan los paquetes a utilizar
library(dplyr)
library(tools)

# Se establece la carpeta donde están los datos de las señales crudas con sus puntajes de variación, provenientes del script anterior
data_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/raw_data_mean_replicas_w_variation"
setwd(data_folder)

# Se establece la carpeta donde se guardará el resultado
output_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/raw_data_mean_replicas_filtered_w_formula"

# Se obtiene la lista de archivos de muestras a analizar
samples<-list.files(data_folder)

# Se itera sobre cada muestra
for(sample in samples){
  # Se cargan los datos del archivo
  raw_data <- read.table(sample, header=TRUE, sep = "", na.strings = NULL)
  # Se filtran los outliers utilizando el puntaje de variación, y la variación absoluta y relativa
  raw_data_filtered <- raw_data %>%
    filter(variation_score < 0.16 | relative_variation < 1.5 | absolute_variation < 1000)
  
  # Se obtiene el nombre del archivo sin la extensión
  file_name_without_ext <- file_path_sans_ext(sample)
  
  # Se establece el nombre del archivo de salida filtrado y se guarda
  raw_data_filtered_output_file <- sprintf("%s/%s_filtered.tsv", output_folder, file_name_without_ext)
  write.table(raw_data_filtered, file = raw_data_filtered_output_file, sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)
}

