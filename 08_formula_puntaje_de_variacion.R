# Se calcula el puntaje de variación con el uso de la fórmula y de esta forma se pueden eliminar outliers, considerando error relativo y absoluto,
# y normalizando cada métrica con el máximo de la misma
# Fórmula final: Variación absoluta/Máximo Variación absoluta  * α + Variación relativa/Variación relativa * β


# Se cargan los paquetes a utilizar
library(dplyr)

# Se establece la carpeta donde están los datos de las señales crudas con sus variaciones, provenientes del script anterior
data_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/raw_data_mean_replicas_w_variation/no_formula"
setwd(data_folder)

# Se establece la carpeta donde se guardará el resultado
output_folder <- "C:/Users/ramik/OneDrive - Universidad de San Martin/Tesina/R final/Data_processing/raw_data_mean_replicas_w_variation"


# Se obtienen los máximos de varación relativa y absoluta para relativizar
samples<-list.files(data_folder)
max_relative_variation <- data.frame()
max_absolute_variation <- data.frame()
for(sample in samples){
  raw_data <- read.table(sample, header=TRUE, sep = "", na.strings = NULL)


  max_relative_variation_aux <- max(raw_data$relative_variation)
  max_relative_variation <- rbind(max_relative_variation_aux,max_relative_variation)

  max_absolute_variation_aux <- max(raw_data$absolute_variation)
  max_absolute_variation <- rbind(max_absolute_variation_aux,max_absolute_variation)
}
print(max(max_relative_variation[,1]))
print(max(max_absolute_variation[,1]))

# Los máximos obtenidos fueron
max_absolute_variation <- 65525.84
max_relative_variation <- 9704.382


# Se establecen a y b deseados
# a = Peso de variación absoluta
# b = Peso de variación relativa
a <- 1
b <- 50

# Se obtiene la lista de archivos de muestras a analizar
samples<-list.files(data_folder)

# Se itera sobre cada muestra
for(sample in samples){
  # Se cargan los datos crudos del archivo
  raw_data <- read.table(sample, header=TRUE, sep = "", na.strings = NULL)
  
  # Se calcula el puntaje de variación
  raw_data_formula <- raw_data %>%
    group_by(sequence,source,time) %>%
    mutate(variation_score = (absolute_variation/max_absolute_variation)*a + (relative_variation/max_relative_variation)*b) %>%
    ungroup()
  
  # Se guarda el resultado en un archivo de salida para cada muestra
  raw_data_formula_output_file <- sprintf("%s/%s", output_folder, sample)
  write.table(raw_data_formula, file = raw_data_formula_output_file, sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)
}




