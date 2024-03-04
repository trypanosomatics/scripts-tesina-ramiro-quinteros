# Este script permite la normalizacion por cuantiles o quantile normalization de todas las muestras presentes en el ensayo. 
# Se normalizan en conjunto todas las señales pertenecientes a un paciente, abarcado los 3 tiempos de la misma y ambas réplicas.
# También permite obtener estadísticas globales de estas nuevas señales (desviación estándar, moda).

#Se cargan los paquetes a utilizar
library(data.table)
library(preprocessCore) #quantile.normalization

#Se establece la carpeta donde se encuentran los datos de las señales crudas
data_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/Array Roche E1224"
setwd(data_folder)

#Se establece la carpeta donde se guardarán los datos normalizados
output_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R/Datos Normalizados/todos"

################-
#### CONFIG ####
################-
# Se da el path al diseño y se define si se quiere dejar todos o algunos de los orígenes de péptidos en la normalización
design_data_file <- sprintf("%s/00_inputs/01_design/01_individual_serum_design_v3.tsv", data_folder)
# design_origins_to_normalize <- c() #si esta vacio usa todos los peptidos
design_origins_to_normalize <- c("analyzed_proteome_cutoff_4_SD_2_pep",
                                 "peaks_from_cardiacos_cutoff_4_SD_2_pep",
                                 "peaks_from_serodiscordantes_cutoff_8_SD_2_pep",
                                 "peaks_from_BLAST_against_231_pident_95_length_similarity_80",
                                 "peaks_from_other_pathogens_cutoff_4_SD_2_pep") #si esta vacio usa todos los peptidos

#Se da el path a donde estan los archivos con las señales
raw_data_folder <- sprintf("%s/00_inputs/02_serums/raw_data", data_folder)
raw_data_suffix <- "_raw_data"

#Se obtiene la lista de sueros de pacientes
layout_data_file <- sprintf("%s/00_inputs/02_serums/parsed_layout_266_2.tsv", data_folder)
times_to_normalize <- c() #si esta vacio usa los tres tiempos, TMZ, D65 y 12M

#Se asigna la cantidad de decimales para las señales
output_signal_decimals <- 2

#Se define el sufijo de cada uno de los archivos a crear (uno por suero con su señal normalizada)
# output_suffix <- "_normalized_data.tsv"
output_suffix <- "_cruzi_normalized_data.tsv"

#Se define el nombre del archivo con la informacion estadística
# statistics_output_file <- "global_statistics.tsv"
# statistics_output_file <- "global_statistics_onlyCruziPeptides.tsv"
statistics_output_file <- sprintf("%s/global_statistics_onlyCruziPeptides.tsv",output_folder)
output_statistics_mode_decimals <- 2
output_statistics_sd_decimals <- 2

#############################-
#### AUXILIARY FUNCTIONS ####
#############################-
calculateMode <- function(x, decimals = 0) {
  x <- round(x, decimals)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

##############-
#### MAIN ####
##############-
###Se consigue la lista de sueros
layout_data <- fread(layout_data_file, header = T, sep = "\t", na.strings = NULL)

#Se filtra los tiempos si hace falta
if (length(times_to_normalize) > 0) {
  layout_data <- layout_data[time %in% times_to_normalize]
}

##Se ve si se tiene que filtrar los péptidos
if (length(design_origins_to_normalize) > 0) {
    design_data <- fread(design_data_file, header = T, sep = "\t", na.strings = NULL)
    design_data <- design_data[origin %in% design_origins_to_normalize]
    peptides_to_normalize <- unique(design_data$truncated_peptide) #uso los truncated_peptide aca porque son los que realmente se testearon
    
    rm(design_data) 
    gc()
} else {
    peptides_to_normalize <- c()
}

#Se extraen las sources
unique_parent_sources <- unique(layout_data$parent_source)

#Se aplica un for para recorrer cada uno de los sueros parentales
for (parent_source_for in unique_parent_sources) {
    # parent_source_for <- unique_parent_sources[1]
    
    unique_sources <- unique(layout_data[parent_source == parent_source_for]$source)

    ###Se recorre cada uno de los sueros y se agrega sus datos a una gran matriz que voy a usar al momento de normalizar
    for (source_for in unique_sources) {
        # source_for <- unique_sources[1]
        
        writeLines(source_for) #para ver por consola que el for avanza
        
        #Se carga el archivo de raw data
        raw_data_file <- sprintf("%s/%s%s", raw_data_folder, source_for, raw_data_suffix)
        raw_data <- fread(raw_data_file, header = T, sep = "\t", na.strings = NULL)
        
        #Si tengo que filtrar péptidos se hace
        if (length(peptides_to_normalize) > 0) {
            raw_data <- raw_data[sequence %in% peptides_to_normalize]    
        }
        
        #Cada combinacion de source + type + replica se normaliza independientemente, se necesita entonces ponerle una etiqueta única a cada grupo
        #type equivale al tiempo analizado
        raw_data[, group_aux := sprintf("%s_%s_%s", source, type, replica)]
        unique_groups <- unique(raw_data$group_aux) #creo una lista con cada grupo
        
        #Se Saca las columnas que ya no necesito y se ordena por secuencia dentro de cada grupo
        raw_data <- raw_data[, .(sequence, signal, group_aux)][order(group_aux, sequence)]
        
        #Si es la primera vez que se entra en este ciclo se tiene que crear la matriz de salida
        if (source_for == unique_sources[1]) {
            full_raw_data <- raw_data[group_aux == unique_groups[1], .(sequence)]
        }
        
        #Por cada grupo se agrega a una matriz grande que se va a usar para normalizar
        #El nombre de la columna va a ser el grupo
        for (group_for in unique_groups) {
            full_raw_data[[group_for]] <- raw_data[group_aux == group_for]$signal
        }
    }

    ###Aca se tiene una matriz con una columna sequence y después una columna por suero (source + type) + replica
    ##Se normaliza la matriz
    #normalize.quantiles es la funcion y necesita algo del tipo matrix, por lo que se tiene que convertir el data.table a matrix
    #Se saca "sequence" ya que sino normalize.quantiles querría normalizar esa columna tambien (y es una variable string)
    matrix_aux <- normalize.quantiles(as.matrix(full_raw_data[, -c("sequence")]))
    normalized_data <- as.data.table(matrix_aux)
    colnames(normalized_data) <- colnames(full_raw_data[, -c("sequence")])    
    normalized_data$sequence <- full_raw_data$sequence
    rm(full_raw_data) #Se saca lo que ya no es necesario para vaciar memoria
    rm(matrix_aux)
    gc()
    
    ###Se guardan los datos
    unique_times <- unique(layout_data$time) 
    unique_replicas <- unique(layout_data$replica)
    normalized_signals <- c() # Se inicializa un vector vacio
    
    for (time_for in unique_times) {
        # time_for <- unique_times[1]
        
        for (replica_for in unique_replicas) {
            # replica_for <- unique_replicas[1]
            
            group_for <- sprintf("%s_%s_%s", parent_source_for, time_for, replica_for) #recreo el grupo (la columna con los datos)
            setnames(normalized_data, group_for, "signal") #cambio el nombre de dicha columna a "signal" para poder trabajar con ella mas fácil
            
            normalized_data[, signal := round(signal, output_signal_decimals)] #redondeo signal
            
            if (replica_for == unique_replicas[1]) {
                output_normalized_data <- normalized_data[, .(sequence, source = parent_source_for, time = time_for, replica = replica_for, signal)]
            } else {
                output_normalized_data <- rbindlist(list(output_normalized_data, 
                                                         normalized_data[, .(sequence, source = parent_source_for, time = time_for, replica = replica_for, signal)]))
            }
            
            #Se guardan las signals para hacer estadística más adelante y renombro las columna
            normalized_signals <- c(normalized_signals, normalized_data$signal)
            setnames(normalized_data, "signal", group_for)
        }
        
        output_file <- sprintf("%s/%s_%s%s", output_folder, parent_source_for, time_for, output_suffix)
        write.table(output_normalized_data, file = output_file, col.names = T, row.names = F, sep = "\t", quote = T)
    }
}

###Se calculan estadísticas globales
mode_aux <- calculateMode(normalized_signals, decimals = output_statistics_mode_decimals)
sd_aux <- round(sd(normalized_signals), output_statistics_sd_decimals)

global_statistics <- data.table(mode = mode_aux,
                                sd = sd_aux)

write.table(global_statistics, file = statistics_output_file, col.names = T, row.names = F, sep = "\t", quote = T)
