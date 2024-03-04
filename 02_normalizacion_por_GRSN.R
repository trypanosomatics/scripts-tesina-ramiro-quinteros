# Este script permite la normalizacion por GRSN de todas las muestras presentes en el ensayo. Utiliza el mismo procedimiento que la
# normalización por cuantiles, utilizando la función de normalización por GRSN brindada por Carl Pelz
# Se normalizan en conjunto todas las señales pertenecientes a un paciente, abarcado los 3 tiempos de la misma y ambas réplicas.
# También permite obtener estadísticas globales de estas nuevas señales (desviación estándar, moda).

#Se cargan los paquetes a utilizar
library(data.table)

#Se establece la carpeta donde se encuentran los datos de las señales crudas
data_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/Array Roche E1224"
setwd(data_folder)

#Se establece la carpeta donde se guardarán los datos normalizados
output_folder <- "C:/Users/Ramiro/OneDrive - Universidad de San Martin/Tesina/R/Datos Normalizados/GRSN"

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
# output_suffix <- "GRSN_normalized_data.tsv"
output_suffix <- "_cruzi_GRSN_normalized_data.tsv"

#Se define el nombre del archivo con la informacion estadística
# statistics_output_file <- "global_statistics.tsv"
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
  peptides_to_normalize <- unique(design_data$truncated_peptide) #tengo que usar los truncated aca porque son los que realmente se testearon
  
  rm(design_data) #esto saca a design_data de memoria para hacer mas lugar para despues
  gc()
} else {
  peptides_to_normalize <- c()
}


#Se extraen las sources
unique_parent_sources <- unique(layout_data$parent_source)
# unique_parent_sources <- c("1-005-101") #DEBUG
# unique_parent_sources <- unique_parent_sources[1:3] #DEBUG

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
    
    # Cada combinacion de source + type + replica se normaliza independientemente, se necesita entonces ponerle una etiqueta única a cada grupo
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
  
  ###Aca ya se tiene una matriz con una columna sequence y despues una columna por suero (source + type) + replica
  # Main GRSN implementation script obtenido de Carl Pelz

  #############################-
  #### AUXILIARY FUNCTIONS ####
  #############################-
  GRSN <- function(data,           # exprSet from affy package or matrix with column for each sample.
                   width=15,       # Width (inches) of diagnostic plot.
                   height=5,       # Height (inches) of diagnostic plot.
                   pointsize=16,   # Point size for text of diagnostic plot.
                   filetype="png", # Diagnostic plot format "png", "wmf", or "postscript".
                   ReportName=NA,  # Name for Diagnostic plots.
                   count=5000,     # Size of Global Rank-invariant Set to use.
                   f=0.25)         # Smoother parameter for lowess.
  {
    # Check the class of the input data.
    
    if (max(grepl("matrix",class(data))) == 1)
    {
      rawData <- data
    } else
    {
      stop("data parameter is not a valid type!")
    }
    
    # } 
    
    
    # Make sure that the input data is not log scale.
    if (max(rawData) > 31)
    { # Assume that input data is not log scale.
      isItLogScaled <- FALSE
    } else
    { # Assume that input data is log scale and convert it.
      isItLogScaled <- TRUE
      rawData = 2^rawData
    }
    
    # Find Affymetrix(R) control probe sets by looking for 
    # probe set IDs starting in "AFFY".
    affyIdx <- grep ("^AFFX", attr(rawData, "dimnames")[[1]])
    
    # Data to normalize.
    adjust <- max(0, (0.25 - min(rawData)))
    M1 <- log2(rawData + adjust)
    
    # Get the average of the reference set.
    # Do a trimmed mean to be robust, but eliminate the "artifact" that 
    # shows up when doing median on an odd number of samples.
    Mavg <- apply(M1[, ], 1, mean, trim=0.25)
    
    # New method for a global invariant set.
    total <- dim(M1)[[1]]
    idx <- 1:total
    subSet <- 1:total
    
    # Exclude Affy control probe sets from 
    # approximate global rank invaraint set (GRiS).
    if (length(affyIdx) > 0)
    {
      total <- total - length(affyIdx)
      idx <- idx[-affyIdx]
      subSet <- subSet[-affyIdx]
    }
    
    # Calculate number of probe sets to exclude at each iteration.
    discardNumber <- (total - count) / 4
    
    ### Main iteration loop to get approximate GRiS. ###
    while (TRUE)
    {
      total <- floor(max(total - discardNumber, count))
      M2 <- cbind(apply(M1[idx, ], 2, rank))
      V2 <- apply(M2, 1, var)
      subSet <- order(V2, decreasing=FALSE)[1:total]     
      idx <- idx[subSet]
      if (total == count) break
    }
    invariantIdx <- idx
    
    # Use invariant set to normalize all samples to the average.
    Mnew <- NULL
    x <- Mavg
    for (b in 1:dim(M1)[[2]])
    {
      if (!is.na(ReportName))
      {
        PrintToFile <- GraphSetup(width=width, height=height, 
                                  pointsize=pointsize, 
                                  ReportName=paste(ReportName, b, sep=""),
                                  filetype=filetype)
      }
      
      if (PrintToFile)
      {
        # Plot three graphs side-by-side.
        layout(matrix(c(1,2,3), nrow=1, ncol=3, byrow = TRUE))
        par(mar=c(5,5,2,1)) # c(bottom, left, top, right).
      }
      
      y <- M1[,b]
      
      ### M vs. A transformed data.  ###
      M <- y-x
      A <- (y+x)/2
      
      ### Lowess curve based on M vs. A transformed data. ###
      curve <- lowess(x=A[invariantIdx], y=M[invariantIdx], f=f)
      
      ### Create evenly space lookup from calibration curve. ###
      aCurve <- curve[[1]]
      mCurve <- curve[[2]]
      steps <- 1000
      sampleMin <- min(A)
      sampleMax <- max(A)
      step <- (sampleMax - sampleMin) / steps
      position <- seq(sampleMin, sampleMax, length=steps + 1)
      adjust <- array(0,c(steps+1))
      count <- length(aCurve)
      
      idxL <- 1
      idxR <- 2
      for (i in 1:(steps + 1))
      {
        while (idxR < count && position[i] > aCurve[idxR])
        {
          idxR <- idxR + 1
        }
        while ((idxL + 1) < idxR && position[i] > aCurve[idxL + 1])
        {
          idxL <- idxL + 1
        }
        while (idxR < count && aCurve[idxL] >= aCurve[idxR])
        {
          idxR <- idxR + 1
        }
        if (aCurve[idxL] < aCurve[idxR])
        {
          adjust[i] <- (((mCurve[idxR] - mCurve[idxL])/(aCurve[idxR] - aCurve[idxL]))
                        *(position[i] - aCurve[idxL]) + mCurve[idxL])
        }
      }
      
      ### Apply lookup to data.  Can be applied to transformed or untransformed data. ###
      yPrime <- y - adjust[(A - sampleMin) / step + 1.5]
      mPrime <- yPrime - x
      
      Mnew <- cbind(Mnew, yPrime)
      if (PrintToFile)
      {
        sampleName <- attr(rawData,"dimnames")[[2]][b]
        
        plot(x=A, y=M, pch=".", col="blue", 
             main= paste("Scatter Plot ", b, " vs. Reference Set", sep=""), 
             sub="", 
             xlab=paste("(log2(", sampleName, ")+log2(ref))/2", sep = ""), 
             ylab=paste("log2(", sampleName, ")-log2(ref)", sep = ""))
        lines(x=c(-5, 20), y=c(0, 0), col="black")
        
        plot(x=A[invariantIdx], y=M[invariantIdx], pch=".", col="blue", 
             main= paste("Invariant Scatter Plot ", b, " vs. Reference Set ",sep=""), 
             sub="", 
             xlab=paste("(log2(", sampleName, ")+log2(ref))/2", sep = ""), 
             ylab=paste("log2(", sampleName, ")-log2(ref)", sep = ""))
        points(x=A[invariantIdx], y=mPrime[invariantIdx], pch=".", col="red")
        lines(x=c(-5, 20), y=c(0, 0), col="black")
        
        lines(x=position, y=adjust, lwd=2, col="green")
        
        plot(x=A, y=mPrime, pch=".", col="red", 
             main= paste("Normalized Scatter Plot ", b, " vs. Reference Set ", sep=""), 
             sub="", 
             xlab=paste("(log2(", sampleName, ")+log2(ref))/2", sep = ""), 
             ylab=paste("log2(", sampleName, ")-log2(ref)", sep = ""))
        lines(x=c(-5, 20), y=c(0, 0), col="black")
      }
      
      if (PrintToFile)
      {
        # Stop outputting to file.
        dev.off()
      }
    }
    
    if (isItLogScaled == FALSE)
    {
      # Convert.
      Mnew <- 2^Mnew
    }
    
    # Check the class of the input data.
    # if (class(data) == "exprSet")
    # {
    #   # Update expression set values.
    #   attr(data, "exprs")[,] <- Mnew[,]
    # } else if (class(data) == "matrix")
    # {
    # Update matrix values.
    data[,] <- Mnew[,]
    # } else if (class(data)[[1]] == "ExpressionSet")
    # {
    #   exprs(data) <- Mnew[,]
    # }
    
    return(data)
  }
  # Terminamos de definir función GRSN
  
  ##Se normaliza la matriz
  #Se aplica normalizacion GRSN  
  matrix_aux <- as.matrix(full_raw_data[, -c("sequence")])
  GRSN_data <- GRSN(matrix_aux)
  normalized_data <- as.data.table(GRSN_data)
  colnames(normalized_data) <- colnames(full_raw_data[, -c("sequence")])    
  normalized_data$sequence <- full_raw_data$sequence  

  ###Se guardan los datos
  unique_times <- unique(layout_data$time) 
  unique_replicas <- unique(layout_data$replica)
  normalized_signals <- c() #inicializo un vector vacio

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

###Se calculan estadísticas globaless
statistics_output_file <- sprintf("%s/global_statistics_onlyCruziPeptides.tsv",output_folder)
output_statistics_mode_decimals <- 2
output_statistics_sd_decimals <- 2


mode_aux <- calculateMode(normalized_signals, decimals = output_statistics_mode_decimals)
sd_aux <- round(sd(normalized_signals), output_statistics_sd_decimals)

global_statistics <- data.table(mode = mode_aux,
                                sd = sd_aux)

write.table(global_statistics, file = statistics_output_file, col.names = T, row.names = F, sep = "\t", quote = T)




