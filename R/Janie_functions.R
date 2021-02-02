#'Read a .wig file
#'
#'This function reads a .wig file into R
#'
#'@param data_file Path to the data file
#'@param start Start of the nucleotide range you would like to pull out
#'@param end End of the nucleotide range you would like to pull out
#'@return A data frame with the nucleotide number and the reactivity
#' @export
read.wif = function(data_file, start, end){
  con = file(data_file)

  lines = readLines(con)

  close(con)

  lines = lines[-c(1:2)]

  lines = lines[c(start:end)]

  splitting.function = function(x){
    output = data.frame("N" = as.numeric(as.character(strsplit(x, split = "\t")[[1]][1])),
                        "Stops" = as.numeric(as.character(strsplit(x, split = "\t")[[1]][2])))
  }

  df.list = lapply(lines, FUN = splitting.function)

  df = dplyr::bind_rows(df.list)

  print(head(df))

  df = df
}

#'Quickly compiles stops data for analysis in other programs
#'
#'This function reads a .wig file and a fasta file into R, pulls out a specified subset of the data, cleans stop counts associated
#'with Gs and Us, Normalizes the data using the Winsorization method, calculates sliding window averages, and writes a csv file
#'for subsequent analysis. It also plots the data as a png for convenient assessments of the result.
#'
#'@param file.wif Path to the .wif formatted file
#'@param file.fasta Path to the fasta file
#'#'@param start End of the nucleotide range you would like to pull out
#'@param end End of the nucleotide range you would like to pull out
#'@param strand Strand of the genome you want to analyze. Options: Default = "+" for the forward strand and "-" for the reverse strand
#'@param Name Name of the RNA. (This will be the prefix on the file that is saved). Default = "RNA".
#'@param Nucleic_acid The nucleic acid you want the output to be. Options: Default = "RNA" to print Us in the final file and "DNA" to print "Ts".
#'@param Winsorization Use the Winsorization method to normalize the data. Options: Default = TRUE to do the normalization or FALSE to skip the normalization.
#'@param Normalization_quantile Quartile you want to do the Winsorization normalization to.
#'@param Window Window type you want to use when calculating mean reactivities in sliding windows. Options: Default = "center" to apply the mean to the nucleotide at the center of the window, or "first" to apply the mean to the first nucleotide in the window.
#'@param Window_size Width of the sliding window that you want to apply. Default = 10 nucleotides
#'@param custum_windows A list of custum windows you want to calculate the mean for. Default = FALSE or a list of vectors in the format custom_windows = list(c(window1.start:window1.end), c(window2.start:window2.end)). Example: custom_windows = list(c(31:42), c(50:200)). Make sure that you use the relative location of the window on the RNA, not the absolute location on the genome.
#'@return A csv file, a plot, and a dataframe
#' @export
stops.extract = function(file.wif,
                         file.fasta,
                         start,
                         end,
                         strand = "+",
                         Name = "RNA",
                         Nucleic_acid = "RNA",
                         Winsorization = TRUE,
                         Normalization_quantile = 0.95,
                         Window = "center",
                         Window_size = 10,
                         Use_custom_windows = FALSE,
                         custom_windows = list()){
  ####Read in .wig file####

  list.files()

  df = JFlabR::read.wif(file.wif , start, end)

  print("Read in .wig file")

  ####Load in RNA sequence data####

  Nucleotide = seqinr::read.fasta(file.fasta)[[1]][c(start:end)]

  print("Load in RNA sequence data")

  ####Make reverse strand####

  df$Nucleotide = Nucleotide

  if (strand == "-"){
    Reverse = c()
    for (i in 1:length(Nucleotide)){
      if (Nucleotide[i] == "a"){Reverse[i] <- "t"}
      if (Nucleotide[i] == "c"){Reverse[i] <- "g"}
      if (Nucleotide[i] == "g"){Reverse[i] <- "c"}
      if (Nucleotide[i] == "t"){Reverse[i] <- "a"}
    }
    Nucleotide = rev(Reverse)
    print("Make reverse strand")
  }

  ####Capitalize the letters and convert Ts to Us####

  if (Nucleic_acid == "RNA"){
    Nucleotide[which(Nucleotide == "t")] <- "U"
  }else{
    Nucleotide[which(Nucleotide == "t")] <- "T"
  }

  Nucleotide[which(Nucleotide == "a")] <- "A"
  Nucleotide[which(Nucleotide == "c")] <- "C"
  Nucleotide[which(Nucleotide == "g")] <- "G"

  print("Capitalize the letters and convert Ts to Us")

  ####Add the sequence strand to df####

  df$Nucleotide = Nucleotide

  print("Add the sequence strand to df")

  ####Remove stop counts associated with Gs and Us####

  df[which(df$Nucleotide == "U"), 2] <- NA
  df[which(df$Nucleotide == "G"), 2] <- NA

  print("Remove stop counts associated with Gs and Us")

  ####Normalize with 95% Winsorization####

  if (Winsorization){
    data.95 = quantile(df$Stops, Normalization_quantile, na.rm = TRUE)

    #hist(df$Stops)
    #abline(v = data.95)

    df$Stops[which(df$Stops > data.95)] <- data.95

    #hist(df$Stops)
    #abline(v = data.95)

    df$Stops <- df$Stops/data.95

    print(paste("Performing", data.95, "quartile normalization"))
  }

  ####Calculate average for a sliding widow####

  if (Window == "center"){
    Window.mean = c()
    edges = ceiling(0.5*Window_size)
    for (i in Window_size:(length(df$N)-(Window_size))){
      Window.mean[i] <- mean(df$Stops[(i - edges):(i + edges)][-which(is.na(df$Stops[(i - edges):(i + edges)]))])
    }
    Window.mean = c(Window.mean, rep(NA, Window_size))
    df$Window.mean <- Window.mean
  }

  #plot(df$N, df$Window.mean)

  if (Window == "first"){
    Window.mean = c()
    for (i in 1:(length(df$N)-(Window_size))){
      Window.mean[i] <- mean(df$Stops[(i):(i + Window_size -1)][-which(is.na(df$Stops[(i):(i + Window_size -1)]))])
    }
    Window.mean = c(Window.mean, rep(NA, Window_size))
    df$Window.mean <- Window.mean
  }

  #plot(df$N, df$Window.mean)

  print("Calculate average for a sliding widow")

  ####Calculate mean transcript reactivity####

  df$Transcripthead.mean <- mean(df$Stops[-which(is.na(df$Stops))])

  print("Calculate mean transcript reactivity")

  ####Make N start at 1####

  df$N <- 1:length(df$N)

  print("Make N start at 1")

  ####Calculate means for custom windows####

  #custom_windows = list(c(31:42), c(50:200))

  if (Use_custom_windows){
    Labels <- rep(NA, length(df$N))
    Custom.window.mean <- rep(NA, length(df$N))
    for (i in 1:length(custom_windows)){
      a <- mean(df$Stops[custom_windows[[i]]][-which(is.na(df$Stops[custom_windows[[i]]]))])
      Custom.window.mean[custom_windows[[i]]] <- a
      Labels[custom_windows[[i]]] <- paste("N =", min(custom_windows[[i]]), "to", max(custom_windows[[i]]))
    }
    df$Custom.window.mean <- Custom.window.mean
    df$Labels <- Labels
    print("Calculate means for custom windows")
  }

  ####Plot data####

  head(df)

  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = df$Transcripthead.mean[1], size = 1, color = "blue") +
    ggplot2::geom_point(data = df, mapping = ggplot2::aes(x = N, y = Stops), color = "grey") +
    ggplot2::geom_line(data = df, mapping = ggplot2::aes(x = N, y = Window.mean), size = 1) +
    ggplot2::theme_classic()

  if (Use_custom_windows){
    ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = df$Transcripthead.mean[1], size = 1, color = "blue") +
      ggplot2::geom_point(data = df, mapping = ggplot2::aes(x = N, y = Stops), color = "grey") +
      ggplot2::geom_line(data = df, mapping = ggplot2::aes(x = N, y = Window.mean), size = 1) +
      ggplot2::geom_line(data = df, mapping = ggplot2::aes(x = N, y = Custom.window.mean, color = Labels), size = 1) +
      ggplot2::theme_classic()
  }

  print("Plot data")

  ####Save plot####

  ggplot2::ggsave(paste(Name, "_", start, "_to_", end, ".png", sep = ""), width = 5, height = 2, dpi = 300, scale = 2)

  print("Save plot")

  ####Save data####

  write.csv(df, paste(Name, "_", start, "_to_", end, ".csv", sep = ""), row.names = FALSE)

  print("Write csv")

  ####Make a data frame be an output####

  output = df

  print("Done")

  ####End####
}
