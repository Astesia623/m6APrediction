#' DNA Sequence Encoding for m6A Prediction
#'
#' This function converts DNA sequences into a numerical feature matrix
#' by splitting the sequence into individual nucleotide positions.
#'
#' @param dna_strings A character vector of DNA sequences (must be same length)
#' @return A data.frame with nucleotide positions as factor columns
#' @examples
#' sequences <- c("ATCGT", "GCTAC", "CCGGA", "TTAAG")
#' encoded_data <- dna_encoding(sequences)
#' print(encoded_data)
#' @export
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A Sites for Multiple Sequences
#'
#' This function takes feature data and a trained model to predict m6A modification sites
#' for multiple sequences.
#'
#' @param ml_fit A trained random forest model object
#' @param feature_df A data.frame containing the required features:
#'   gc_content, RNA_type, RNA_region, exon_length, distance_to_junction,
#'   evolutionary_conservation, DNA_5mer
#' @param positive_threshold Probability threshold for positive prediction (default: 0.5)
#' @return A data.frame with original features plus prediction results
#' @examples
#' rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' example_data <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
#' predictions <- prediction_multiple(ml_fit = rf_model, feature_df = example_data)
#' head(predictions)
#' @importFrom stats predict
#' @import randomForest
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  required_cols <- c("gc_content", "RNA_type", "RNA_region", "exon_length",
                     "distance_to_junction", "evolutionary_conservation", "DNA_5mer")

  if(!all(required_cols %in% colnames(feature_df))) {
    missing_cols <- setdiff(required_cols, colnames(feature_df))
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  result_df <- feature_df

  if(any(nchar(result_df$DNA_5mer) != 5)) {
    stop("All DNA_5mer sequences must be exactly 5 characters long")
  }

  dna_encoded <- dna_encoding(result_df$DNA_5mer)
  result_df <- cbind(result_df, dna_encoded)

  result_df$RNA_type <- factor(result_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  result_df$RNA_region <- factor(result_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  nucleotide_levels <- c("A", "T", "C", "G")
  for(i in 1:5) {
    col_name <- paste0("nt_pos", i)
    if(col_name %in% colnames(result_df)) {
      result_df[[col_name]] <- factor(result_df[[col_name]], levels = nucleotide_levels)
    }
  }

  predicted_prob <- predict(ml_fit, newdata = result_df, type = "prob")
  result_df$predicted_m6A_prob <- predicted_prob[, "Positive"]
  result_df$predicted_m6A_status <- ifelse(result_df$predicted_m6A_prob > positive_threshold, "Positive", "Negative")
  result_df$predicted_m6A_status <- factor(result_df$predicted_m6A_status, levels = c("Negative", "Positive"))

  return(result_df)
}

#' Predict m6A Site for a Single Sequence
#'
#' This function takes individual feature values and a trained model to predict
#' m6A modification for a single sequence.
#'
#' @param ml_fit A trained random forest model object
#' @param gc_content Numeric, GC content value
#' @param RNA_type Character, RNA type (mRNA, lincRNA, lncRNA, pseudogene)
#' @param RNA_region Character, RNA region (CDS, intron, 3'UTR, 5'UTR)
#' @param exon_length Numeric, exon length
#' @param distance_to_junction Numeric, distance to junction
#' @param evolutionary_conservation Numeric, conservation score
#' @param DNA_5mer Character, 5-mer DNA sequence (must be exactly 5 characters)
#' @param positive_threshold Probability threshold for positive prediction (default: 0.5)
#' @return A named vector with prediction probability and status
#' @examples
#' rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' result <- prediction_single(
#'   ml_fit = rf_model,
#'   gc_content = 0.45,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 200,
#'   distance_to_junction = 50,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "ATCGT"
#' )
#' print(result)
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5){

  if(nchar(DNA_5mer) != 5) {
    stop("DNA_5mer must be exactly 5 characters long")
  }

  feature_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer
  )

  prediction_result <- prediction_multiple(ml_fit, feature_df, positive_threshold)

  returned_vector <- c(
    predicted_m6A_prob = prediction_result$predicted_m6A_prob,
    predicted_m6A_status = as.character(prediction_result$predicted_m6A_status)
  )

  return(returned_vector)
}
