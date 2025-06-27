#' @title Perform adaptive-weighted transcriptomic meta-analysis using AWmeta
#' @description This function automates a meta-analysis workflow for differential
#'              gene expression. It loads multiple datasets, automatically detects
#'              data type (continuous/discrete) to assign the appropriate DE analysis
#'              method, combines p-values using AW-Fisher, and calculates a summary
#'              fold-change using AW-REM model.
#'
#' @param raw.data.dir A path to the raw expression data files.
#' @param raw.clin.dir A path to the clinical/phenotype data files.
#' @raw.sep The field separator character. Values on each line of the raw expression and clinical/phenotype data file are separated by this character.
#' @param DE.method A character vector specifying the DE analysis method(s).
#'                  - **Single Method (length 1):** e.g., `"limma"`. Applies this method to all studies.
#'                  - **Automatic Detection (length 2):** e.g., `c("limma", "DESeq2")`. Provide one
#'                    continuous method ("limma", "sam") and one discrete method ("edgeR", "DESeq2", "limmaVoom").
#'                    The function will auto-detect if a study's data is integer-based (discrete) or
#'                    decimal-based (continuous) and apply the corresponding method.
#' @param compare.group A character vector of length 2 specifying the names of the two groups to compare
#'                      in the clinical data (e.g., c("control", "PD")).
#' @param ref.level A character string specifying which of the two groups in `compare.group` is the
#'                  reference or baseline level (e.g., "control").
#' @param paired A logical value (TRUE/FALSE) indicating whether the samples are paired.
#' @param core.num An integer specifying the number of CPU cores to use for parallel computation.
#'
#'
#' @examples
#' \dontrun{
#' # --- Create dummy files for demonstration ---
#' # Study 1: Microarray-like data (continuous)
#' expr_data1 <- data.frame(g1=rnorm(4,10,2), g2=rnorm(4,8,1), g3=rnorm(4,12,3))
#' rownames(expr_data1) <- paste0("s", 1:4)
#' write.csv(t(expr_data1), "GSE001_expr.csv")
#' clin_data1 <- data.frame(row.names = paste0("s", 1:4), group = c("control","control","PD","PD"))
#' write.csv(clin_data1, "GSE001_pheno.csv")
#'
#' # Study 2: RNA-Seq-like data (discrete counts)
#' expr_data2 <- data.frame(g1=rnbinom(4,mu=100,size=10), g2=rnbinom(4,mu=50,size=10), g4=rnbinom(4,mu=200,size=10))
#' rownames(expr_data2) <- paste0("s", 5:8)
#' write.csv(t(expr_data2), "GSE002_expr.csv")
#' clin_data2 <- data.frame(row.names = paste0("s", 5:8), group = c("control","control","PD","PD"))
#' write.csv(clin_data2, "GSE002_pheno.csv")
#'
#' # --- Define file paths and run meta-analysis ---
#' data_paths <- c("GSE001_expr.csv", "GSE002_expr.csv")
#' clin_paths <- c("GSE001_pheno.csv", "GSE002_pheno.csv")
#'
#' # Run with auto-detection of data types
#' meta_results <- AWmeta(raw.data.dir = data_paths,
#'                        raw.clin.dir = clin_paths,
#'                        DE.method = c("limma", "edgeR"), # Provide one continuous, one discrete
#'                        compare.group = c("control", "PD"),
#'                        ref.level = "control",
#'                        paired = FALSE,
#'                        core.num = 4)
#'
#' # View the first few results
#' head(meta_results)
#' }


AWmeta <- function(raw.data.dir, raw.clin.dir, raw.sep, DE.method, compare.group, ref.level, paired, core.num) {
  
  # 1. Check for required packages
  required_packages <- c("MetaDE", "plyr", "AWFisher", "MetaVolcanoR", "parallel")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed. Please install it."), call. = FALSE)
    }
  }
  
  # 2. Input validation
  if (length(raw.data.dir) != length(raw.clin.dir)) {
    stop("The number of data files must match the number of clinical files.", call. = FALSE)
  }
  if (!(length(DE.method) %in% c(1, 2))) {
    stop("`DE.method` must have 1 or 2 elements.", call. = FALSE)
  }
  
  # 3. Generate internal variable names and load data
  study.index <- tools::file_path_sans_ext(list.files(raw.data.dir))
  study.clin.index <- tools::file_path_sans_ext(list.files(raw.clin.dir))
  raw.data.dir_full <- list.files(raw.data.dir, full.names = T)
  raw.clin.dir_full <- list.files(raw.clin.dir, full.names = T)
  num.studies <- length(study.index)
  
  message("Loading data...")
  data.list <- lapply(raw.data.dir_full, function(f) read.table(f, sep = raw.sep, header = TRUE, row.names = 1, check.names = FALSE))
  clin.data.list <- lapply(raw.clin.dir_full, function(f) read.table(f, sep = raw.sep, header = TRUE, row.names = 1, check.names = FALSE))
  names(data.list) <- study.index
  names(clin.data.list) <- study.clin.index
  
  # 4. Set up DE analysis parameters based on DE.method
  message("Setting up DE analysis parameters...")
  if (length(DE.method) == 1) {
    # Apply one method to all studies
    ind.method.index <- rep(DE.method, num.studies)
    message(paste("Applying ", DE.method, " to all studies.", sep=""))
  } else {
    # Auto-detect data type and assign method
    continuous_methods <- c("limma", "sam")
    discrete_methods <- c("edgeR", "DESeq2", "limmaVoom")
    
    # Identify which user-provided method is for which data type
    method1_is_cont <- DE.method[1] %in% continuous_methods
    method2_is_disc <- DE.method[2] %in% discrete_methods
    method1_is_disc <- DE.method[1] %in% discrete_methods
    method2_is_cont <- DE.method[2] %in% continuous_methods
    
    if (method1_is_cont && method2_is_disc) {
      cont_method_choice <- DE.method[1]
      disc_method_choice <- DE.method[2]
    } else if (method1_is_disc && method2_is_cont) {
      cont_method_choice <- DE.method[2]
      disc_method_choice <- DE.method[1]
    } else {
      stop("When DE.method has two elements, one must be for continuous data ('limma', 'sam') and one for discrete data ('edgeR', 'DESeq2', 'limmaVoom').", call. = FALSE)
    }
    
    ind.method.index <- character(num.studies) # pre-allocate
    message("Auto-detecting data types to assign DE methods:")
    for (i in 1:num.studies) {
      current_df_matrix <- as.matrix(data.list[[i]])
      # Check if all non-NA values are whole numbers
      is_discrete <- all(current_df_matrix[!is.na(current_df_matrix)] %% 1 == 0)
      
      if (is_discrete) {
        ind.method.index[i] <- disc_method_choice
        message(paste0("  - Study '", study.index[i], "': Detected discrete data."))
      } else {
        ind.method.index[i] <- cont_method_choice
        message(paste0("  - Study '", study.index[i], "': Detected continuous data."))
      }
    }
  }
  
  data.type.map <- c("limma" = "continuous", "sam" = "continuous", "edgeR" = "discrete", "DESeq2" = "discrete", "limmaVoom" = "discrete")
  data.type.index <- data.type.map[ind.method.index]
  
  if (any(is.na(data.type.index))) {
    stop("Invalid DE.method assigned. Check your `DE.method` parameter. Valid methods are 'limma', 'sam', 'edgeR', 'DESeq2', 'limmaVoom'.", call. = FALSE)
  }
  
  # 5. Perform individual DE analysis for each study
  message("Performing differential expression analysis for each study...")
  indi.de.results <- list()
  for (i in 1:num.studies) {
    message(paste("  - Analyzing study:", study.index[i]))
    indi.de.results[[i]] <- MetaDE::Indi.DE.Analysis(
      data = list(data.list[[i]]),
      clin.data = list(clin.data.list[[i]]),
      data.type = data.type.index[i],
      resp.type = "twoclass",
      ind.method = ind.method.index[i],
      select.group = compare.group,
      ref.level = ref.level,
      paired = paired
    )
  }
  names(indi.de.results) <- paste0("indi.de.", study.index)
  
  # 6. Aggregate all unique gene names
  gene.name.list <- lapply(indi.de.results, function(res) rownames(res$log2FC))
  gene.name.all <- data.frame(gene.name = sort(Reduce(union, gene.name.list)))
  
  # 7. Format p-values into a single data frame
  message("Formatting p-values for meta-analysis...")
  meta.p.df <- gene.name.all
  for (i in 1:num.studies) {
    res <- indi.de.results[[i]]
    temp.p <- data.frame(gene.name = rownames(res$p), p_value = res$p)
    meta.p.df <- plyr::join(x = meta.p.df, y = temp.p, by = "gene.name", type = "left")
  }
  rownames(meta.p.df) <- meta.p.df$gene.name
  meta.p.df <- meta.p.df[, -1, drop = FALSE]
  colnames(meta.p.df) <- paste0(study.index, "_Pvalue")
  
  # 8. Perform AWFisher meta-analysis
  message("Running AW-Fisher to combine p-values...")
  meta.pvalue.results <- t(apply(meta.p.df, 1, function(g) {
    valid.indices <- which(!is.na(g))
    if (length(valid.indices) > 1) {
      aw.fisher <- AWFisher::AWFisher_pvalue(g[valid.indices])
      aw.weights <- rep(NA, length(g))
      aw.weights[valid.indices] <- aw.fisher$weights
      return(c(aw.weights, aw.fisher$pvalues, 1)) # 1 is a placeholder for FDR calculation
    } else if (length(valid.indices) == 1) {
      aw.weights <- rep(NA, length(g))
      aw.weights[valid.indices] <- 1 # Weight is 1 for a single study
      return(c(aw.weights, g[valid.indices], 0)) # 0 means do not apply FDR correction later
    } else {
      return(c(rep(NA, length(g) + 1), 0)) # No data, no FDR
    }
  }))
  
  colnames(meta.pvalue.results) <- c(paste0(study.index, "_AWmeta_Weight"), "AWmeta_P_value", "FDR_flag")
  
  # 9. Calculate FDR for the meta-analysis p-values
  message("Calculating FDR...")
  FDR.index <- which(meta.pvalue.results[, "FDR_flag"] == 1)
  no.FDR.index <- which(meta.pvalue.results[, "FDR_flag"] == 0)
  
  # 10. Initialize FDR column
  meta.pvalue.results <- cbind(meta.pvalue.results, AWmeta_FDR = NA_real_)
  
  if (length(FDR.index) > 0) {
    meta.pvalue.results[FDR.index, "AWmeta_FDR"] <- p.adjust(meta.pvalue.results[FDR.index, "AWmeta_P_value"], method = "fdr")
  }
  if (length(no.FDR.index) > 0) {
    meta.pvalue.results[no.FDR.index, "AWmeta_FDR"] <- meta.pvalue.results[no.FDR.index, "AWmeta_P_value"]
  }
  meta.pvalue.results <- meta.pvalue.results[, -which(colnames(meta.pvalue.results) == "FDR_flag")] # Remove placeholder
  
  # 11. Format individual study fold-changes
  ind.FC.df <- gene.name.all
  for (i in 1:num.studies) {
    res <- indi.de.results[[i]]
    temp.fc <- data.frame(gene.name = rownames(res$log2FC), fc = res$log2FC)
    ind.FC.df <- plyr::join(x = ind.FC.df, y = temp.fc, by = "gene.name", type = "left")
  }
  rownames(ind.FC.df) <- ind.FC.df$gene.name
  ind.FC.df <- ind.FC.df[, -1, drop = FALSE]
  colnames(ind.FC.df) <- paste0(study.index, "_FC")
  
  # 12. Combine all intermediate results
  meta.aw.fisher.output <- cbind(ind.FC.df, meta.p.df, meta.pvalue.results)
  
  # 13. Prepare to calculate summary FC
  message("Preparing data for AWmeta fold-change calculation...")
  metaFC.input <- gene.name.all
  aw.weights <- meta.pvalue.results[, 1:num.studies, drop=FALSE]
  aw.weights[aw.weights == 0] <- NA # AWFisher can give 0 weights, which should be treated as no contribution
  
  col.name.vec <- "gene.name"
  for (i in 1:num.studies) {
    res <- indi.de.results[[i]]
    temp.data <- data.frame(
      gene.name = rownames(res$log2FC),
      log2fc = as.numeric(res$log2FC),
      vi = as.numeric(res$lfcSE)**2
    )
    metaFC.input <- plyr::join(x = metaFC.input, y = temp.data, by = "gene.name", type = "left")
    
    # Apply AW weights
    fc_col_idx <- ncol(metaFC.input) - 1
    vi_col_idx <- ncol(metaFC.input)
    metaFC.input[, fc_col_idx] <- metaFC.input[, fc_col_idx] * aw.weights[, i]
    metaFC.input[, vi_col_idx] <- metaFC.input[, vi_col_idx] * aw.weights[, i] # Note: weighting variance is non-standard but follows original script logic
    
    col.name.vec <- c(col.name.vec, paste0("Log2FC_", i), paste0("vi_", i))
  }
  colnames(metaFC.input) <- col.name.vec
  
  # 14. Calculate summary FC using a Random Effects Model in parallel
  message(paste("Running AW-REM to calculate AWmeta fold-change using", core.num, "CPU cores..."))

  metaFC.results <- do.call(rbind,
                            parallel::mclapply(split(metaFC.input, metaFC.input$gene.name),
                                               function(g) MetaVolcanoR::remodel(g, "Log2FC", "vi"),
                                               mc.cores = core.num))
  
  # 15. Finalize the output data frame
  metaFC.results$gene.name <- rownames(metaFC.results)
  final.metaFC <- plyr::join(gene.name.all, metaFC.results, by = "gene.name", type = "left")
  rownames(final.metaFC) <- final.metaFC$gene.name
  final.metaFC <- final.metaFC[, "randomSummary", drop = TRUE]
  
  # Combine everything for the final output
  final.results <- cbind(meta.aw.fisher.output, AWmeta_FC = final.metaFC)
  
  message("Sorting results by P-value...")
  # Use dplyr::arrange for clarity and robustness against NAs
  # NAs in FDR will be sorted to the end by default
  final.results <- dplyr::arrange(final.results, .data$AWmeta_P_value)
  
  message("Meta-analysis complete.")
  
  # 16. Return the final result
  return(final.results)
}
