#!/usr/bin/env Rscript
# Standalone script to test biology-batch confounding
# Tests association between biological variables and batch labels

#' Test biology-batch confounding
#' 
#' Analyzes the relationship between biological variables and batch labels
#' to detect potential confounding that could affect batch correction.
#' 
#' @param batch factor/character vector of batch labels
#' @param biology factor/character/numeric vector of biological labels or values
#' @param plot logical; create visualization plots
#' @return list with test results and statistics
#' @examples
#' \dontrun{
#' # Test with categorical biology
#' batch <- factor(rep(c("A", "B", "C"), each = 10))
#' biology <- factor(rep(c("Type1", "Type2"), 15))
#' result <- test_biology_batch_confounding(batch, biology)
#' 
#' # Test with continuous biology
#' biology_cont <- rnorm(30, mean = c(5, 6, 7)[as.numeric(batch)])
#' result <- test_biology_batch_confounding(batch, biology_cont)
#' }
test_biology_batch_confounding <- function(batch, biology, plot = TRUE) {
  batch <- as.factor(batch)
  n_batches <- nlevels(batch)
  n_samples <- length(batch)
  
  cat("═══════════════════════════════════════════════════════════════════════════════\n")
  cat("  Biology-Batch Confounding Analysis\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n")
  cat(sprintf("\nSamples: %d\n", n_samples))
  cat(sprintf("Batches: %d (%s)\n", n_batches, paste(levels(batch), collapse = ", ")))
  
  # Detect biology type
  if (is.numeric(biology)) {
    biology_type <- "continuous"
    cat(sprintf("Biology: Continuous variable\n"))
    cat(sprintf("  Mean: %.3f, SD: %.3f\n", mean(biology, na.rm = TRUE), sd(biology, na.rm = TRUE)))
  } else {
    biology <- as.factor(biology)
    biology_type <- "categorical"
    n_levels <- nlevels(biology)
    cat(sprintf("Biology: Categorical variable (%d levels: %s)\n", 
                n_levels, paste(levels(biology), collapse = ", ")))
  }
  
  # Test association
  cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("  Association Tests\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  
  if (biology_type == "categorical") {
    # Contingency table
    contingency <- table(batch, biology)
    cat("\nContingency Table:\n")
    print(contingency)
    
    # Chi-square test
    chi_test <- chisq.test(contingency)
    cat(sprintf("\nChi-square test of independence:\n"))
    cat(sprintf("  X-squared = %.4f\n", chi_test$statistic))
    cat(sprintf("  df = %d\n", chi_test$parameter))
    cat(sprintf("  p-value = %.4e\n", chi_test$p.value))
    
    # Cramér's V (effect size)
    n <- sum(contingency)
    cramers_v <- sqrt(chi_test$statistic / (n * (min(nrow(contingency), ncol(contingency)) - 1)))
    cat(sprintf("  Cramér's V = %.4f\n", cramers_v))
    
    # Interpretation
    if (chi_test$p.value < 0.001) {
      cat("\n  ⚠️  STRONG ASSOCIATION: Batch and biology are highly confounded!\n")
      cat("     Batch correction may remove biological signal.\n")
    } else if (chi_test$p.value < 0.05) {
      cat("\n  ⚠️  MODERATE ASSOCIATION: Some confounding detected.\n")
      cat("     Use anchor-based methods if available.\n")
    } else {
      cat("\n  ✓  NO SIGNIFICANT ASSOCIATION: Batch and biology appear independent.\n")
    }
    
    # Per-batch biology distribution
    cat("\nBiology distribution by batch:\n")
    for (b in levels(batch)) {
      idx <- batch == b
      bio_counts <- table(biology[idx])
      bio_props <- prop.table(bio_counts)
      cat(sprintf("  Batch %s (n=%d):\n", b, sum(idx)))
      for (bio in names(bio_counts)) {
        cat(sprintf("    %s: %d (%.1f%%)\n", bio, bio_counts[bio], bio_props[bio] * 100))
      }
    }
    
    test_result <- list(
      type = "categorical",
      contingency = contingency,
      chi_square = chi_test$statistic,
      p_value = chi_test$p.value,
      cramers_v = cramers_v,
      confounded = chi_test$p.value < 0.05
    )
    
  } else {
    # Continuous biology: ANOVA
    aov_fit <- aov(biology ~ batch)
    aov_summary <- summary(aov_fit)
    
    cat("\nANOVA: Biology ~ Batch\n")
    print(aov_summary[[1]])
    
    f_stat <- aov_summary[[1]]$`F value`[1]
    p_value <- aov_summary[[1]]$`Pr(>F)`[1]
    
    # Eta-squared (effect size)
    ss_between <- aov_summary[[1]]$`Sum Sq`[1]
    ss_total <- sum(aov_summary[[1]]$`Sum Sq`)
    eta_squared <- ss_between / ss_total
    
    cat(sprintf("\nEffect size:\n"))
    cat(sprintf("  Eta-squared = %.4f\n", eta_squared))
    
    # Per-batch statistics
    cat("\nBiology statistics by batch:\n")
    for (b in levels(batch)) {
      idx <- batch == b
      bio_vals <- biology[idx]
      cat(sprintf("  Batch %s (n=%d): Mean=%.3f, SD=%.3f\n", 
                  b, sum(idx), mean(bio_vals, na.rm = TRUE), sd(bio_vals, na.rm = TRUE)))
    }
    
    # Interpretation
    if (p_value < 0.001) {
      cat("\n  ⚠️  STRONG ASSOCIATION: Batch explains %.1f%% of biology variance!\n", eta_squared * 100)
      cat("     Batch correction may remove biological signal.\n")
    } else if (p_value < 0.05) {
      cat(sprintf("\n  ⚠️  MODERATE ASSOCIATION: Batch explains %.1f%% of biology variance.\n", eta_squared * 100))
      cat("     Use anchor-based methods if available.\n")
    } else {
      cat("\n  ✓  NO SIGNIFICANT ASSOCIATION: Batch and biology appear independent.\n")
    }
    
    test_result <- list(
      type = "continuous",
      f_statistic = f_stat,
      p_value = p_value,
      eta_squared = eta_squared,
      confounded = p_value < 0.05
    )
  }
  
  # Visualization
  if (plot) {
    cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
    cat("  Creating plots...\n")
    cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
    
    if (biology_type == "categorical") {
      # Bar plot: biology distribution by batch
      par(mfrow = c(1, 2))
      
      # Stacked bar plot
      barplot(t(contingency), 
              col = rainbow(ncol(contingency)),
              xlab = "Batch", ylab = "Count",
              main = "Biology Distribution by Batch",
              legend.text = colnames(contingency),
              args.legend = list(x = "topright"))
      
      # Grouped bar plot
      barplot(t(contingency), 
              beside = TRUE,
              col = rainbow(ncol(contingency)),
              xlab = "Batch", ylab = "Count",
              main = "Biology Distribution by Batch (Grouped)",
              legend.text = colnames(contingency),
              args.legend = list(x = "topright"))
      
    } else {
      # Box plots and violin plots
      par(mfrow = c(1, 2))
      
      # Box plot
      boxplot(biology ~ batch,
              xlab = "Batch", ylab = "Biology Value",
              main = "Biology Distribution by Batch",
              col = rainbow(n_batches))
      
      # Stripchart
      stripchart(biology ~ batch,
                 method = "jitter",
                 vertical = TRUE,
                 xlab = "Batch", ylab = "Biology Value",
                 main = "Biology Values by Batch",
                 col = rainbow(n_batches),
                 pch = 19)
    }
    
    par(mfrow = c(1, 1))
  }
  
  cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
  cat("  Analysis Complete\n")
  cat("═══════════════════════════════════════════════════════════════════════════════\n")
  
  invisible(test_result)
}

# Example usage (commented out)
if (FALSE) {
  # Example 1: Confounded case (categorical)
  set.seed(123)
  batch <- factor(rep(c("A", "B", "C"), each = 20))
  biology <- factor(c(rep("Type1", 30), rep("Type2", 30)))  # Confounded!
  result <- test_biology_batch_confounding(batch, biology)
  
  # Example 2: Independent case (categorical)
  batch <- factor(rep(c("A", "B", "C"), each = 20))
  biology <- factor(sample(c("Type1", "Type2"), 60, replace = TRUE))  # Independent
  result <- test_biology_batch_confounding(batch, biology)
  
  # Example 3: Confounded case (continuous)
  batch <- factor(rep(c("A", "B", "C"), each = 20))
  biology <- rnorm(60, mean = c(5, 6, 7)[as.numeric(batch)])  # Confounded!
  result <- test_biology_batch_confounding(batch, biology)
  
  # Example 4: Independent case (continuous)
  batch <- factor(rep(c("A", "B", "C"), each = 20))
  biology <- rnorm(60, mean = 6)  # Independent
  result <- test_biology_batch_confounding(batch, biology)
}


