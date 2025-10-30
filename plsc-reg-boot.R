library("lmerTest")
library("readr")
library("ggplot2")
library("dplyr")

# Use sum-to-zero contrasts so Type III ANOVA p-values are meaningful for factors
options(contrasts = c("contr.sum", "contr.poly"))

# Base directory for regression analyses
reg_dir <- "/Users/chloehampson/Desktop/abcd-plsc/derivatives/none-reduced/regression"
dimensions <- c("dim1", "dim3")
score <- "score"  # This is our bootstrap ratio outcome variable

# Level-1 Predictors
categorical_vars <- c("demo_sex_v2", "demo_prnt_gender_id_v2", "demo_origin_v2", "mri_info_manufacturer")
numerical_vars <- c("interview_age", "demo_prnt_age_v2", "demo_prnt_ed_v2_2yr_l", "demo_prtnr_ed_v2_2yr_l", "demo_comb_income_v2", "rsfmri_meanmotion")
phyhealth_vars <- c(
    "BMI", "mctq_sdweek_calc", "mctq_msfsc_calc", "resp_wheeze_yn_y", "resp_pmcough_yn_y",
    "resp_diagnosis_yn_y", "resp_bronch_yn_y", "blood_pressure_sys_mean", "blood_pressure_dia_mean",
    "physical_activity1_y", "cbcl_scr_syn_internal_t", "cbcl_scr_syn_external_t"
)

phyhealth_cats <- c("resp_wheeze_yn_y", "resp_pmcough_yn_y", "resp_diagnosis_yn_y", "resp_bronch_yn_y")

for (dim in dimensions) {
    message(sprintf("\nProcessing %s...", dim))
    
    # Set up results collector for this dimension
    results <- data.frame(
        phyhealth_var = character(),
        var_type = character(),
        N = integer(),
        p_value = numeric(),
        stringsAsFactors = FALSE
    )
    
    # Load the bootstrap ratio data for this dimension
    data_path <- file.path(reg_dir, dim, sprintf("phyhealth_%s_data.csv", dim))
    if (!file.exists(data_path)) {
        message(sprintf("Data file not found for %s: %s", dim, data_path))
        next
    }
    
    data <- read.table(file = data_path, sep = ",", header = TRUE)
    
    for (phyhealth_var in phyhealth_vars) {
        message(sprintf("  Analyzing %s...", phyhealth_var))
        
        all_columns <- c(score, categorical_vars, numerical_vars, phyhealth_var, "site_id_l", "rel_family_id")
        # Check if all required columns exist
        missing_cols <- setdiff(all_columns, colnames(data))
        if (length(missing_cols) > 0) {
            message(sprintf("    Missing columns for %s: %s", phyhealth_var, paste(missing_cols, collapse = ", ")))
            next
        }
        
        sub_data <- data[, all_columns]
        sub_data <- na.omit(sub_data)
        
        # Convert categorical variables to factors
        for (var in categorical_vars) {
            sub_data[[var]] <- factor(sub_data[[var]])
        }
        
        # Scale continuous predictors
        for (var in numerical_vars) {
            sub_data[[var]] <- scale(sub_data[[var]], center = TRUE, scale = TRUE)
        }
        
        if (phyhealth_var %in% phyhealth_cats) {
            sub_data[[phyhealth_var]] <- factor(sub_data[[phyhealth_var]])
        } else {
            sub_data[[phyhealth_var]] <- scale(sub_data[[phyhealth_var]], center = TRUE, scale = TRUE)
        }
        
        # Fixed effects formula
        fixed_effects <- paste(c(numerical_vars, categorical_vars), collapse = " + ")
        equation_lme <- paste(score, "~", phyhealth_var, "+", fixed_effects, "+ (1|site_id_l/rel_family_id)")
        
        # Run the model (guard against failures)
        fit_ok <- TRUE
        model <- tryCatch({
            lmer(as.formula(equation_lme), data = sub_data)
        }, error = function(e) {
            message(sprintf("    Model failed: %s", e$message))
            fit_ok <<- FALSE
            return(NULL)
        })
        
        if (fit_ok && !is.null(model)) {
            model_summary <- summary(model)
            print(model_summary)
            
            # Save per-model coefficient table
            model_table <- as.data.frame(coef(summary(model)))
            out_file <- file.path(reg_dir, dim, sprintf("phyhealth_%s_%s_table.csv", dim, phyhealth_var))
            write.csv(model_table, file = out_file, row.names = TRUE)
            
            # Extract p-value using Type III ANOVA
            p_val <- NA_real_
            an_tab <- tryCatch({ anova(model, type = 3) }, error = function(e) NULL)
            
            if (!is.null(an_tab)) {
                # Find the p-value column (usually "Pr(>F)")
                p_col <- grep("^Pr\\(>F\\)$", colnames(an_tab), value = TRUE)
                if (length(p_col) == 1 && phyhealth_var %in% rownames(an_tab)) {
                    p_val <- suppressWarnings(as.numeric(an_tab[phyhealth_var, p_col]))
                } else {
                    # Fallback for continuous variable: take coefficient p-value directly
                    cs <- model_table
                    if (phyhealth_var %in% rownames(cs)) {
                        if ("Pr(>|t|)" %in% colnames(cs)) {
                            p_val <- suppressWarnings(as.numeric(cs[phyhealth_var, "Pr(>|t|)"]))
                        }
                    }
                }
            }
            
            # Record result
            results <- rbind(
                results,
                data.frame(
                    phyhealth_var = phyhealth_var,
                    var_type = ifelse(phyhealth_var %in% phyhealth_cats, "categorical", "continuous"),
                    N = nrow(sub_data),
                    p_value = p_val,
                    stringsAsFactors = FALSE
                )
            )
        }
    }
    
    # Add significance flags and write dimension summary
    if (!is.null(results) && nrow(results) > 0) {
        results$significant_p05 <- ifelse(!is.na(results$p_value) & results$p_value < 0.05, TRUE, FALSE)
        results$significant_p01 <- ifelse(!is.na(results$p_value) & results$p_value < 0.01, TRUE, FALSE)
        
        summary_out <- file.path(reg_dir, sprintf("plsc-reg-boot-%s-results.csv", dim))
        write.csv(results, file = summary_out, row.names = FALSE)
        message(sprintf("Wrote summary p-values table: %s", summary_out))
    } else {
        message(sprintf("No results to write for %s", dim))
    }
} # End of for (dim in dimensions) loop