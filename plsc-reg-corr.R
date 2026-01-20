library("lmerTest")
library("readr")

# Use sum-to-zero contrasts so Type III ANOVA p-values are meaningful for factors
options(contrasts = c("contr.sum", "contr.poly"))

data_dir <- "/Users/chloehampson/Desktop/abcd-plsc/derivatives/none-reduced/regression/dim1/" # Make sure to leave the slash at the end
networks <- c("cgc-dt", "dt-dla", "dt-dt", "dt-vs", "vs-vs")
#networks <- c("dt-smm", "vta-vs")
roi <- "rsfc"

# Level-1 Predictors
categorical_vars <- c("demo_sex_v2", "demo_prnt_gender_id_v2", "demo_origin_v2", "mri_info_manufacturer")
numerical_vars <- c("interview_age", "demo_prnt_age_v2", "demo_prnt_ed_v2_2yr_l", "demo_prtnr_ed_v2_2yr_l", "demo_comb_income_v2", "rsfmri_meanmotion")
phyhealth_vars <- c(
    "BMI",
    "mctq_sdweek_calc",
    "mctq_msfsc_calc",
    "resp_wheeze_yn_y",
    "resp_pmcough_yn_y",
    "resp_diagnosis_yn_y",
    "resp_bronch_yn_y",
    "blood_pressure_sys_mean",
    "blood_pressure_dia_mean",
    "physical_activity1_y",
    "cbcl_scr_syn_internal_t",
    "cbcl_scr_syn_external_t"
)

phyhealth_cats <- c("resp_wheeze_yn_y", "resp_pmcough_yn_y", "resp_diagnosis_yn_y", "resp_bronch_yn_y")

# Collector for p-values across all networks and phyhealth variables
results <- data.frame(
    network = character(),
    phyhealth_var = character(),
    var_type = character(),
    N = integer(),
    p_value = numeric(),
    stringsAsFactors = FALSE
)

for (network in networks) {
    data_path <- paste0(data_dir, "phyhealth_", network, "_data.csv")
    data <- read.table(file = data_path, sep = ",", header = TRUE)

    for (phyhealth_var in phyhealth_vars) {
        all_columns <- c(roi, categorical_vars, numerical_vars, phyhealth_var, "site_id_l", "rel_family_id")
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
            sub_data[[phyhealth_var]] <- scale(sub_data[[phyhealth_var]],
                center = TRUE, scale = TRUE
            )
        }

        # Fixed effects formula
        fixed_effects <- paste(c(numerical_vars, categorical_vars), collapse = " + ")

        # Full model
        equation_lme <- paste(roi, "~", phyhealth_var, "+", fixed_effects, "+ (1|site_id_l/rel_family_id)") # "+ (1|site_id_l/rel_family_id)")

        # Run the full model (guard against failures on edge cases)
        fit_ok <- TRUE
        model <- tryCatch({
            lmer(as.formula(equation_lme), data = sub_data)
        }, error = function(e) {
            message(sprintf("Model failed for %s in %s: %s", phyhealth_var, network, e$message))
            fit_ok <<- FALSE
            return(NULL)
        })

        if (fit_ok && !is.null(model)) {
            model_summary <- summary(model)
            print(model_summary)
        }

        if (fit_ok && !is.null(model)) {
            # Create formula for plotting
            plot_formula <- as.formula(paste(roi, "~ fitted(.) |", phyhealth_var))
            # Generate and print plot
            p <- plot(model, plot_formula, abline = c(0, 1))
            print(p)
        }

        # Save per-model coefficient table if model fit
        if (fit_ok && !is.null(model)) {
            model_table <- as.data.frame(coef(summary(model)))
            out_file <- paste0(data_dir, "phyhealth_", network, "_", phyhealth_var, "_table.csv")
            write.csv(model_table, file = out_file, row.names = TRUE)
        }

        # Extract a single p-value for the phyhealth_var using Type III ANOVA (overall effect for factors)
        p_val <- NA_real_
        if (fit_ok && !is.null(model)) {
            an_tab <- tryCatch({
                anova(model, type = 3)
            }, error = function(e) NULL)

            if (!is.null(an_tab)) {
                # Find the p-value column (usually "Pr(>F)")
                p_col <- grep("^Pr\\(>F\\)$", colnames(an_tab), value = TRUE)
                if (length(p_col) == 1 && phyhealth_var %in% rownames(an_tab)) {
                    p_val <- suppressWarnings(as.numeric(an_tab[phyhealth_var, p_col]))
                } else {
                    # Fallback for continuous variable: take coefficient p-value directly
                    cs <- tryCatch(as.data.frame(coef(summary(model))), error = function(e) NULL)
                    if (!is.null(cs) && phyhealth_var %in% rownames(cs)) {
                        if ("Pr(>|t|)" %in% colnames(cs)) {
                            p_val <- suppressWarnings(as.numeric(cs[phyhealth_var, "Pr(>|t|)"]))
                        }
                    }
                }
            }
        }

        # Record result row
        results <- rbind(
            results,
            data.frame(
                network = network,
                phyhealth_var = phyhealth_var,
                var_type = ifelse(phyhealth_var %in% phyhealth_cats, "categorical", "continuous"),
                N = nrow(sub_data),
                p_value = p_val,
                stringsAsFactors = FALSE
            )
        )
    }
}

# Add significance flag and write summary table to the regression folder (parent of data_dir)
if (!is.null(results) && nrow(results) > 0) {
    results$significant_p05 <- ifelse(!is.na(results$p_value) & results$p_value < 0.05, TRUE, FALSE)
    results$significant_p01 <- ifelse(!is.na(results$p_value) & results$p_value < 0.01, TRUE, FALSE)
    reg_root <- dirname(data_dir)  # parent directory of the dimension folder
    summary_out <- file.path(reg_root, "plsc-reg-corr-dim1-results.csv")
    write.csv(results, file = summary_out, row.names = FALSE)
    message(sprintf("Wrote summary p-values table to regression folder: %s", summary_out))
}
