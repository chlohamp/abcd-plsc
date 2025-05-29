library("lmerTest")
library("readr")
library("ggplot2")
library("dplyr")

data_dir <- "/Users/chloehampson/Desktop/projects/abcd-plsc/derivatives/none-reduced/regression/"
dimensions <- c("dim1", "dim3")
score <- "score"

# Only plot for these variables
significant_vars <- c("resp_wheeze_yn_y")

# Level-1 Predictors
categorical_vars <- c("demo_sex_v2", "demo_prnt_gender_id_v2", "demo_origin_v2", "mri_info_manufacturer")
numerical_vars <- c("interview_age", "demo_prnt_age_v2", "demo_prnt_ed_v2_2yr_l", "demo_prtnr_ed_v2_2yr_l", "demo_comb_income_v2", "rsfmri_meanmotion")
phyhealth_vars <- c(
    "mctq_sdweek_calc", "mctq_msfsc_calc", "resp_wheeze_yn_y", "resp_pmcough_yn_y",
    "resp_diagnosis_yn_y", "resp_bronch_yn_y", "blood_pressure_sys_mean", "blood_pressure_dia_mean",
    "physical_activity1_y", "cbcl_scr_syn_internal_t", "cbcl_scr_syn_external_t"
)

phyhealth_cats <- c("resp_wheeze_yn_y", "resp_pmcough_yn_y", "resp_diagnosis_yn_y", "resp_bronch_yn_y")

for (dimension in dimensions) {
    data_path <- paste0(data_dir, "phyhealth_", dimension, "_data.csv")
    data <- read.table(file = data_path, sep = ",", header = TRUE)

    for (phyhealth_var in phyhealth_vars) {
        all_columns <- c(score, categorical_vars, numerical_vars, phyhealth_var, "site_id_l", "rel_family_id")
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

        # Run the model
        model <- lmer(as.formula(equation_lme), data = sub_data)
        model_summary <- summary(model)
        print(model_summary)

        # Save results to CSV
        model_table <- as.data.frame(coef(summary(model)))
        out_file <- paste0(data_dir, "phyhealth_", dimension, "_", phyhealth_var, "_table.csv")
        write.csv(model_table, file = out_file, row.names = TRUE)

        # Only plot if it's a significant variable and categorical
        if (phyhealth_var %in% significant_vars) {
            sub_data$predicted <- predict(model)

            p <- ggplot(sub_data, aes_string(x = phyhealth_var, y = "predicted")) +
                geom_jitter(width = 0.2, alpha = 0.4) +
                stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "blue") +
                stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
                labs(
                    x = phyhealth_var, y = "Predicted Score",
                    title = paste("Predicted", score, "by", phyhealth_var, "(", dimension, ")")
                ) +
                theme_classic(base_size = 10)

            plot_file <- paste0(data_dir, "plot_phyhealth_", dimension, "_", phyhealth_var, ".png")
            ggsave(plot_file, plot = p, width = 6, height = 4, dpi = 300)
        }
    }
}
