library("lmerTest")
library("readr")

data_dir <- "/Users/chloehampson/Desktop/projects/abide-analysis/dset/group/habenula/pheno/" # Make sure to leave the slash at the end
clusters <- c("1", "1a", "2", "2a")
roi <- "RSFC"

# Level-1 Predictors
group_var <- "Group"
categorical_vars <- c("Sex")
numerical_vars <- c("Age")
phen_vars <- c("Phen1", "Phen2", "Phen3", "Phen4")

for (cluster in clusters) {
    data_path <- paste0(data_dir, "cluster-", cluster, "_data.csv")
    data <- read.table(file = data_path, sep = ",", header = TRUE)

    for (phen_var in phen_vars) {
        all_columns <- c(roi, categorical_vars, numerical_vars, phen_var, group_var, "Site")
        sub_data <- data[, all_columns]
        sub_data <- na.omit(sub_data)

        # Convert categorical variables to factors
        for (var in categorical_vars) {
            sub_data[[var]] <- factor(sub_data[[var]])
        }

        # Relevel 'Group' to make NT the reference category
        sub_data[[group_var]] <- relevel(factor(sub_data[[group_var]]), ref = "NT")

        # Scale continuous predictors
        for (var in numerical_vars) {
            sub_data[[var]] <- scale(sub_data[[var]], center = TRUE, scale = TRUE)
        }

        sub_data[[phen_var]] <- scale(sub_data[[phen_var]], center = TRUE, scale = TRUE)

        # Fixed effects formula
        fixed_effects <- paste(c(numerical_vars, categorical_vars), collapse = " + ")

        # Full model with interaction
        equation_lme <- paste(roi, "~", fixed_effects, "+", group_var, "*", phen_var, "+ (1|Site)")

        # Run the full model
        model <- lmer(as.formula(equation_lme), data = sub_data)
        model_summary <- summary(model)
        print(model_summary)

        # Save results
        model_table <- as.data.frame(coef(summary(model)))
        out_file <- paste0(data_dir, "cluster-", cluster, "_", phen_var, "_table.csv")
        write.csv(model_table, file = out_file, row.names = TRUE)

        # ASD-Only Model (Main Effects Only, No Interaction)
        asd_data <- sub_data[sub_data[[group_var]] == "ASD", ]
        equation_lme_asd <- paste(roi, "~", paste(c(fixed_effects, phen_var), collapse = " + "), "+ (1|Site)")

        # Run ASD model
        model_asd <- lmer(as.formula(equation_lme_asd), data = asd_data)
        model_asd_summary <- summary(model_asd)
        print(model_asd_summary)

        # Save ASD-only results
        asd_phen_effect <- as.data.frame(coef(summary(model_asd)))
        out_file_asd <- paste0(data_dir, "cluster-", cluster, "_", phen_var, "_ASD_only.csv")
        write.csv(asd_phen_effect, file = out_file_asd, row.names = TRUE)
    }
}
