# server_gmm.R
# This module contains the logic for the "Subpopulation Detection (GMM)" tab.

# =========================================================================
# UTILITY FUNCTIONS FOR GMM ANALYSIS
# =========================================================================

# Standardizes a numeric vector to have a mean of 0 and a standard deviation of 1
z_transform <- function(x) {
  if (sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  }
  return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

# Applies a Yeo-Johnson transformation if the data is highly skewed
apply_conditional_yeo_johnson <- function(data_vector, skewness_threshold = 0.5) {
  transformed_data <- data_vector
  transformation_applied <- FALSE
  skew <- moments::skewness(data_vector, na.rm = TRUE)

  if (abs(skew) > skewness_threshold) {
    tryCatch({
      pt_result <- car::powerTransform(data_vector)
      lambda <- pt_result$lambda
      transformed_data <- car::yjPower(data_vector, lambda)
      transformation_applied <- TRUE
      message(paste("Yeo-Johnson transformation applied (skewness=", round(skew, 2), ")"))
    }, error = function(e) {
      warning(paste("Could not apply Yeo-Johnson transformation:", e$message))
    })
  } else {
    message(paste("Yeo-Johnson transformation not needed (skewness=", round(skew, 2), ")"))
  }

  return(list(transformed_data = transformed_data, transformation_applied = transformation_applied))
}

# Runs the GMM analysis using the mclust package
run_gmm_with_criterion <- function(data_mat, G_range, modelNames) {
  if (!is.matrix(data_mat) && !is.data.frame(data_mat)) {
    stop("Input data_mat must be a matrix or data frame for GMM analysis.")
  }
  if (!all(sapply(data_mat, is.numeric))) {
    stop("All columns in data_mat must be numeric.")
  }
  if (any(is.na(data_mat))) {
    stop("Input data_mat contains NA values. Please remove or impute before clustering.")
  }

  tryCatch({
    gmm_model <- mclust::Mclust(data_mat, G = G_range, modelNames = modelNames)
    return(gmm_model)
  }, error = function(e) {
    warning(paste("GMM run failed:", e$message))
    return(NULL)
  })
}

# Assigns cluster labels to the original data frame
assign_clusters <- function(df, gmm_model) {
  if (is.null(gmm_model) || is.null(gmm_model$classification)) {
    warning("GMM model or classification is NULL. Cannot assign clusters.")
    return(df)
  }
  df$cluster <- gmm_model$classification
  return(df)
}

# Generates a scatter plot of age vs. value, colored by cluster
plot_value_age <- function(df, value_col_name, age_col_name) {
  if (is.null(df) || nrow(df) == 0) {
    return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No GMM data available for plotting.", size = 6, color = "grey50"))
  }

  plot_title <- paste(value_col_name, "vs", age_col_name, "by Subpopulation Cluster")
  
  cluster_means <- df %>%
    dplyr::group_by(Gender, cluster) %>%
    dplyr::summarise(mean_Age = mean(Age, na.rm = TRUE),
                     mean_Value = mean(Value, na.rm = TRUE),
                     .groups = 'drop')

  ggplot2::ggplot(df, ggplot2::aes(x = Age, y = Value, color = factor(cluster))) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.2, height = 0.2), alpha = 0.6) +
    ggplot2::stat_ellipse(geom = "polygon", ggplot2::aes(fill = factor(cluster)), alpha = 0.2, show.legend = FALSE, level = 0.95) +
    ggplot2::geom_point(data = cluster_means, ggplot2::aes(x = mean_Age, y = mean_Value), shape = 4, size = 5, color = "red", stroke = 2) +
    ggplot2::facet_wrap(~Gender, labeller = as_labeller(function(x) paste(x, "Population"))) +
    ggplot2::labs(title = plot_title,
                  x = age_col_name, y = value_col_name, color = "Cluster") +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::theme_light() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 20)),
      strip.text = ggplot2::element_text(size = 14, face = "bold", color = "black"),
      strip.background = ggplot2::element_rect(fill = "#EEEEEE", color = NA),
      axis.title.x = ggplot2::element_text(size = 14, margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_text(size = 14, margin = ggplot2::margin(r = 10)),
      axis.text = ggplot2::element_text(size = 11, color = "black"),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 10, color = "black"),
      legend.position = "bottom",
      legend.background = ggplot2::element_rect(fill = "white", color = "grey90", size = 0.5, linetype = "solid"),
      panel.border = ggplot2::element_rect(colour = "#CCCCCC", fill = NA, size = 1),
      panel.grid.major = ggplot2::element_line(color = "#F0F0F0", size = 0.5),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# Runs the GMM analysis on a specific data subset
run_gmm_analysis_on_subset <- function(data_subset, gender_label, value_col_name, age_col_name, message_rv, progress_increment, model_names, g_range) {
    if (nrow(data_subset) > 0) {
        yj_result <- apply_conditional_yeo_johnson(data_subset$Value)
        data_subset$Value_transformed <- yj_result$transformed_data
        
        data_subset$Value_z <- z_transform(data_subset$Value_transformed)
        data_subset$Age_z <- z_transform(data_subset$Age)
        
        incProgress(progress_increment, detail = paste("Running GMM for", gender_label, "data (BIC)..."))
        
        gmm_model <- tryCatch({
            run_gmm_with_criterion(data_subset %>% dplyr::select(Value = Value_z, Age = Age_z), G_range = g_range, modelNames = model_names)
        }, error = function(e) {
            message_rv(list(text = paste("Error running BIC GMM for", gender_label, "data:", e$message), type = "error"))
            NULL
        })

        if (is.null(gmm_model)) {
            stop(paste("GMM model for", gender_label, "data could not be generated."))
        }
        
        data_clustered <- assign_clusters(data_subset, gmm_model)
        data_clustered$cluster <- as.factor(data_clustered$cluster)

        return(list(model = gmm_model, clustered_data = data_clustered, transformed_flag = yj_result$transformation_applied))
    } else {
        message_rv(list(text = paste("Error: No", tolower(gender_label), "data found after filtering. Please check the gender column and selection."), type = "error"))
        return(list(model = NULL, clustered_data = NULL, transformed_flag = FALSE))
    }
}

# =========================================================================
# MAIN SERVER LOGIC
# =========================================================================

gmmServer <- function(input, output, session, gmm_uploaded_data_rv, gmm_processed_data_rv, gmm_transformation_details_rv, message_rv, analysis_running_rv) {

  # Reactive value to hold models for BIC criterion
  gmm_models_bic_rv <- reactiveVal(list(male = NULL, female = NULL))

  # Helper function to guess column names
  guess_column <- function(cols_available, common_names) {
    for (name in common_names) {
      match_idx <- grep(paste0("^", name, "$"), cols_available, ignore.case = TRUE)
      if (length(match_idx) > 0) {
        return(cols_available[match_idx[1]])
      }
    }
    return("")
  }

  # Observer for GMM file upload
  observeEvent(input$gmm_file_upload, {
    req(input$gmm_file_upload)
    tryCatch({
      data <- readxl::read_excel(input$gmm_file_upload$datapath)
      gmm_uploaded_data_rv(data)
      message_rv(list(text = "GMM data uploaded successfully.", type = "success"))

      col_names <- colnames(data)
      all_col_choices_with_none <- c("None" = "", col_names)

      updateSelectInput(session, "gmm_value_col", choices = all_col_choices_with_none, selected = guess_column(col_names, c("Value", "Result", "Measurement", "Waarde", "HGB", "hgb", "HB", "hb")))
      updateSelectInput(session, "gmm_age_col", choices = all_col_choices_with_none, selected = guess_column(col_names, c("Age", "age", "leeftijd")))
      updateSelectInput(session, "gmm_gender_col", choices = all_col_choices_with_none, selected = guess_column(col_names, c("Gender", "gender", "Sex", "sex", "geslacht")))

    }, error = function(e) {
      message_rv(list(text = paste("Error reading GMM file:", e$message), type = "error"))
      gmm_uploaded_data_rv(NULL)
    })
  })

  # Renders the gender choice radio buttons if a gender column is selected
  output$gmm_gender_choice_ui <- renderUI({
    req(input$gmm_gender_col)
    if (input$gmm_gender_col != "") {
      radioButtons(inputId = "gmm_gender_choice", label = "Select Gender Analysis:", choices = c("Male" = "Male", "Female" = "Female", "Both" = "Both"), selected = "Both", inline = TRUE)
    }
  })

  # Renders the manual model selection UI
  output$gmm_manual_model_ui <- renderUI({
    if (input$gmm_model_selection_choice == "Manual Selection") {
      model_names <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EVE", "VEE", "VVV")
      
      tagList(
        selectInput("gmm_manual_model", "Select GMM Model:", choices = model_names, selected = "EII"),
        selectInput("gmm_component_number", "Select Component Number:", choices = c("Auto", as.character(2:10)), selected = "Auto")
      )
    }
  })

  # Core observer for the GMM analysis button
  observeEvent(input$run_gmm_analysis_btn, {
    if (is.null(gmm_uploaded_data_rv())) {
      message_rv(list(text = "Please upload an Excel file first.", type = "error"))
      return(NULL)
    }
    
    if (input$gmm_value_col == "" || input$gmm_age_col == "") {
      message_rv(list(text = "Please select the columns from the dropdown menus.", type = "error"))
      return(NULL)
    }

    data_check <- gmm_uploaded_data_rv()

    if (!(input$gmm_value_col %in% colnames(data_check))) {
      message_rv(list(text = paste0("Error: The selected column for Values ('", input$gmm_value_col, "') was not found in your uploaded data. Please select a valid column."), type = "error"))
      return(NULL)
    }
    if (!(input$gmm_age_col %in% colnames(data_check))) {
      message_rv(list(text = paste0("Error: The selected column for Age ('", input$gmm_age_col, "') was not found in your uploaded data. Please select a valid column."), type = "error"))
      return(NULL)
    }
    if (input$gmm_gender_col != "" && !(input$gmm_gender_col %in% colnames(data_check))) {
      message_rv(list(text = paste0("Error: The selected gender column '", input$gmm_gender_col, "' was not found in your uploaded data. Please select a valid column."), type = "error"))
      return(NULL)
    }

    if (input$gmm_value_col == input$gmm_age_col) {
      message_rv(list(text = "Error: The same column cannot be selected for both Values and Age. Please choose a different column for one of the inputs.", type = "error"))
      return(NULL)
    }

    if (input$gmm_gender_col != "") {
      req(input$gmm_gender_choice)
    }

    if (input$gmm_model_selection_choice == "Manual Selection") {
      req(input$gmm_manual_model)
      model_names_valid <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EVE", "VEE", "VVV")
      if (!input$gmm_manual_model %in% model_names_valid) {
        message_rv(list(text = paste0("Error: Invalid GMM model name selected. Please choose from the list."), type = "error"))
        return(NULL)
      }
    }


    if (analysis_running_rv()) {
      message_rv(list(text = "An analysis is already running. Please wait.", type = "warning"))
      return(NULL)
    }

    gmm_processed_data_rv(NULL)
    gmm_models_bic_rv(list(male=NULL, female=NULL))
    message_rv(list(text = "Starting new GMM analysis...", type = "info"))


    analysis_running_rv(TRUE)
    shinyjs::disable("run_gmm_analysis_btn")
    shinyjs::runjs("$('#run_gmm_analysis_btn').text('Analyzing...');")
    session$sendCustomMessage('analysisStatus', TRUE)

    model_names_to_use <- if (input$gmm_model_selection_choice == "Auto-select") {
      c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EVE", "VEE", "VVV")
    } else {
      input$gmm_manual_model
    }

    g_range_to_use <- if (input$gmm_model_selection_choice == "Auto-select") {
      2:10
    } else {
      if (input$gmm_component_number == "Auto") {
        2:10
      } else {
        2:as.numeric(input$gmm_component_number)
      }
    }

    tryCatch({
      withProgress(message = 'Running GMM Analysis', value = 0, {
        incProgress(0.1, detail = "Loading data...")
  
        data <- gmm_uploaded_data_rv()
        value_col <- input$gmm_value_col
        age_col <- input$gmm_age_col
        gender_col <- input$gmm_gender_col
        gender_choice <- if (gender_col == "") "None" else input$gmm_gender_choice
        
        value_col_sym <- rlang::sym(value_col)
        age_col_sym <- rlang::sym(age_col)
  
        if (gender_col != "") {
          gender_col_sym <- rlang::sym(gender_col)
          gmm_data <- data %>%
            dplyr::select(Value = !!value_col_sym, Age = !!age_col_sym, Gender_orig = !!gender_col_sym)
        } else {
          gmm_data <- data %>%
            dplyr::select(Value = !!value_col_sym, Age = !!age_col_sym)
        }
        
        original_rows_count <- nrow(gmm_data)
        
        gmm_data <- gmm_data %>%
          mutate(
            Value = as.numeric(as.character(Value)),
            Age = as.numeric(as.character(Age))
          ) %>%
          na.omit()
        
        removed_rows_count <- original_rows_count - nrow(gmm_data)
  
        if (nrow(gmm_data) == 0) {
          message_rv(list(text = "No complete rows for GMM after data cleaning and NA removal. Check your data or selections.", type = "error"))
          return(NULL)
        }
  
        incProgress(0.2, detail = "Splitting data by gender and transforming...")
  
        if (gender_col != "") {
          gmm_data <- gmm_data %>%
            mutate(Gender = case_when(
              grepl("male|m|man|jongen(s)?|heren|mannelijk(e)?", Gender_orig, ignore.case = TRUE) ~ "Male",
              grepl("female|f|vrouw(en)?|v|meisje(s)?|dame|mevr|vrouwelijke", Gender_orig, ignore.case = TRUE) ~ "Female",
              TRUE ~ "Other"
            )) %>%
            filter(Gender %in% c("Male", "Female"))
        } else {
          gmm_data <- gmm_data %>%
            mutate(Gender = "Combined")
        }
        
        if (nrow(gmm_data) == 0) {
          message_rv(list(text = "Filtered dataset is empty after gender selection. Please check the data or gender column.", type = "warning"))
          return(NULL)
        }
  
        combined_clustered_data <- tibble()
        male_value_transformed_flag <- FALSE
        female_value_transformed_flag <- FALSE
        combined_gmm_model_bic <- NULL
        male_gmm_model_bic <- NULL
        female_gmm_model_bic <- NULL
        
        results <- list(male = NULL, female = NULL, combined = NULL)
  
        if (gender_col == "" || gender_choice == "Both") {
          if (gender_col == "") {
            message("Running GMM on combined data...")
            results$combined <- run_gmm_analysis_on_subset(
              data_subset = gmm_data,
              gender_label = "Combined",
              value_col_name = value_col,
              age_col_name = age_col,
              message_rv = message_rv,
              progress_increment = 0.6,
              model_names = model_names_to_use,
              g_range = g_range_to_use
            )
            if (!is.null(results$combined$model)) {
              combined_gmm_model_bic <- results$combined$model
              male_value_transformed_flag <- results$combined$transformed_flag
              combined_clustered_data <- dplyr::bind_rows(combined_clustered_data, results$combined$clustered_data)
            }
          } else {
            male_data <- gmm_data %>% dplyr::filter(Gender == "Male")
            results$male <- run_gmm_analysis_on_subset(
              data_subset = male_data,
              gender_label = "Male",
              value_col_name = value_col,
              age_col_name = age_col,
              message_rv = message_rv,
              progress_increment = 0.3,
              model_names = model_names_to_use,
              g_range = g_range_to_use
            )
            if (!is.null(results$male$model)) {
              male_gmm_model_bic <- results$male$model
              male_value_transformed_flag <- results$male$transformed_flag
              combined_clustered_data <- dplyr::bind_rows(combined_clustered_data, results$male$clustered_data)
            }
  
            female_data <- gmm_data %>% dplyr::filter(Gender == "Female")
            results$female <- run_gmm_analysis_on_subset(
              data_subset = female_data,
              gender_label = "Female",
              value_col_name = value_col,
              age_col_name = age_col,
              message_rv = message_rv,
              progress_increment = 0.3,
              model_names = model_names_to_use,
              g_range = g_range_to_use
            )
            if (!is.null(results$female$model)) {
              female_gmm_model_bic <- results$female$model
              female_value_transformed_flag <- results$female$transformed_flag
              combined_clustered_data <- dplyr::bind_rows(combined_clustered_data, results$female$clustered_data)
            }
          }
        } else if (gender_choice == "Male") {
          male_data <- gmm_data %>% dplyr::filter(Gender == "Male")
          results$male <- run_gmm_analysis_on_subset(
            data_subset = male_data,
            gender_label = "Male",
            value_col_name = value_col,
            age_col_name = age_col,
            message_rv = message_rv,
            progress_increment = 0.6,
            model_names = model_names_to_use,
            g_range = g_range_to_use
          )
          if (!is.null(results$male$model)) {
            male_gmm_model_bic <- results$male$model
            male_value_transformed_flag <- results$male$transformed_flag
            combined_clustered_data <- dplyr::bind_rows(combined_clustered_data, results$male$clustered_data)
          }
        } else if (gender_choice == "Female") {
          female_data <- gmm_data %>% dplyr::filter(Gender == "Female")
          results$female <- run_gmm_analysis_on_subset(
            data_subset = female_data,
            gender_label = "Female",
            value_col_name = value_col,
            age_col_name = age_col,
            message_rv = message_rv,
            progress_increment = 0.6,
            model_names = model_names_to_use,
            g_range = g_range_to_use
          )
          if (!is.null(results$female$model)) {
            female_gmm_model_bic <- results$female$model
            female_value_transformed_flag <- results$female$transformed_flag
            combined_clustered_data <- dplyr::bind_rows(combined_clustered_data, results$female$clustered_data)
          }
        }
        
        gmm_models_bic_rv(list(
          combined = combined_gmm_model_bic,
          male = male_gmm_model_bic,
          female = female_gmm_model_bic
        ))
        gmm_transformation_details_rv(list(
          male_value_transformed = male_value_transformed_flag, 
          female_value_transformed = female_value_transformed_flag
        ))
  
        if (nrow(combined_clustered_data) > 0) {
          gmm_processed_data_rv(list(bic = combined_clustered_data, removed_rows = removed_rows_count))
          message_rv(list(text = "GMM analysis complete!", type = "success"))
        } else {
          message_rv(list(text = "No data available after GMM processing for plotting/summary.", type = "error"))
          gmm_processed_data_rv(NULL)
        }
  
        incProgress(0.1, detail = "Generating plots and summaries...")
      })
    }, error = function(e) {
      message_rv(list(text = paste("Analysis Error:", e$message), type = "danger"))
      gmm_processed_data_rv(NULL)
    }, finally = {
      analysis_running_rv(FALSE)
      shinyjs::enable("run_gmm_analysis_btn")
      shinyjs::runjs("$('#run_gmm_analysis_btn').text('Analyze');")
      session$sendCustomMessage('analysisStatus', FALSE)
    })
  })

  # Observer for the Reset button on the GMM tab
  observeEvent(input$reset_gmm_analysis_btn, {
    # Reset all reactive values, allowing the UI to update automatically
    gmm_uploaded_data_rv(NULL)
    gmm_processed_data_rv(NULL)
    gmm_transformation_details_rv(list(male_value_transformed = FALSE, female_value_transformed = FALSE))
    gmm_models_bic_rv(list(male = NULL, female = NULL))
    shinyjs::reset("gmm_file_upload")
    message_rv(list(text = "GMM data and results reset.", type = "info"))
    
    updateSelectInput(session, "gmm_value_col", choices = c("None" = ""), selected = "")
    updateSelectInput(session, "gmm_age_col", choices = c("None" = ""), selected = "")
    updateSelectInput(session, "gmm_gender_col", choices = c("None" = ""), selected = "")
  })

  # Renders the UI for GMM results
  output$gmm_results_ui <- renderUI({
    results <- gmm_processed_data_rv()

    if (is.null(results) || nrow(results$bic) == 0) {
      return(NULL)
    }

    tagList(
      div(class = "output-box",
          card(
            h4(class = "gmm-title", "BIC Criterion Results"),
            plotOutput("gmm_bic_plots", height = "400px"),
            card_footer(
              "Plot Description:",
              tooltip(
                bs_icon("info-circle"),
                "The primary purpose of this plot is to visually identify the optimal GMM model for the given dataset. A GMM analysis aims to cluster patient data into different subpopulations based on their Value and Age. The model with the highest BIC value is chosen as the best fit, as it balances model complexity (number of components) with how well the model explains the data."
              )
            )
          ),
          div(class = "spacing-div"),
          plotOutput("plot_output_gmm_bic", height = "600px"),
          div(class = "spacing-div"),
          verbatimTextOutput("gmm_summary_output_bic")
      )
    )
  })
  
  # Renders the BIC plots
  output$gmm_bic_plots <- renderPlot({
    models <- gmm_models_bic_rv()

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mar = c(5.1, 4.1, 1.1, 2.1), mgp = c(2.5, 1, 0))
    
    has_combined_model <- !is.null(models$combined)
    has_male_model <- !is.null(models$male)
    has_female_model <- !is.null(models$female)

    if (has_combined_model) {
      if (!is.null(models$combined) && !inherits(models$combined, "try-error")) {
        plot_title <- paste0("BIC Plot for ", input$gmm_value_col, " and ", input$gmm_age_col, " (Combined)")
        plot(models$combined, what = "BIC", main = plot_title)
      } else {
        return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "GMM model for combined data was not generated.", size = 6, color = "grey50"))
      }
    } else if (has_male_model || has_female_model) {
      par(mfrow = c(1, 2))
      if (has_male_model && !inherits(models$male, "try-error")) {
        plot_title <- paste0("BIC Plot for ", input$gmm_value_col, " and ", input$gmm_age_col, " (Male)")
        plot(models$male, what = "BIC", main = plot_title)
      } else {
        plot.new()
        text(0.5, 0.5, "GMM model for male data was not generated.")
      }
      if (has_female_model && !inherits(models$female, "try-error")) {
        plot_title <- paste0("BIC Plot for ", input$gmm_value_col, " and ", input$gmm_age_col, " (Female)")
        plot(models$female, what = "BIC", main = plot_title)
      } else {
        plot.new()
        text(0.5, 0.5, "GMM model for female data was not generated.")
      }
      par(mfrow = c(1, 1))
    } else {
      return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No GMM models available for plotting.", size = 6, color = "grey50"))
    }
  })

  # Renders the clustered scatter plot
  output$plot_output_gmm_bic <- renderPlot({
    plot_data <- gmm_processed_data_rv()$bic
    if (is.null(plot_data) || nrow(plot_data) == 0) {
      return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No GMM data available for plotting.", size = 6, color = "grey50"))
    }
    plot_value_age(plot_data,
                value_col_name = input$gmm_value_col,
                age_col_name = input$gmm_age_col)
  })  
  
  # A helper function to generate the summary for a given subpopulation
  generate_subpop_summary <- function(model, plot_data, gender_label, value_col, age_col) {
    cat(paste0("\n--- ", gender_label, " Subpopulations ---\n"))
    print(summary(model))
    num_clusters <- model$G
    cat("\n")
    for (i in 1:num_clusters) {
      cat(paste0("Cluster ", i, ":\n"))
      cat(paste0("  Proportion: ", round(model$parameters$pro[i], 3), "\n"))
      
      cluster_data <- plot_data %>% dplyr::filter(Gender == gender_label, cluster == i)
      mean_value <- mean(cluster_data$Value, na.rm = TRUE)
      mean_age <- mean(cluster_data$Age, na.rm = TRUE)
      sd_value <- sd(cluster_data$Value, na.rm = TRUE)
      sd_age <- sd(cluster_data$Age, na.rm = TRUE)
      
      cat(paste0("  Mean ", value_col, ": ", round(mean_value, 3), "\n"))
      cat(paste0("  Mean ", age_col, ": ", round(mean_age, 3), "\n"))
      cat(paste0("  Std Dev ", value_col, ": ", round(sd_value, 3), "\n"))
      cat(paste0("  Std Dev ", age_col, ": ", round(sd_age, 3), "\n"))
      
      if (!is.na(sd_age)) {
        lower_age <- round(mean_age - 2 * sd_age, 1)
        upper_age <- round(mean_age + 2 * sd_age, 1)
        cat(paste0("  Estimated ", age_col, " Range (Mean +/- 2SD): [", max(0, lower_age), " to ", upper_age, "] years\n"))
      } else {
        cat(paste0("  Estimated ", age_col, " Range: N/A (Std Dev Age problematic)\n"))
      }
      cat("\n")
    }
  }

  # Renders the GMM summary text
  output$gmm_summary_output_bic <- renderPrint({
    results <- gmm_processed_data_rv()

    if (is.null(results) || nrow(results$bic) == 0) {
      return("No GMM analysis results to display.")
    }

    cat("--- GMM Analysis Summary (BIC Criterion) ---\n")
    
    if (!is.null(results$removed_rows)) {
      cat(paste0("Note: ", results$removed_rows, " rows were removed due to missing data.\n"))
    }
    
    plot_data <- results$bic
    models <- gmm_models_bic_rv()
    
    if (!is.null(models$combined) && !inherits(models$combined, "try-error")) {
        generate_subpop_summary(models$combined, plot_data, "Combined", input$gmm_value_col, input$gmm_age_col)
    } else {
        if (!is.null(models$male) && !inherits(models$male, "try-error")) {
            generate_subpop_summary(models$male, plot_data, "Male", input$gmm_value_col, input$gmm_age_col)
        } else {
            cat("No male subpopulations detected.\n")
        }
        
        if (!is.null(models$female) && !inherits(models$female, "try-error")) {
            generate_subpop_summary(models$female, plot_data, "Female", input$gmm_value_col, input$gmm_age_col)
        } else {
            cat("No female subpopulations detected.\n")
        }
    }

    if (gmm_transformation_details_rv()$male_value_transformed || gmm_transformation_details_rv()$female_value_transformed) {
      cat("\nNote: ", input$gmm_value_col, " values were transformed (Yeo-Johnson) for GMM input due to skewness. Reported ", input$gmm_value_col, " values are original.\n")
    }
  })
}