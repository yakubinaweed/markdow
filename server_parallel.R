# server_parallel.R
# This module contains the logic for the "Parallel Analysis" tab.

# Load all necessary libraries.
library(shiny)
library(readxl)
library(tidyverse)
library(refineR)
library(shinyjs)
library(shinyWidgets)
library(bslib)
library(ggplot2)
library(future)
library(future.apply)
library(moments)

# =========================================================================
# UTILITY FUNCTIONS FOR PARALLEL ANALYSIS
# =========================================================================

# Helper function to guess column names based on common keywords
guess_column <- function(cols_available, common_names) {
  for (name in common_names) {
    match_idx <- grep(paste0("^", name, "$"), cols_available, ignore.case = TRUE)
    if (length(match_idx) > 0) {
      return(cols_available[match_idx[1]])
    }
  }
  return("")
}

# Filters data based on gender, age, and column names
filter_data <- function(data, gender_choice, age_min, age_max, col_gender, col_age) {
  if (col_age == "") {
    stop("Age column not found in data.")
  }

  filtered_data <- data %>%
    filter(!!rlang::sym(col_age) >= age_min & !!rlang::sym(col_age) <= age_max)

  if (col_gender != "" && col_gender %in% names(data)) {
    filtered_data <- filtered_data %>%
      mutate(Gender_Standardized = case_when(
        grepl("male|m|man|jongen(s)?|heren|mannelijk(e)?", !!rlang::sym(col_gender), ignore.case = TRUE) ~ "Male",
        grepl("female|f|vrouw(en)?|v|meisje(s)?|dame|mevr|vrouwelijke", !!rlang::sym(col_gender), ignore.case = TRUE) ~ "Female",
        TRUE ~ "Other"
      ))

    if (gender_choice != "Both") {
      filtered_data <- filtered_data %>%
        filter(Gender_Standardized == case_when(
          gender_choice == "M" ~ "Male",
          gender_choice == "F" ~ "Female",
          TRUE ~ NA_character_
        )) %>%
        filter(!is.na(Gender_Standardized))
    }
  } else {
    if (gender_choice %in% c("M", "F")) {
      return(data[FALSE, ])
    } else {
      filtered_data <- filtered_data %>%
        mutate(Gender_Standardized = "Combined")
    }
  }

  return(filtered_data)
}

# Parses the age ranges from text input
parse_age_ranges <- function(ranges_text, gender) {
  ranges <- strsplit(ranges_text, ",")[[1]]
  parsed_ranges <- list()
  for (range in ranges) {
    parts <- trimws(strsplit(range, "-")[[1]])
    if (length(parts) == 2) {
      age_min <- as.numeric(parts[1])
      age_max <- as.numeric(parts[2])
      if (!is.na(age_min) && !is.na(age_max) && age_min <= age_max) {
        parsed_ranges <- c(parsed_ranges, list(list(gender = gender, age_min = age_min, age_max = age_max)))
      } else {
        warning(paste("Invalid numeric range for", gender, ":", range))
      }
    } else if (trimws(range) != "") {
      warning(paste("Invalid format for", gender, ":", range))
    }
  }
  return(parsed_ranges)
}

# Wrapper for a single refineR analysis with data filtering and cleaning
run_single_refiner_analysis <- function(subpopulation, data, col_value, col_age, col_gender, model_choice, nbootstrap_value) {
  gender <- subpopulation$gender
  age_min <- subpopulation$age_min
  age_max <- subpopulation$age_max
  label <- paste0(gender, " (", age_min, "-", age_max, ")")

  tryCatch({
    filter_gender_choice <- ifelse(gender == "Male", "M", ifelse(gender == "Female", "F", "Both"))

    filtered_data_for_refiner <- filter_data(data,
                                 gender_choice = filter_gender_choice,
                                 age_min = age_min,
                                 age_max = age_max,
                                 col_gender = col_gender,
                                 col_age = col_age)

    value_col_name <- col_value
    
    if (!value_col_name %in% names(filtered_data_for_refiner)) {
      stop("Selected value column not found after filtering.")
    }

    original_rows_count <- nrow(filtered_data_for_refiner)

    cleaned_data <- filtered_data_for_refiner %>%
      mutate(!!rlang::sym(value_col_name) := as.numeric(!!rlang::sym(value_col_name))) %>%
      filter(!is.na(!!rlang::sym(value_col_name)))

    # Corrected: Set Gender_Standardized to the explicit gender of the subpopulation
    raw_subpopulation_data <- cleaned_data %>%
                                    rename(Age = !!rlang::sym(col_age), Value = !!rlang::sym(col_value)) %>%
                                    mutate(label = label, Gender_Standardized = gender)
    
    if (nrow(cleaned_data) == 0) {
      stop(paste("No data found for subpopulation:", label, "after cleaning."))
    }
    
    removed_rows_count <- original_rows_count - nrow(cleaned_data)

    final_model_choice <- if (model_choice == "AutoSelect") {
      skew <- moments::skewness(cleaned_data[[value_col_name]], na.rm = TRUE)
      if (abs(skew) > 0.5) {
        "modBoxCox"
      } else {
        "BoxCox"
      }
    } else {
      model_choice
    }

    model <- refineR::findRI(Data = cleaned_data[[value_col_name]],
                             NBootstrap = nbootstrap_value,
                             model = final_model_choice)

    if (is.null(model) || inherits(model, "try-error")) {
      stop(paste("RefineR model could not be generated for subpopulation:", label))
    }
    
    ri_data_fulldata <- getRI(model, RIperc = c(0.025, 0.975), pointEst = "fullDataEst")
    ri_data_median <- getRI(model, RIperc = c(0.025, 0.975), pointEst = "medianBS")

    ri_low_fulldata <- ri_data_fulldata$PointEst[ri_data_fulldata$Percentile == 0.025]
    ri_high_fulldata <- ri_data_fulldata$PointEst[ri_data_fulldata$Percentile == 0.975]
    
    ci_low_low <- ri_data_median$CILow[ri_data_median$Percentile == 0.025]
    ci_low_high <- ri_data_median$CIHigh[ri_data_median$Percentile == 0.025]
    ci_high_low <- ri_data_median$CILow[ri_data_median$Percentile == 0.975]
    ci_high_high <- ri_data_median$CIHigh[ri_data_median$Percentile == 0.975]


    list(
      label = label,
      model = model,
      raw_data = raw_subpopulation_data,
      removed_rows = removed_rows_count,
      age_min = age_min,
      age_max = age_max,
      ri_low_fulldata = ri_low_fulldata,
      ri_high_fulldata = ri_high_fulldata,
      ci_low_low = ci_low_low,
      ci_low_high = ci_low_high,
      ci_high_low = ci_high_low,
      ci_high_high = ci_high_high,
      status = "success",
      message = "Analysis complete.",
      final_model = final_model_choice
    )

  }, error = function(e) {
    list(
      label = label,
      status = "error",
      message = paste("Error:", e$message)
    )
  })
}

# Main server logic for the parallel tab
parallelServer <- function(input, output, session, parallel_data_rv, parallel_results_rv, parallel_message_rv, analysis_running_rv) {
  
  # Reactive value to store all raw data from successful analyses
  combined_raw_data_rv <- reactiveVal(tibble())
  
  # Observer for file upload, updating column selectors
  observeEvent(input$parallel_file, {
    req(input$parallel_file)
    tryCatch({
      data <- readxl::read_excel(input$parallel_file$datapath)
      parallel_data_rv(data)
      parallel_message_rv(list(type = "success", text = "Data file uploaded successfully."))

      col_names <- colnames(data)
      all_col_choices_with_none <- c("None" = "", col_names)

      updateSelectInput(session, "parallel_col_value", choices = all_col_choices_with_none, selected = guess_column(col_names, c("HB_value", "Value", "Result", "Measurement", "Waarde")))
      updateSelectInput(session, "parallel_col_age", choices = all_col_choices_with_none, selected = guess_column(col_names, c("leeftijd", "age", "AgeInYears", "Years")))
      updateSelectInput(session, "parallel_col_gender", choices = all_col_choices_with_none, selected = guess_column(col_names, c("geslacht", "gender", "sex", "Gender", "Sex")))
    }, error = function(e) {
      parallel_message_rv(list(type = "error", text = paste("Error loading file:", e$message)))
      parallel_data_rv(NULL)
    })
  })

  # Observer for the Run Parallel Analysis button
  observeEvent(input$run_parallel_btn, {
    if (analysis_running_rv()) {
      parallel_message_rv(list(text = "An analysis is already running. Please wait or reset.", type = "warning"))
      return()
    }

    # Stage 1: Pre-run checks and setup
    req(parallel_data_rv(), input$parallel_col_value, input$parallel_col_age)
    if (input$parallel_col_value == "" || input$parallel_col_age == "") {
      parallel_message_rv(list(text = "Please select the value and age columns.", type = "error"))
      return()
    }
    
    subpopulations <- c(
      parse_age_ranges(input$male_age_ranges, "Male"),
      parse_age_ranges(input$female_age_ranges, "Female"),
      parse_age_ranges(input$combined_age_ranges, "Combined")
    )

    if (length(subpopulations) == 0) {
      parallel_message_rv(list(text = "Please enter valid age ranges in the format 'min-max' for at least one gender.", type = "error"))
      return()
    }

    parallel_message_rv(list(text = "Starting parallel analysis...", type = "info"))
    analysis_running_rv(TRUE)
    shinyjs::disable("run_parallel_btn")
    shinyjs::runjs("$('#run_parallel_btn').text('Analyzing...');")
    session$sendCustomMessage('analysisStatus', TRUE)

    # Stage 2: Attempt to set up the future plan and run the parallel analysis
    results_list <- NULL
    
    # New tryCatch block to handle errors during the plan() setup
    tryCatch({
      plan(multisession, workers = input$cores)

      # Run the parallel analysis using `future_lapply`
      data_to_analyze <- parallel_data_rv()
      col_value_input <- input$parallel_col_value
      col_age_input <- input$parallel_col_age
      col_gender_input <- input$parallel_col_gender
      model_choice_input <- input$parallel_model_choice
      nbootstrap_value_input <- input$parallel_nbootstrap_speed

      results_list <- future_lapply(subpopulations, function(sub) {
        run_single_refiner_analysis(
          subpopulation = sub,
          data = data_to_analyze,
          col_value = col_value_input,
          col_age = col_age_input,
          col_gender = col_gender_input,
          model_choice = model_choice_input,
          nbootstrap_value = nbootstrap_value_input
        )
      }, future.seed = TRUE)
      
      # If plan and lapply succeed, process results and set final success message
      if (!is.null(results_list)) {
        parallel_results_rv(results_list)

        raw_data_list <- lapply(results_list, function(r) {
          if (r$status == "success") {
            return(r$raw_data)
          } else {
            return(NULL)
          }
        })
        
        combined_raw_data_rv(bind_rows(raw_data_list))

        if (all(sapply(parallel_results_rv(), function(r) r$status == "error"))) {
          # This case handles when the plan worked, but all individual analyses failed.
          parallel_message_rv(list(text = "Parallel analysis failed for all subpopulations. See the summary for complete error.", type = "error"))
        } else {
          parallel_message_rv(list(text = "Parallel analysis complete!", type = "success"))
        }
      }

    }, error = function(e) {
      # This block now handles all errors, including those from plan().
      error_message <- e$message
      
      # Check if the error message is about too many cores
      if (grepl("too many cores", error_message, ignore.case = TRUE) ||
          grepl("load is set to", error_message, ignore.case = TRUE) ||
          grepl("connections", error_message, ignore.case = TRUE)) {
        
        parallel_message_rv(list(text = paste("Failed to start parallel processes: The number of cores requested (", input$cores, ") exceeds the available resources or the hard limit. Please reduce the number of cores.", sep = ""), type = "error"))
      } else {
        # Generic error message for other failures
        parallel_message_rv(list(text = paste("An unexpected error occurred during parallel analysis:", error_message), type = "error"))
      }

      # Clear out any potential partial results to avoid downstream errors
      parallel_results_rv(list())
      combined_raw_data_rv(tibble())
      
    }, finally = {
      # This block runs regardless of success or failure to clean up the UI
      analysis_running_rv(FALSE)
      shinyjs::enable("run_parallel_btn")
      shinyjs::runjs("$('#run_parallel_btn').text('Run Parallel Analysis');")
      session$sendCustomMessage('analysisStatus', FALSE)
    })
  })

  # Observer for the Reset button
  observeEvent(input$reset_parallel_btn, {
    parallel_data_rv(NULL)
    parallel_results_rv(list())
    combined_raw_data_rv(tibble())
    parallel_message_rv(list(type = "", text = ""))
    shinyjs::reset("parallel_file")
    updateSelectInput(session, "parallel_col_value", choices = c("None" = ""), selected = "")
    updateSelectInput(session, "parallel_col_age", choices = c("None" = ""), selected = "")
    updateSelectInput(session, "parallel_col_gender", choices = c("None" = ""), selected = "")
    updateRadioButtons(session, "parallel_model_choice", selected = "BoxCox")
    updateSliderInput(session, "parallel_nbootstrap_speed", value = 50)
    updateTextAreaInput(session, "male_age_ranges", value = "")
    updateTextAreaInput(session, "female_age_ranges", value = "")
    updateTextAreaInput(session, "combined_age_ranges", value = "")
  })
  
  # Reactive expression to create a single combined table
  combined_summary_table <- reactive({
    results <- parallel_results_rv()
    if (is.null(results) || length(results) == 0) {
      return(NULL)
    }
    
    # Filter results based on the gender selection
    filtered_results <- Filter(function(r) {
        r$status == "success" && (str_extract(r$label, "^\\w+") %in% input$parallel_gender_filter)
    }, results)

    if (length(filtered_results) == 0) {
        return(NULL)
    }

    table_rows <- lapply(filtered_results, function(result) {
      tibble(
        label = result$label,
        Gender = str_extract(result$label, "^\\w+"),
        `Age Range` = paste0(result$age_min, "-", result$age_max),
        age_min = result$age_min,
        age_max = result$age_max,
        `CI Lower (Lower)` = round(result$ci_low_low, 3),
        `RI Lower` = round(result$ri_low_fulldata, 3),
        `CI Lower (Upper)` = round(result$ci_low_high, 3),
        `CI Upper (Lower)` = round(result$ci_high_low, 3),
        `RI Upper` = round(result$ri_high_fulldata, 3),
        `CI Upper (Upper)` = round(result$ci_high_high, 3)
      )
    })

    bind_rows(table_rows)
  })
  
  # Reactive expression to filter and prepare data for all combined plots
  filtered_plot_data_rv <- reactive({
    plot_data_from_raw <- combined_raw_data_rv()
    if (is.null(plot_data_from_raw) || nrow(plot_data_from_raw) == 0) {
      return(tibble())
    }
    
    # Filter data based on selected genders here
    req(input$parallel_gender_filter)
    plot_data_from_raw <- plot_data_from_raw %>%
      filter(Gender_Standardized %in% input$parallel_gender_filter)
    
    return(plot_data_from_raw)
  })

  # Dynamic UI for results
  output$parallel_results_ui <- renderUI({
    results <- parallel_results_rv()
    
    combined_table_data <- combined_summary_table()
    
    ui_elements <- tagList()
    
    if (is.null(combined_table_data)) {
      return(NULL)
    }

    ui_elements <- tagList(
      ui_elements,
      div(class = "summary-box", h5("Summary Table of Reference Intervals")),
      renderTable(combined_table_data %>% dplyr::select(-label, -age_min, -age_max), striped = TRUE, bordered = TRUE),
      br()
    )


    if (length(results) > 0) {
      individual_elements <- lapply(seq_along(results), function(i) {
        result <- results[[i]]
        if (result$status == "success") {
          label_parts <- unlist(strsplit(result$label, " "))
          gender_part <- label_parts[1]
          age_range_part <- gsub("[()]", "", label_parts[2])

          tagList(
            plotOutput(paste0("parallel_plot_", i)),
            verbatimTextOutput(paste0("parallel_summary_", i)),
            hr()
          )
        } else {
          div(class = "alert alert-danger", result$message)
        }
      })
      ui_elements <- tagList(ui_elements, do.call(tagList, individual_elements))
    }

    if (length(ui_elements) > 0) {
      ui_elements
    } else {
      NULL
    }
  })

  # Renders the combined dumbbell plot
  output$combined_dumbbell_plot <- renderPlot({
    plot_data_summary <- combined_summary_table()
    
    if (is.null(plot_data_summary) || nrow(plot_data_summary) == 0) {
      return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No successful reference intervals to plot for the selected genders.", size = 6, color = "grey50"))
    }
    
    # Sort the data by gender and age range to ensure a logical plot order
    plot_data_summary$label <- factor(plot_data_summary$label, 
                                     levels = unique(plot_data_summary$label[order(plot_data_summary$Gender, plot_data_summary$age_min)]))

    unit_label <- if (!is.null(input$parallel_unit_input) && input$parallel_unit_input != "") {
      paste0(input$parallel_col_value, " [", input$parallel_unit_input, "]")
    } else {
      input$parallel_col_value
    }

    gender_colors <- c("Male" = "steelblue", "Female" = "darkred", "Combined" = "darkgreen")

    ggplot2::ggplot(plot_data_summary, ggplot2::aes(y = label)) +
      ggplot2::geom_segment(ggplot2::aes(x = `CI Lower (Lower)`,
                                        xend = `CI Lower (Upper)`,
                                        y = label,
                                        yend = label,
                                        color = Gender),
                            linewidth = 10,
                            alpha = 0.3,
                            lineend = "square") +
      
      ggplot2::geom_segment(ggplot2::aes(x = `CI Upper (Lower)`,
                                        xend = `CI Upper (Upper)`,
                                        y = label,
                                        yend = label,
                                        color = Gender),
                            linewidth = 10,
                            alpha = 0.3,
                            lineend = "square") +

      ggplot2::geom_errorbarh(ggplot2::aes(xmin = `RI Lower`,
                                          xmax = `RI Upper`,
                                          color = Gender),
                              height = 0.1, linewidth = 1.2) +            

      ggplot2::geom_point(ggplot2::aes(x = `RI Lower`, color = Gender), shape = 18, size = 4) +
      ggplot2::geom_point(ggplot2::aes(x = `RI Upper`, color = Gender), shape = 18, size = 4) +
      
      ggplot2::facet_wrap(~ Gender, ncol = 1, scales = "free_y", strip.position = "right") +
      
      ggplot2::labs(
        title = "Estimated Reference Intervals by Subpopulation",
        x = unit_label,
        y = NULL,
        color = "Gender"
      ) +
      ggplot2::scale_color_manual(values = gender_colors, name = "Gender") +
      ggplot2::scale_fill_manual(values = gender_colors, guide = "none") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title.x = ggplot2::element_text(size = 14, margin = ggplot2::margin(t = 10)),
        axis.text = ggplot2::element_text(size = 12),
        strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
        strip.text = ggplot2::element_text(size = 12, face = "bold", color = "black"),
        legend.title = ggplot2::element_text(size = 12, face = "bold"),
        legend.text = ggplot2::element_text(size = 10),
        legend.position = "none"
      )
  })

  # Renders the combined RI plot
  output$combined_ri_plot <- renderPlot({
    plot_data_summary <- combined_summary_table()
    if (is.null(plot_data_summary) || nrow(plot_data_summary) == 0) {
       return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No successful reference intervals to plot for the selected genders.", size = 6, color = "grey50"))
    }
    
    # Sort the data by gender and age range to ensure a logical plot order
    plot_data_summary$label <- factor(plot_data_summary$label, 
                                     levels = unique(plot_data_summary$label[order(plot_data_summary$Gender, plot_data_summary$age_min)]))
    
    unit_label <- if (!is.null(input$parallel_unit_input) && input$parallel_unit_input != "") {
      paste0(input$parallel_col_value, " [", input$parallel_unit_input, "]")
    } else {
      input$parallel_col_value
    }

    gender_colors <- c("Male" = "steelblue", "Female" = "darkred", "Combined" = "darkgreen")

    ggplot2::ggplot(plot_data_summary) +
      ggplot2::geom_rect(ggplot2::aes(xmin = age_min, xmax = age_max, ymin = `CI Lower (Lower)`, ymax = `CI Lower (Upper)`, fill = Gender),
                         alpha = 0.2) +
      ggplot2::geom_rect(ggplot2::aes(xmin = age_min, xmax = age_max, ymin = `CI Upper (Lower)`, ymax = `CI Upper (Upper)`, fill = Gender),
                         alpha = 0.2) +
      ggplot2::geom_segment(ggplot2::aes(x = age_min, xend = age_max, y = `RI Lower`, yend = `RI Lower`, color = Gender),
                            linewidth = 1.2, linetype = "solid") +
      ggplot2::geom_segment(ggplot2::aes(x = age_min, xend = age_max, y = `RI Upper`, yend = `RI Upper`, color = Gender),
                            linewidth = 1.2, linetype = "solid") +
      ggplot2::geom_point(ggplot2::aes(x = age_min, y = `RI Lower`, color = Gender), size = 2) +
      ggplot2::geom_point(ggplot2::aes(x = age_max, y = `RI Lower`, color = Gender), size = 2) +
      ggplot2::geom_point(ggplot2::aes(x = age_min, y = `RI Upper`, color = Gender), size = 2) +
      ggplot2::geom_point(ggplot2::aes(x = age_max, y = `RI Upper`, color = Gender), size = 2) +
      ggplot2::labs(
        title = "Reference Intervals by Age and Gender",
        x = "Age",
        y = unit_label,
        color = "Gender",
        fill = "Gender (95% CI)"
      ) +
      ggplot2::scale_x_continuous(limits = c(0, 120)) +
      ggplot2::scale_color_manual(values = gender_colors) +
      ggplot2::scale_fill_manual(values = gender_colors, guide = "none") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = ggplot2::element_text(size = 14),
        axis.text = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 12, face = "bold"),
        legend.text = ggplot2::element_text(size = 10),
        legend.position = "bottom"
      )
  })
    
  # Renders the faceted density plot
  output$combined_density_plot <- renderPlot({
    plot_data <- filtered_plot_data_rv()
    results <- parallel_results_rv()
    
    if (is.null(plot_data) || nrow(plot_data) == 0) {
      return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data available for plotting.", size = 6, color = "grey50"))
    }
    
    req(input$parallel_gender_filter)
    plot_data <- plot_data %>%
      filter(Gender_Standardized %in% input$parallel_gender_filter)
    
    if (length(unique(plot_data$label)) < 1) {
       return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Insufficient data to create a faceted plot.", size = 6, color = "grey50"))
    }

    unit_label <- if (!is.null(input$parallel_unit_input) && input$parallel_unit_input != "") {
      paste0(input$parallel_col_value, " [", input$parallel_unit_input, "]")
    } else {
      input$parallel_col_value
    }
    
    ri_lines <- tibble()
    for (result in results) {
      if (result$status == "success") {
        ri_lines <- bind_rows(ri_lines, tibble(
          label = result$label,
          ri_low = result$ri_low_fulldata,
          ri_high = result$ri_high_fulldata
        ))
      }
    }
    
    ri_lines <- ri_lines %>%
        filter(str_extract(label, "^\\w+") %in% input$parallel_gender_filter)

    custom_colors <- c(
      "Male" = "steelblue", 
      "Female" = "darkred", 
      "Combined" = "darkgreen"
    )
    
    fill_colors <- setNames(custom_colors[str_extract(unique(plot_data$label), "^\\w+")], unique(plot_data$label))

    ggplot2::ggplot(plot_data, ggplot2::aes(x = Value, fill = label)) +
      ggplot2::geom_density(alpha = 0.6) +
      ggplot2::geom_vline(data = ri_lines, ggplot2::aes(xintercept = ri_low), linetype = "dashed", color = "darkred", size = 1) +
      ggplot2::geom_vline(data = ri_lines, ggplot2::aes(xintercept = ri_high), linetype = "dashed", color = "darkred", size = 1) +
      ggplot2::facet_wrap(~label, scales = "free_y") +
      ggplot2::labs(title = "Faceted Density Plot by Subpopulation",
                    x = unit_label,
                    y = "Density",
                    fill = "Subpopulation") +
      ggplot2::scale_fill_manual(values = fill_colors) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
        strip.text = ggplot2::element_text(size = 12),
        axis.title = ggplot2::element_text(size = 14),
        axis.text = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 12, face = "bold"),
        legend.text = ggplot2::element_text(size = 10),
        legend.position = "bottom"
      )
  })

  # Renders a single density plot for all selected subpopulations
  output$single_density_plot <- renderPlot({
    plot_data <- filtered_plot_data_rv()
    results <- parallel_results_rv()
    
    if (is.null(plot_data) || nrow(plot_data) == 0) {
      return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data available for plotting.", size = 6, color = "grey50"))
    }
    
    req(input$parallel_gender_filter)
    plot_data <- plot_data %>%
      filter(Gender_Standardized %in% input$parallel_gender_filter)
    
    unit_label <- if (!is.null(input$parallel_unit_input) && input$parallel_unit_input != "") {
      paste0(input$parallel_col_value, " [", input$parallel_unit_input, "]")
    } else {
      input$parallel_col_value
    }
    
    ri_lines <- tibble()
    for (result in results) {
      if (result$status == "success") {
        ri_lines <- bind_rows(ri_lines, tibble(
          label = result$label,
          ri_low = result$ri_low_fulldata,
          ri_high = result$ri_high_fulldata
        ))
      }
    }
    
    ri_lines <- ri_lines %>%
        filter(str_extract(label, "^\\w+") %in% input$parallel_gender_filter)
    
    custom_colors <- c(
      "Male" = "steelblue", 
      "Female" = "darkred", 
      "Combined" = "darkgreen"
    )
    
    fill_colors <- setNames(custom_colors[str_extract(unique(plot_data$label), "^\\w+")], unique(plot_data$label))

    ggplot2::ggplot(plot_data, ggplot2::aes(x = Value, fill = label)) +
      ggplot2::geom_density(alpha = 0.6) +
      ggplot2::geom_vline(data = ri_lines, ggplot2::aes(xintercept = ri_low, color = label), linetype = "dashed", size = 1, show.legend = FALSE) +
      ggplot2::geom_vline(data = ri_lines, ggplot2::aes(xintercept = ri_high, color = label), linetype = "dashed", size = 1) +
      ggplot2::labs(title = "Density Plot of Value Distribution",
                    x = unit_label,
                    y = "Density",
                    fill = "Subpopulation") +
      ggplot2::scale_fill_manual(values = fill_colors) +
      ggplot2::scale_color_manual(values = fill_colors, guide = "none") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = ggplot2::element_text(size = 14),
        axis.text = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 12, face = "bold"),
        legend.text = ggplot2::element_text(size = 10),
        legend.position = "bottom"
      )
  })

  # Renders the grouped box plot
  output$combined_box_plot <- renderPlot({
    plot_data <- filtered_plot_data_rv()
    
    if (is.null(plot_data) || nrow(plot_data) == 0) {
      return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data available for plotting.", size = 6, color = "grey50"))
    }
    
    req(input$parallel_gender_filter)
    plot_data <- plot_data %>%
      filter(Gender_Standardized %in% input$parallel_gender_filter)
    
    if (nrow(plot_data) == 0) {
      return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data available for plotting for the selected genders.", size = 6, color = "grey50"))
    }
    
    male_labels <- plot_data %>% dplyr::filter(grepl("Male", label)) %>% dplyr::arrange(label) %>% dplyr::pull(label) %>% unique()
    female_labels <- plot_data %>% dplyr::filter(grepl("Female", label)) %>% dplyr::arrange(label) %>% dplyr::pull(label) %>% unique()
    combined_labels <- plot_data %>% dplyr::filter(grepl("Combined", label)) %>% dplyr::arrange(label) %>% dplyr::pull(label) %>% unique()

    custom_order <- c(male_labels, female_labels, combined_labels)
    
    remaining_labels <- setdiff(unique(plot_data$label), custom_order)
    if(length(remaining_labels) > 0) {
        custom_order <- c(custom_order, remaining_labels)
    }
    
    plot_data$label <- factor(plot_data$label, levels = custom_order)

    unit_label <- if (!is.null(input$parallel_unit_input) && input$parallel_unit_input != "") {
      paste0(input$parallel_col_value, " [", input$parallel_unit_input, "]")
    } else {
      input$parallel_col_value
    }
    
    custom_colors <- c(
      "Male" = "steelblue", 
      "Female" = "darkred", 
      "Combined" = "darkgreen"
    )
    
    fill_colors <- setNames(custom_colors[str_extract(unique(plot_data$label), "^\\w+")], unique(plot_data$label))

    ggplot2::ggplot(plot_data, ggplot2::aes(x = label, y = Value, fill = label)) +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.colour = "red", outlier.shape = 8) +
      ggplot2::labs(
        title = "Summary of Value Distribution by Subpopulation",
        x = "Subpopulation",
        y = unit_label,
        fill = "Subpopulation"
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = fill_colors) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = ggplot2::element_text(size = 14),
        axis.text = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 12, face = "bold"),
        legend.text = ggplot2::element_text(size = 10)
      )
  })
  
  # Renders the combined text summary for all successful subpopulations
  output$combined_summary <- renderPrint({
    results <- parallel_results_rv()
    if (is.null(results) || length(results) == 0) {
      cat("No parallel analysis results to summarize yet.")
      return(NULL)
    }
    
    cat("--- Combined Summary of Reference Intervals ---\n\n")
    
    cat("Table Column Key:\n")
    cat("  RI Lower/Upper: The estimated Reference Interval limits.\n")
    cat("  CI Lower (Lower)/CI Lower (Upper): The Confidence Interval for the RI Lower limit.\n")
    cat("  CI Upper (Lower)/CI Upper (Upper): The Confidence Interval for the RI Upper limit.\n\n")

    has_successful_results <- FALSE
    for (r in results) {
      if (r$status == "success") {
        has_successful_results <- TRUE
        
        ri_low <- r$ri_low_fulldata
        ri_high <- r$ri_high_fulldata
        
        ci_low_low <- r$ci_low_low
        ci_low_high <- r$ci_low_high
        ci_high_low <- r$ci_high_low
        ci_high_high <- r$ci_high_high
        
        cat(paste0("Subpopulation: ", r$label, "\n"))
        cat(paste0("  Sample Size ", nrow(r$raw_data), "\n"))
        cat(paste0("  Rows Removed: ", r$removed_rows, "\n"))
        cat(paste0("  Estimated RI Lower Limit: ", round(ri_low, 3), "\n"))
        cat(paste0("  Confidence Interval for Lower Limit: [", round(ci_low_low, 3), ", ", round(ci_low_high, 3), "]\n"))
        cat(paste0("  Estimated RI Upper Limit: ", round(ri_high, 3), "\n"))
        cat(paste0("  Confidence Interval for Upper Limit: [", round(ci_high_low, 3), ", ", round(ci_high_high, 3), "]\n"))
        
        if(r$final_model != input$parallel_model_choice) {
          cat(paste0("  Transformation Model: Auto-selected ", r$final_model, " (from user's '", input$parallel_model_choice, "' choice)\n"))
        } else {
          cat(paste0("  Transformation Model: ", r$final_model, "\n"))
        }

        if (!is.null(input$parallel_unit_input) && input$parallel_unit_input != "") {
          cat(paste0("  Unit of Measurement: ", input$parallel_unit_input, "\n"))
        }
        cat("\n")
      } else if (r$status == "error") {
        cat(paste0("Subpopulation: ", r$label, " (Analysis Failed)\n"))
        cat(paste0("  Reason: ", r$message, "\n\n"))
      }
    }

    if (!has_successful_results) {
      cat("No successful reference intervals were found to summarize.\n")
    }
  })
# =========================================================================
# DYNAMICALLY RENDER INDIVIDUAL SUBPOPULATION PLOTS AND SUMMARIES
# =========================================================================
  # A function to render the plot and summary for a single subpopulation result
  render_subpop_results <- function(i, result) {
    output_id_plot <- paste0("parallel_plot_", i)
    output_id_summary <- paste0("parallel_summary_", i)
    model <- result$model
    
    label_parts <- unlist(strsplit(result$label, " "))
    gender_part <- label_parts[1]
    age_range_part <- gsub("[()]", "", label_parts[2])
    
    output[[output_id_plot]] <- renderPlot({
      req(model)
      
      value_col_name <- input$parallel_col_value
      model_type <- switch(result$final_model,
                           "BoxCox" = " (BoxCox Transformed)",
                           "modBoxCox" = " (modBoxCox Transformed)")
      
      plot_title <- paste0("Estimated Reference Intervals for ", value_col_name, 
                           model_type, " (Gender: ", gender_part, ", Age: ", age_range_part, ")")
      
      xlab_text <- if (!is.null(input$parallel_unit_input) && input$parallel_unit_input != "") {
        paste0(value_col_name, " ", "[", input$parallel_unit_input, "]")
      } else {
        value_col_name
      }
      
      plot(model, showCI = TRUE, RIperc = c(0.025, 0.975), showPathol = FALSE,
           title = plot_title,
           xlab = xlab_text)
    })

    output[[output_id_summary]] <- renderPrint({
        req(model)
        cat("--- RefineR Summary for ", input$parallel_col_value, " (Gender: ", gender_part, ", Age: ", age_range_part, ") ---\n")
        cat(paste0("Note: ", result$removed_rows, " rows were removed due to missing or invalid data.\n"))
        print(model)
    })
  }

  # Main observer to trigger the rendering of all individual results
  observe({
    results <- parallel_results_rv()
    if (length(results) > 0) {
      for (i in seq_along(results)) {
        result <- results[[i]]
        if (result$status == "success") {
          render_subpop_results(i, result)
        }
      }
    }
  })
}