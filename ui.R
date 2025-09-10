# ui.R
library(shiny)
library(bslib)
library(refineR)
library(readxl)
library(moments)
library(shinyjs)
library(shinyWidgets)
library(shinyFiles)
library(ggplot2)
library(bsicons)

# Defines the main UI using a navigation bar with three tabs
ui <- navbarPage(
  title = "RefineR Reference Interval Estimation",
  id = "tabs",
  theme = bs_theme(version = 5, base_font = font_google("Inter"), heading_font = font_google("Rethink Sans"), font_scale = 1.1, bootswatch = "default"),

  # First tab for the main RefineR analysis
  tabPanel(
    title = "Main Analysis",
    useShinyjs(),
    # JavaScript to manage tab disabling during analysis
    tags$head(
      includeCSS("www/styles.css"),
      tags$script(HTML("
        var analysisRunning = false;
        Shiny.addCustomMessageHandler('analysisStatus', function(status) {
          analysisRunning = status;
          if (status) {
            $('a[data-toggle=\"tab\"]').each(function() {
              if (!$(this).parent().hasClass('active')) {
                $(this).addClass('disabled-tab-link');
              }
            });
          } else {
            $('a.disabled-tab-link').removeClass('disabled-tab-link');
          }
        });
        $(document).on('click', 'a.disabled-tab-link', function(event) {
          event.preventDefault();
          Shiny.setInputValue('tab_switch_blocked', new Date().getTime());
          return false;
        });
      "))
    ),
    sidebarLayout(
      sidebarPanel(
        style = "padding-right: 15px;",
        
        div(class = "card", style = "border: 1px solid #ccc; border-radius: 8px;",
          div(class = "card-header", style = "background-color: #f7f7f7; padding: 10px; border-bottom: 1px solid #ccc; border-top-left-radius: 8px; border-top-right-radius: 8px;",
            h5(
              tooltip(
                trigger = list("Main Analysis Inputs", bs_icon("info-circle")),
                "This section contains the core inputs for filtering data and running the main RefineR analysis, including gender, age range, and model selection."
              ),
              style = "margin-top: 0; margin-bottom: 0;"
            )
          ),
          div(class = "card-body", style = "padding: 15px;",
            selectInput(inputId = "gender_choice", label = "Select Gender:", choices = c("Male" = "M", "Female" = "F", "Both" = "Both"), selected = "Both"),
            sliderInput(inputId = "age_range", label = "Age Range:", min = 0, max = 100, value = c(0, 100), step = 1)
          )
        ),
        
        br(),
        fileInput(inputId = "data_file", label = "Upload Data (Excel File)", accept = c(".xlsx")),
        hr(),
        selectInput(inputId = "col_value", label = "Select Column for Values:", choices = c("None" = ""), selected = ""),
        selectInput(inputId = "col_age", label = "Select Column for Age:", choices = c("None" = ""), selected = ""),
        selectInput(inputId = "col_gender", label = "Select Column for Gender:", choices = c("None" = ""), selected = ""),
        hr(),
        sliderInput(
          inputId = "nbootstrap_speed",
          label = tags$span(
            tooltip(
              trigger = list(tags$span(bs_icon("info-circle"))),
              "Higher values mean more bootstrap iterations for increased accuracy, but will result in slower analysis times (1 = Fast, 50 = Medium, 200 = Slow)."
            ),
            "Select Computation Speed:"
          ),
          min = 1, max = 200, value = 50, step = 1
        ),
        radioButtons(inputId = "model_choice",
                     label = tags$span(
                       tooltip(
                         trigger = list(tags$span(bs_icon("info-circle"))),
                         "BoxCox: For positive-valued data with light to moderate skewness. modBoxCox: For data with high skewness or values close to zero. Auto-select: Automatically chooses the optimal transformation based on data skewness."
                       ),
                       "Select Transformation Model:"
                     ),
                     choices = c("BoxCox" = "BoxCox", "modBoxCox" = "modBoxCox", "Auto-select" = "AutoSelect"),
                     selected = "AutoSelect", inline = TRUE),
        
        actionButton("analyze_btn", "Run Analysis", class = "btn-primary"),
        actionButton("reset_btn", "Reset File", class = "btn-secondary"),
        shinyFiles::shinyDirButton(id = "select_dir_btn", label = "Select Output Directory", title = "Select a directory to save plots", style = "margin-top: 5px;"),
        div(style = "margin-top: 5px; display: flex; align-items: center; justify-content: flex-start; width: 100%;",
            prettySwitch(inputId = "enable_directory", label = "Auto-Save Graph", status = "success", fill = TRUE, inline = TRUE)
        ),
        uiOutput("main_message"),
        hr(),
        numericInput("ref_low", "Reference Lower Limit:", value = NA),
        numericInput("ref_high", "Reference Upper Limit:", value = NA),
        textInput(inputId = "unit_input", label = "Unit of Measurement", value = "mmol/L", placeholder = "ex. g/L")
      ),
      mainPanel(
        plotOutput("result_plot"),
        verbatimTextOutput("result_text")
      )
    )
  ),

  # Second tab for Gaussian Mixture Model (GMM) subpopulation detection
  tabPanel(
    title = "Subpopulation Detection (GMM)",
    useShinyjs(),
    sidebarLayout(
      sidebarPanel(
        div(class = "card", style = "border: 1px solid #ccc; border-radius: 8px;",
          div(class = "card-header", style = "background-color: #f7f7f7; padding: 10px; border-bottom: 1px solid #ccc; border-top-left-radius: 8px; border-top-right-radius: 8px;",
            h5(
              tooltip(
                trigger = list("GMM Analysis", bs_icon("info-circle")),
                "Gaussian Mixture Models (GMM) detect hidden subpopulations. The mclust package selects the best model and components using BIC. Data is preprocessed with Yeo-Johnson transformation (if skewed) and standardization for values and age."
              ),
              style = "margin-top: 0; margin-bottom: 0;"
            )
          ),
          div(class = "card-body", style = "padding: 15px;",
            fileInput(inputId = "gmm_file_upload", label = "Upload Data (Excel File)", accept = c(".xlsx")),
            hr(),
            selectInput(inputId = "gmm_value_col", label = "Select Column for Values:", choices = c("None" = ""), selected = ""),
            selectInput(inputId = "gmm_age_col", label = "Select Column for Age:", choices = c("None" = ""), selected = ""),
            selectInput(
              inputId = "gmm_gender_col",
              label = tags$span(
                tooltip(
                  trigger = list(tags$span(bs_icon("info-circle")), "Select Column for Gender:"),
                  "Optional. If not selected, analysis will be run on combined data."
                )
              ),
              choices = c("None" = ""),
              selected = ""
            ),
            hr(),
            uiOutput("gmm_gender_choice_ui"),
            radioButtons(
              inputId = "gmm_model_selection_choice",
              label = "Select BIC Model Option:",
              choices = c("Auto-select", "Manual Selection"),
              selected = "Auto-select"
            ),
            uiOutput("gmm_manual_model_ui"),
            actionButton("run_gmm_analysis_btn", "Run Analysis", class = "btn-primary"),
            actionButton("reset_gmm_analysis_btn", "Reset File", class = "btn-secondary"),
            div(style = "margin-top: 15px;", uiOutput("app_message"))
          )
        )
      ),
      mainPanel(
        uiOutput("gmm_results_ui")
      )
    )
  ),

  # Third tab for Parallel RefineR Analysis
  tabPanel(
    title = "Parallel Analysis",
    useShinyjs(),
    sidebarLayout(
      sidebarPanel(
        style = "padding-right: 15px;",
        div(class = "card", style = "border: 1px solid #ccc; border-radius: 8px;",
          div(class = "card-header", style = "background-color: #f7f7f7; padding: 10px; border-bottom: 1px solid #ccc; border-top-left-radius: 8px; border-top-right-radius: 8px;",
            h5(
              tooltip(
                trigger = list("Parallel Analysis", bs_icon("info-circle")),
                "This tool runs multiple RefineR analyses for different subpopulations simultaneously using parallel processing, significantly speeding up computation time."
              ),
              style = "margin-top: 0; margin-bottom: 0;"
            )
          ),
          div(class = "card-body", style = "padding: 15px;",
            textAreaInput(
              inputId = "male_age_ranges",
              label = tags$span(tooltip(tags$span(bs_icon("info-circle")), "Enter age ranges for the male subpopulation. Use commas to separate multiple ranges."), "Male Age Ranges:"),
              rows = 1,
              placeholder = "e.g., 0-10, 10-20"
            ),
            textAreaInput(
              inputId = "female_age_ranges",
              label = tags$span(tooltip(tags$span(bs_icon("info-circle")), "Enter age ranges for the female subpopulation. Use commas to separate multiple ranges."), "Female Age Ranges:"),
              rows = 1,
              placeholder = "e.g., 0-10, 10-20"
            ),
            textAreaInput(
              inputId = "combined_age_ranges",
              label = tags$span(tooltip(tags$span(bs_icon("info-circle")), "Enter age ranges for the both subpopulations. Use commas to separate multiple ranges."), "All Genders Age Ranges:"),
              rows = 1,
              placeholder = "e.g., 0-10, 10-20"
            )
          )
        ),
        
        br(),
        fileInput(inputId = "parallel_file", label = "Upload Data (Excel File)", accept = c(".xlsx")),
        hr(),
        selectInput(inputId = "parallel_col_value", label = "Select Column for Values:", choices = c("None" = ""), selected = ""),
        selectInput(inputId = "parallel_col_age", label = "Select Column for Age:", choices = c("None" = ""), selected = ""),
        selectInput(inputId = "parallel_col_gender", label = "Select Column for Gender:", choices = c("None" = ""), selected = ""),
        hr(),
        sliderInput(
          inputId = "parallel_nbootstrap_speed",
          label = tags$span(tooltip(trigger = list(tags$span(bs_icon("info-circle"))), "Higher values mean more bootstrap iterations for increased accuracy, but will result in slower analysis times (1 = Fast, 50 = Medium, 200 = Slow)."), "Select Computation Speed:"),
          min = 1, max = 200, value = 50, step = 1
        ),
        radioButtons(inputId = "parallel_model_choice",
                     label = tags$span(tooltip(trigger = list(tags$span(bs_icon("info-circle"))), "BoxCox: For positive-valued data with light to moderate skewness. modBoxCox: For data with high skewness or values close to zero. Auto-select: Automatically chooses the optimal transformation based on data skewness."), "Select Transformation Model:"),
                     choices = c("BoxCox" = "BoxCox", "modBoxCox" = "modBoxCox", "Auto-select" = "AutoSelect"),
                     selected = "AutoSelect", inline = TRUE),
        
        div(class = "parallel-buttons",
            actionButton("run_parallel_btn", "Run Parallel Analysis", class = "btn-primary"),
            actionButton("reset_parallel_btn", "Reset File", class = "btn-secondary")
        ),
        div(style = "margin-top: 15px;", uiOutput("parallel_message")),
        hr(),
        numericInput("cores", "Number of Cores:", value = 1, min = 1),
        textInput(inputId = "parallel_unit_input", label = "Unit of Measurement", value = "", placeholder = "ex. g/L")
      ),
      mainPanel(
        tabsetPanel(
          type = "pills", id = "my-nav",
          tabPanel("Individual Results",
                   div(style = "margin-top: 15px;"),
                   uiOutput("parallel_results_ui")
          ),
          tabPanel("Combined Summary",
                   div(style = "margin-top: 15px;"),
                   div(class = "gender-filter-container",
                       tags$span("Select Genders to Display:"),
                       checkboxGroupInput(inputId = "parallel_gender_filter", label = NULL, choices = c("Male", "Female", "Combined"), selected = c("Male", "Female", "Combined"), inline = TRUE)
                   ),
                   plotOutput("combined_dumbbell_plot"),
                   div(class = "spacing-div"),
                   plotOutput("combined_ri_plot"),
                   div(class = "spacing-div"),
                   plotOutput("combined_density_plot"),
                   div(class = "spacing-div"),
                   plotOutput("single_density_plot"),
                   div(class = "spacing-div"),
                   card(
                    plotOutput("combined_box_plot"),
                    card_footer(
                      "Plot Description:",
                      tooltip(bs_icon("info-circle"), "The square in each box plot represents the middle 50% of the data, also known as the interquartile range (IQR). The line inside the box is the median, which is the midpoint of the HGB data. The whiskers extending from the box show the normal range of the data that is not considered an outlier. The red dots are outliers, which are values significantly different from the rest of their subpopulation and fall outside of the whiskers.")
                    )
                   ),
                   div(class = "spacing-div"),
                   verbatimTextOutput("combined_summary")
          )
        )
      )
    )
  ),

  footer = tags$footer(
    HTML('© 2025 <a href="https://github.com/yakubinaweed/refineR-reference-interval" target="_blank">Naweed Yakubi</a> • All rights reserved.'),
    style = "
      position: bottom;
      bottom: 0;
      width: 100%;
      text-align: center;
      padding: 10px 0;
      color: #777;
      font-size: 0.8em;"
  )
)