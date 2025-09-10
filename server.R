# server.R
library(shiny)
library(readxl)
library(tidyverse)
library(mclust)
library(moments)
library(shinyjs)
library(car)
library(refineR)
library(shinyFiles)
library(shinyWidgets)
library(bslib)
library(ggplot2)

# Source the individual server modules
source("server_main.R")
source("server_gmm.R")
source("server_parallel.R")

server <- function(input, output, session) {

  # --- Reactive Values for State Management ---
  # These reactive values store and manage the state of the application across different tabs.
  data_reactive <- reactiveVal(NULL)
  gmm_uploaded_data_rv <- reactiveVal(NULL)
  gmm_processed_data_rv <- reactiveVal(NULL)
  gmm_transformation_details_rv <- reactiveVal(list(male_hgb_transformed = FALSE, female_hgb_transformed = FALSE))
  selected_dir_reactive <- reactiveVal(NULL)
  message_rv <- reactiveVal(list(type = "", text = ""))
  analysis_running_rv <- reactiveVal(FALSE)
  parallel_data_rv <- reactiveVal(NULL)
  parallel_results_rv <- reactiveVal(list())
  parallel_message_rv <- reactiveVal(list(type = "", text = ""))

  # Renders an alert-style message UI element based on a reactive value
  renderMessageUI <- function(rv) {
    renderUI({
      msg <- rv()
      if (is.null(msg) || msg$text == "") {
        return(NULL)
      }
      class_name <- switch(msg$type,
                           "error" = "alert alert-danger",
                           "success" = "alert alert-success",
                           "warning" = "alert alert-warning",
                           "info" = "alert alert-info",
                           "alert alert-secondary")
      div(class = class_name, msg$text)
    })
  }

  # Render message UIs for each tab
  output$app_message <- renderMessageUI(message_rv)
  output$main_message <- renderMessageUI(message_rv)
  output$parallel_message <- renderMessageUI(parallel_message_rv)

  # Observer that prevents switching tabs when an analysis is running
  observeEvent(input$tabs, {
    if (!analysis_running_rv()) {
      message_rv(list(type = "", text = ""))
    } else {
      message_rv(list(text = "Tab switch blocked: An analysis is currently running. Please wait or reset the analysis.", type = "warning"))
    }
  })

  # Event handler for when a blocked tab is clicked
  observeEvent(input$tab_switch_blocked, {
    if (analysis_running_rv()) {
      message_rv(list(text = "Cannot switch tabs while an analysis is running. Please wait or reset the analysis.", type = "warning"))
    }
  })

  # Call the server logic for each module
  mainServer(input, output, session, data_reactive, selected_dir_reactive, message_rv, analysis_running_rv)
  gmmServer(input, output, session, gmm_uploaded_data_rv, gmm_processed_data_rv, gmm_transformation_details_rv, message_rv, analysis_running_rv)
  parallelServer(input, output, session, parallel_data_rv, parallel_results_rv, parallel_message_rv, analysis_running_rv)
}