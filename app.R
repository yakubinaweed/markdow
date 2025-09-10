# app.R

# Source the user interface (UI) and server logic from separate files
source("ui.R")
source("server.R")

# Run the Shiny application
shinyApp(ui = ui, server = server)