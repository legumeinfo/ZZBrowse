packages <- c("shiny", "plyr", "tools", "xtable", "devtools", "markdown", "jsonlite", "shinyjs", "rintrojs", "RCurl", "yaml", "future", "promises", "DT")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

if(!is.element('rCharts', installed.packages()[,1])){
  require(devtools)
  install_github('ramnathv/rCharts')
}

# Specify the port to match the reverse proxy that enables Broadcast Channel communication.
# Default host is 127.0.0.1 or shiny.host, and must be an IPv4, not localhost.
shiny::runApp(launch.browser = TRUE, port = 3838)
