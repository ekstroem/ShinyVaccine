#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library("shiny")
library("deSolve")
library("cowplot")
library("ggplot2")
library("tidyverse")
library("ggrepel")
library("shinydashboard")

## Create an SIR function
sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dI <-  beta * S * I - gamma * I
    dR <-                 gamma * I
    dV <- 0
    return(list(c(dS, dI, dR, dV)))
  })
}

# Define UI for application that draws a histogram
ui <- dashboardPage(
  dashboardHeader(disable = TRUE),
  dashboardSidebar(
    sliderInput("connum",
      "Basic reproductive number (R0, # persons):",
      min = .5, max = 20, value = 5
    ),
    sliderInput("pinf",
      "# infected at outbreak:",
      min = 1, max = 50, value = 2
    ),
    sliderInput("pvac",
      "Proportion vaccinated / immune (%):",
      min = 0, max = 100, value = 75
    ),
    sliderInput("vaceff",
      "Vaccine effectiveness (%):",
      min = 0, max = 100, value = 85
    ),
    sliderInput("infper",
      "Infection period (days):",
      min = 1, max = 30, value = 7
    ),
    sliderInput("timeframe",
      "Time frame (days):",
      min = 1, max = 400, value = 200
    )
    
  ),
  dashboardBody(
    tags$head(tags$style(HTML('
                              /* body */
                              .content-wrapper, .right-side {
                              background-color: #fffff8;
                              }                              
                              '))),
        
    #    mainPanel(
    fluidRow(plotOutput("distPlot")),
    br(),
    fluidRow(
      # Dynamic valueBoxes
      valueBoxOutput("progressBox", width = 6),
      valueBoxOutput("approvalBox", width = 6),
      valueBoxOutput("BRRBox", width = 6),
      valueBoxOutput("HIBox", width = 6)
    ),
    br(),
    br()
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # Create reactive input
  popsize <- 5700000
  dataInput <- reactive({
    init       <-
      c(
        S = 1 - input$pinf / popsize - input$pvac / 100 * input$vaceff / 100,
        I = input$pinf / popsize,
        R = 0,
        V = input$pvac / 100 * input$vaceff / 100
      )
    ## beta: infection parameter; gamma: recovery parameter
    parameters <-
      c(beta = input$connum * 1 / input$infper,
        # * (1 - input$pvac/100*input$vaceff/100),
        gamma = 1 / input$infper)
    ## Time frame
    times      <- seq(0, input$timeframe, by = .2)
    
    ## Solve using ode (General Solver for Ordinary Differential Equations)
    out <- ode(
        y = init,
        times = times,
        func = sir,
        parms = parameters
      )   
    #    out
    as.data.frame(out)
  })
  
  output$distPlot <- renderPlot({
    out <-
      dataInput() %>%
      gather(key, value, -time) %>%
      mutate(
        id = row_number(),
        key2 = recode(
          key,
          S = "Susceptible (S)",
          I = "Infected (I)",
          R = "Recovered (R)",
          V = "Vaccinated / immune (V)"
        ),
        keyleft = recode(
          key,
          S = "Susceptible (S)",
          I = "",
          R = "",
          V = "Vaccinated / immune (V)"
        ),
        keyright = recode(
          key,
          S = "",
          I = "Infected (I)",
          R = "Recovered (R)",
          V = ""
        )
      )
    
    ggplot(data = out,
           aes(
             x = time,
             y = value,
             group = key2,
             col = key2,
             label = key2,
             data_id = id
           )) + # ylim(0, 1) +
      ylab("Proportion of full population") + xlab("Time (days)") +
      geom_line(size = 2) +
      geom_text_repel(
        data = subset(out, time == max(time)),
        aes(label = keyright),
        size = 6,
        segment.size  = 0.2,
        segment.color = "grey50",
        nudge_x = 0,
        hjust = 1,
        direction = "y"
      ) +
      geom_text_repel(
        data = subset(out, time == min(time)),
        aes(label = keyleft),
        size = 6,
        segment.size  = 0.2,
        segment.color = "grey50",
        nudge_x = 0,
        hjust = 0,
        direction = "y"
      ) +
      theme(legend.position = "none") +
      scale_colour_manual(values = c("red", "green4", "black", "blue")) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      theme(
        rect=element_rect(size=0),
        legend.position="none",
        panel.background=element_rect(fill="transparent", colour=NA),
        plot.background=element_rect(fill="transparent", colour=NA),
        legend.key = element_rect(fill = "transparent", colour = "transparent")
      )
    
  })
  
  output$progressBox <- renderValueBox({
    valueBox(
      dataInput() %>% filter(time == max(time)) %>% select(R) %>% mutate(R = round(100 * R, 2)) %>% paste0("%"),
      "Proportion of full population that got the disease by end of time frame",
      icon = icon("thumbs-up", lib = "glyphicon"),
      color = "black"
    )
  })
  
  output$approvalBox <- renderValueBox({
    valueBox(
      paste0(round(
        100 * (dataInput() %>% filter(row_number() == n()) %>% mutate(res = (R + I) / (S + I + R)) %>% pull("res")), 2), "%"),
      "Proportion of susceptibles that will get the disease by end of time frame",
      icon = icon("thermometer-full"),
      color = "black"
    )
  })
  
  output$BRRBox <- renderValueBox({
    valueBox(
      paste0(round(input$connum *
                     (1 - input$pvac / 100 * input$vaceff / 100), 2), ""),
      "Effective R0 (for populationen at outbreak, when immunity is taken into account)",
      icon = icon("arrows-alt"),
      color = "red"
    )
  })
  
  output$HIBox <- renderValueBox({
    valueBox(
      paste0(round(100 * (1 - 1 / (input$connum)), 2), "%"),
      "Proportion of population that needs to be immune for herd immunity",
      icon = icon("medkit"),
      color = "blue"
    )
  })  
}

# Run the application
shinyApp(ui = ui, server = server)
