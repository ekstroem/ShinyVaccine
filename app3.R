# This is a Shiny web application to show a simple SIR model with vaccination. 
# 
# Created by Claus Ekstrøm 2020
#
# You can run the application by clicking
# the 'Run App' button above.
#



#
# First load the necessary packages
#

library("shiny")
library("deSolve")
library("cowplot")
library("ggplot2")
library("tidyverse")
library("ggrepel")
library("shinydashboard")

#
# Create an SIR function for ude with the ordinary differential equations later on
#

seiqr <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * (I + E)
    dE <-  beta * S * (I + E) - (alpha + kappa + delta2)*E
    dI <-  alpha*E - (gamma + delta) * I
    dQ <-  delta * I + delta2*E - gamma * Q
    dR <-  kappa*E + gamma * (I + Q)
    dV <-  0
    return(list(c(dS, dE, dI, dQ, dR, dV)))
  })
}

#
# Define the UI
#

ui <- dashboardPage(
  dashboardHeader(disable = TRUE),
  dashboardSidebar(
  sidebarMenu(
      menuItem("Population", tabName = "menu_1",
    sliderInput(
      "connum",
      "Basic reproductive number (R0, antal personer):",
      min = .5,
      max = 20,
      value = 2.5
    ),
    sliderInput(
      "pinf",
      "Antal inficerede ved start:",
      min = 1,
      max = 50,
      value = 1
    ),
    sliderInput(
      "pvac",
      "Andel naturligt immune (%):",
      min = 0,
      max = 100,
      value = 0
    )
),
      menuItem("Sygdomskarakteristika", tabName = "menu_3",
    sliderInput(
      "infper",
      "Inkubationsperiode (dage):",
      min = 1,
      max = 30,
      value = 14
    ),
    sliderInput(
      "infper2",
      "Sygdomsperiode (dage):",
      min = 1,
      max = 30,
      value = 7
    )
),
      menuItem("Karantæne", tabName = "menu_2",
        sliderInput(
          "q1",
          "Karantænetid for inficerede (dage):",
          min = 0,
          max = 28,
          value = 0
        ),
        sliderInput(
          "q2",
          "Karantænetid for bærere (dage):",
          min = 0,
          max = 28,
          value = 0
        )
      ),
    sliderInput(
      "timeframe",
      "Tidsperiode (dage):",
      min = 1,
      max = 1200,
      value = 400
    )
)
  ),
  dashboardBody(
    tags$head(tags$style(HTML('
                              /* main sidebar 
                              .skin-blue .main-sidebar {
                              background-color: #808080;
                              } */

                              /* body */
                              .content-wrapper, .right-side {
                              background-color: #fffff8;
                              }                              
                              '))
),
        
    fluidRow(plotOutput("distPlot")),
    br(),
    fluidRow(
      # Dynamic valueBoxes
      valueBoxOutput("progressBox", width = 4),
      valueBoxOutput("approvalBox", width = 4),
      valueBoxOutput("BRRBox", width = 4)
#      valueBoxOutput("HIBox", width = 6)
      
    ),
    br(),
    br()
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # Create reactive input
  popsize <- 5700000   # Set a rough population size for Denmark
  dataInput <- reactive({
    init       <-
      c(
        S = 1 - input$pinf / popsize - input$pvac / 100,
        E = input$pinf / popsize,
        I = 0,
        Q = 0,
        R = 0,
        V = input$pvac / 100 
      )

    ## Time frame
    times      <- seq(0, input$timeframe, by = .2)

    ## beta: infection parameter; gamma: recovery parameter
    ## alpha + kappa er dage som exposed = 1/14
    parameters <-
      c(beta   = input$connum * 1 / (input$infper + input$infper2),   # Smitte
        gamma  = 1 / input$infper2,  # Dage 
        alpha  = 1 / input$infper , # 
        delta  = ifelse(input$q1 == 0, 0, 1 / input$q1), # quarraintaine rate
        delta2 = ifelse(input$q2 == 0, 0, 1 / input$q2), # quarraintaine rate
        kappa  = 0/28)   # Fra exposed -> R

    
    ## Solve using ode (General Solver for Ordinary Differential Equations)
    out <-
      ode(
        y = init,
        times = times,
        func = seiqr,
        parms = parameters
      )
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
          S = "Modtagelige (S)",
          E = "Bærere (E)",
          I = "Inficerede (I)",
          Q = "Karantæne (Q)",
          R = "Sygdomsramte (R)",
          V = "Naturligt immune (V)"
        ),
        keyleft = recode(
          key,
          S = "Modtagelige (S)",
          E = "Bærere (E)",
          I = "",
          R = "",
          Q = "",
          V = "Naturligt immune (V)"
        ),
        keyright = recode(
          key,
          S = "",
          I = "Inficerede (I)",
          R = "Sygdomsramte (R)",
          Q = "Karantæne (Q)",
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
      ylab("Procent af hele befolkningen") + xlab("Tid (dage)") +
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
      scale_colour_manual(values = c("khaki3", "red", "orange", "green4", "blue", "black")) +
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
      dataInput() %>% filter(time == max(time)) %>% dplyr::select(R) %>% mutate(R = round(100 *
                                                                                     R, 2)) %>% paste0("%"),
      paste0("Andel af populationen, der har haft virussen efter ", input$timeframe, " dage", collapse=""),
      color = "black",
      icon = icon("thermometer-full")
    )
  })
  
  output$approvalBox <- renderValueBox({
    valueBox(
      paste0(round(
        ifelse(input$q1==0, 0, 100 * ( 1/input$q1 / (1/input$infper2 + 1/input$q1))), 2
      ), "%"),
      "Andel af de inficerede, der kommer i karantæne",
      icon = icon("hospital"),
      color = "red"
    )
  })
  
  output$BRRBox <- renderValueBox({
    valueBox(
      paste0(round(
        ifelse(input$q2==0, 0, 100 * ( 1/input$q2 / (1/input$infper + 1/input$q2))), 2
      ), "%"),
      "Andel af bærere, der kommer i karantæne",
      icon = icon("hospital"),
      color = "yellow"
    )
  })  
}

# Run the application
shinyApp(ui = ui, server = server)
