library(shiny)
library(ggplot2)
library(smatr)
library(shinyjs)

fluidPage(
  useShinyjs(),
  tags$head(
    tags$style(HTML("
                    p {
                    font-size:small
                    }
  "))
  ),
  
  # Application title
titlePanel("Line Fitting"),
  
tabsetPanel(
  #First Panel
   tabPanel("Population Scaling Relationship",
            # Explanatory Text
            fluidRow(
              column(12,
                     p(HTML("The population scaling relationship, based on 500 indivduals.
                            <span style='color: #00FF00;'> The green </span> line is the true scaling relationship (&mu;<sub>k<sub>y</sub></sub>/&mu;<sub>k<sub>x</sub></sub>).
                            The <span style='color: #FF0000;'> red line</span> is the OLS regression.
                            The <span style='color: #0000FF;'>blue line</span> is the MA regression.
                            The <span style='color: #FFA500;'>orange line</span> is the SMA regression.
                            Broken lines show the regression fitted to the population.
                            Solid lines show the true regression, based on the parameter values."))
              
            )
            ),
            # PLOT
            fluidRow(
              #Buffer
              column(2,
                     h4("")),
              #Plot
              column(8,
                     plotOutput("pop.scaling",width="600px",height="350px"))
              )),
   
   #Second Panel
   tabPanel("Best Line Fitting Method",
            # Explanatory Text
            fluidRow(
              column(12,
                     p(HTML("At each point the plot indicates the line-fitting method where the slope of the population scaling relationship
                              is closest to the true slope (&mu;<sub>k<sub>y</sub></sub>/&mu;<sub>k<sub>x</sub></sub>). 
                            Parameter combinations where OLS is the best fit are show in <span style='color: #FF0000;'>red</span>.
                            Parameter combinations where MA is the best fit are show in <span style='color: #0000FF;'>blue</span>.
                            Parameter combinations where SMA is the best fit are show in <span style='color: #FFA500;'>orange</span>.
                            The solid white lines indicate the currently assigned parameter values."))
                     )
                     ),
            # PLOT
            fluidRow(
              #Buffer
              column(4,
                     h4("")),
              #Plot
              column(6, 
                     selectInput("variable", NULL,c("\\u03C3i[x] v \\u03C3i[y]"="si1 v. si2", "\\u03C3k[x] v \\u03C3k[y]"="sk1 v. sk2")),
                     conditionalPanel(condition="input.variable=='si1 v. si2'", plotOutput("si1.v.si2.plot",width="400px",height="300px")),
                     conditionalPanel(condition="input.variable=='sk1 v. sk2'", plotOutput("sk1.v.sk2.plot",width="400px",height="300px"))
              )
            )),
   
   #Third Tab Panel
   tabPanel("Effect of Parameters on Slope",
            # Explanatory Text
            fluidRow(
              column(12,
                     p(HTML("The relationship between the slope of the population scaling relationship, as determined by <span style='color: #FF0000;'>OLS</span>, 
                              <span style='color: #0000FF;'>MA</span>, 
                              and <span style='color: #FFA500;'>SMA</span> regression. 
                              The true slope (&mu;<sub>k<sub>y</sub></sub>/&mu;<sub>k<sub>x</sub></sub>), is indicated by the <span style='color: #00FF00;'> 
                              green </span> line.
                            The solid black line indicates the currently assigned parameter value (e.g.&mu;<sub>k<sub>x</sub></sub>). 
                            The broken black line indicates the assigned vlaue of the same parameter for the other trait (e.g. &mu;<sub>k<sub>y</sub></sub>)."))
                     )),
            # PLOT
            fluidRow(
              #Buffer
              column(2,
                     h4("")),
              #Plot
              column(10,
                     selectInput("xaxis", NULL,c("\\u03BCS"='G',"\\u03C3S"='sG',"\\u03BCk[x]"="k1", "\\u03C3k[x]"="sk1","\\u03BCk[y]"="k2","\\u03C3k[y]"="sk2","\\u03C3i[x]"="si1", "\\u03C3i[y]"="si2")),
                     conditionalPanel(condition="input.xaxis=='k1'",plotOutput("slope.against.k1",width="500px",height="300px")),
                     conditionalPanel(condition="input.xaxis=='k2'",plotOutput("slope.against.k2",width="500px",height="300px")),
                     conditionalPanel(condition="input.xaxis=='si1'",plotOutput("slope.against.si1",width="500px",height="300px")),
                     conditionalPanel(condition="input.xaxis=='si2'",plotOutput("slope.against.si2",width="500px",height="300px")),
                     conditionalPanel(condition="input.xaxis=='sk1'",plotOutput("slope.against.sk1",width="500px",height="300px")),
                     conditionalPanel(condition="input.xaxis=='sk2'",plotOutput("slope.against.sk2",width="500px",height="300px")),
                     conditionalPanel(condition="input.xaxis=='sG'",plotOutput("slope.against.sG",width="500px",height="300px")),
                     conditionalPanel(condition="input.xaxis=='G'",plotOutput("slope.against.G",width="500px",height="300px"))
            ))
   )
   
),
  
  #hr(),
  
  # Input at the bottom
  fluidRow(
    column(3,
           numericInput("G",HTML("&mu;<sub>S</sub> : <small>mean level of systemic growth regulator</small>"), value = 0,step=0.1),
           numericInput("sG",HTML("&sigma;<sub>S</sub> : <small>variation in level of systemic growth regulator</small>"), value = 0.165,step=0.1),
           actionButton("resetAll","Reset Parameter Values")),
    column(4,
           numericInput("k1",HTML("&mu;<sub>k<sub>x</sub></sub> : <small>mean sensitivity to systemic growth regulator for <i>x trait</i></small>"), value = 0.842, step=.1),
           numericInput("sk1",HTML("&sigma;<sub>k<sub>x</sub></sub> : <small>variation in sensitivity to systemic growth regulator for <i>x trait</i></small>"), value = 0.605, step=.1),
           numericInput("k2",HTML("&mu;<sub>k<sub>y</sub></sub> : <small>mean sensitivity to systemic growth regulator for <i>y trait</i></small>"),value = 0.581,step=.1),
           numericInput("sk2",HTML("&sigma;<sub>k<sub>y</sub></sub> : <small>variation in sensitivity to systemic growth regulator for <i>y trait</i></small>"),value = 0.818,step=.1)),
    column(4,
           numericInput("i1",HTML("&mu;<sub>i<sub>x</sub></sub> : <small>mean trait-autonomous growth rate for <i>x trait</i></small>"), value = 13.545,step=.1),
           numericInput("si1",HTML("&sigma;<sub>i<sub>x</sub></sub> : <small>variation in trait-autonomous growth rate for <i>x trait</i></small>"),value = 0.105,step=.1),
           numericInput("i2",HTML("&mu;<sub>i<sub>y</sub></sub> : <small>mean trait-autonomous growth rate for <i>y trait</i></small>"), value = 13.945,step=.1),
           numericInput("si2",HTML("&sigma;<sub>i<sub>y</sub></sub> : <small>variation in trait-autonomous growth rate for <i>y trait</i></small>"), value = 0.109,step=.1))
  ))

