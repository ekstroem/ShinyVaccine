# ShinyVaccine

Shiny app to show the SIR model with an additional vaccine compartment

You can run this locally on your computer with 

```{r}
shiny::runGitHub("ekstroem/ShinyVaccine")
```

it requires that you have the following packages installed:

```{r}
library("shiny")
library("deSolve")
library("cowplot")
library("ggplot2")
library("tidyverse")
library("ggrepel")
library("shinydashboard")
```