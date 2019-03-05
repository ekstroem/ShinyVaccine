# ShinyVaccine

Shiny app to show the SIR model with an additional vaccine compartment

You can run the **Danish** language version locally on your computer with 

```{r}
shiny::runGitHub("ekstroem/ShinyVaccine/app.R")
```

and the English version is called with

```{r}
shiny::runGitHub("ekstroem/ShinyVaccine/app2.R")
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