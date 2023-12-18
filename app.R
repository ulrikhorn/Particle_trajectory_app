library(shiny)
library(tidyverse)
library(magrittr)
library(data.table)
library(purrr)

#flow rate to velocity as a function of the tube diameter
lpm2ms <- function(lpm, diameter) {
  meters_m <- (4 * lpm) / (pi * (diameter ^ 2))
  meters_m / 60 * 10
}

# Relaxation time accoridng to stokes
relaxation_time <- function(vi, pdim) {
  scf = 1 + 0.066 / pdim * (2.514 + 0.8 * exp(-0.55 * pdim / 0.066))
  relax = (0.01 * (pdim ^ 2) * scf) / (18 * 0.000182)
  stop = (vi * relax) / 1000000
  relax2 = stop / vi
  relax2
}

# Particle stopping distance
stopping <- function(vi, pdim) {
  scf = 1 + 0.066 / pdim * (2.514 + 0.8 * exp(-0.55 * pdim / 0.066))
  relax = (0.01 * (pdim ^ 2) * scf) / (18 * 0.000182)
  stop = (vi * relax) / 1000000
  relax2 = stop / vi
  stop * 100
}

#Particle settling velocity
settling_velocity <- function(pdim) {
  scf <- 1 + 0.066 / pdim * (2.514 + 0.8 * exp(-0.55 * pdim / 0.066))
  settling <- (1 * (pdim ^ 2) * 0.00000981 * scf) / (18 * 0.000180711) #cm/s
  settling * 0.1
}

# velocity as a function of time and particle size
velocity_as_f_of_time <- function(pdim1, v0, t) {
  scf <- 1 + 0.066 / pdim1 * (2.514 + 0.8 * exp(-0.55 * pdim1 / 0.066))
  relax = relaxation_time(pdim = pdim1, vi = v0)
  v0 * exp(-t / relax)
}

# travel distance as a function of time
travel_dist_as_f_of_time <- function(pdim1, v0, t) {
  scf <- 1 + 0.066 / pdim1 * (2.514 + 0.8 * exp(-0.55 * pdim1 / 0.066))
  relax = relaxation_time(pdim = pdim1, vi = v0)
  v0 * relax * (1 - exp(-t / relax)) * 100
}

#function to discribe needed concentration for desired drydim
needed_concentration <- function(wetdim, drydim) {

  radius_m_wet <- (wetdim / 2) / (1 * 10 ^ 6) #m
  radius_m_dry <- (drydim / 2) / (1 * 10 ^ 6) #m
  
  volume_wet <- (4 / 3) * pi * (radius_m_wet ^ 3)
  volume_dry <- (4 / 3) * pi * (radius_m_dry ^ 3)
  paste(round(volume_dry / volume_wet, digits = 4), "g/l")
}

#Relative humidity function
rh_function <- function(mass_water, volume_of_space, temp) {
  saturated_water_vapor_pressure = 3.17 # kPa at 25C
  r = 8.314 #j mol-1K-1 gas constant
  n = mass_water / 18 # 18 = g/mol of water, n = mol of water
  v = volume_of_space #m3
  t =  temp #kelvin
  p = (n * r * t) / v #Pa
  kp = p / 1000
  
  RH = (kp / saturated_water_vapor_pressure) * 100
  RH
  
}

# Function to calculate the ml/min needed to generate a certain number of PPL
needed_mlmin <-
  function(wanted_ppl,
           flowrate,
           initial_particle_diameter) {
    droplet_size_um = initial_particle_diameter
    droplet_size_cm = droplet_size_um * 0.0001
    droplet_volume_cm3 = (4 / 3) * pi * (droplet_size_cm / 2) ^ 3
    droplets_per_ml = 1 / droplet_volume_cm3
    
    mlmin = (wanted_ppl * flowrate) / droplets_per_ml
    mlmin
  }

# Evaporation time function
evap_time <-
  function(temp = 293.15,
           particle_density = 1,
           initial_particle_diameter,
           vapor_diff_coeff = 0.24,
           RH = 0,
           vapor_mol_weight = 18) {
    saturation_vapor_pressure = exp(18.72 - 4062 / (temp - 37))
    surface_temp = temp + (6.65 + 0.345 * (temp - 273.15) + 0.0031 * (temp -
                                                                        273.15) ^ 2) * (RH - 1) / (1 + (0.082 + 0.00782 * (temp - 273.15)) * RH)
    
    vapor_press_surface = exp(18.72 - 4062 / (surface_temp - 37))
    
    evap_time = (83140000 * particle_density * (initial_particle_diameter *
                                                  0.0001) ^ 2) / (
                                                    8 * vapor_diff_coeff * vapor_mol_weight * (
                                                      vapor_press_surface * 1333 / surface_temp - saturation_vapor_pressure *
                                                        RH * 1333 / temp
                                                    )
                                                  )
    evap_time
  }


# Particle diameter as a function of evaporation time
pm_as_f_of_evap_time <- function(temp = 293.15,
                                 particle_density = 1,
                                 evap_time,
                                 vapor_diff_coeff = 0.24,
                                 RH = 0,
                                 vapor_mol_weight = 18) {
  saturation_vapor_pressure = exp(18.72 - 4062 / (temp - 37))
  surface_temp = temp + (6.65 + 0.345 * (temp - 273.15) + 0.0031 * (temp -
                                                                      273.15) ^ 2) * (RH - 1) / (1 + (0.082 + 0.00782 * (temp - 273.15)) * RH)
  
  vapor_press_surface = exp(18.72 - 4062 / (surface_temp - 37))
  
  
  PM = sqrt((evap_time * (
    8 * vapor_diff_coeff * vapor_mol_weight * (
      vapor_press_surface * 1333 / surface_temp - saturation_vapor_pressure *
        RH * 1333 / temp
    )
  )) / (83140000 * particle_density)) / 0.0001
  
  PM
}


# Define server logic required to draw a histogram
server <- function(input, output) {
  plot_reactive <- reactive({

    initial_velocity = lpm2ms(lpm = input$n2_lpm, diameter = 0.45)
    # Defining ml/min 
    mlmin1 <-
      needed_mlmin(
        wanted_ppl = input$ppl,
        flowrate = input$flow_lpm1+input$n2_lpm,
        initial_particle_diameter = input$initial_pm
      )
    #Defining system RH
    rh1 = rh_function(
      mass_water = mlmin1 - needed_mlmin(
        input$ppl,
        flowrate = input$flow_lpm1+input$n2_lpm,
        initial_particle_diameter = input$drydim
      ), # - the ml of dried particles that take up volume in a sphere
      volume_of_space = (input$rhflow * 0.001),
      temp = 293.15
    ) / 100
    
    # Making the dataframe
    df <- tibble(time = seq(from = 0, to = 15, by = 0.001)) %>%
      mutate("pm" = input$initial_pm) %>%
      mutate("velocity" = velocity_as_f_of_time(pdim1 = pm, v0 = initial_velocity, time)) %>%
      mutate("flow_lpm" = input$flow_lpm1+input$n2_lpm) %>%
      mutate("flow_velocity" = lpm2ms(flow_lpm, diameter = 10)) %>%
      mutate("mlmin" = mlmin1) %>%
      add_column("RH" = rh1) %>% 
      mutate("um_per_second" = pm / evap_time(initial_particle_diameter =  pm, RH = rh1)) %>% 
      mutate("um_per_millisec" = um_per_second * 0.001)
    
    # Calculating the particle diameter over time
    newpm = max(df$pm)
    for (i in 2:length(rev(df$pm))) {
      newpm[i] <- newpm[i - 1] - df$um_per_millisec[i - 1]
    }
    
    # Defining the time it takes for the particle to reach the desired dry dimention
    # Then defining the particle size to be this dimention after it has reached this time
    length_finaldim <- length(newpm[newpm > input$drydim])
    newpm[newpm < input$drydim] = input$drydim
    
    # Adding new variables to the dataframe
    df <- df %>% mutate("pm" = newpm) %>%
      mutate("settling_velocity" = settling_velocity(pm) * 0.1)
    
    # Make two seperate dataframes where one encapsulates the time before particle velocity matches the flow velocity
    # the second is the time when flow velocity is  greater than the particle velocity
    # Assume that the particle continues at the flow velocity
    df1 <- df %>% filter(velocity >= flow_velocity)
    df2 <-
      df %>% filter(velocity <= flow_velocity) %>% mutate("velocity" = flow_velocity)
    df3 <- bind_rows(df1, df2)
    
    # Defining weather the particle velocity overcomes the settling velocity based on its diameter
    df3 <- df3 %>%
      mutate("overcome" = ifelse(velocity > settling_velocity, "yes", "no")) %>%
      mutate(dist_traveled = c(0, rep(NA, times = nrow(df) - 1)))
    
    # Calculating the traveling distance of the particle based on velocity or the settling velocity
    # depending on weather this is overcome or not
    start <- df3$dist_traveled[1]
    for (i in 2:nrow(df3)) {
      if (df3$overcome[i - 1] == "yes") {
        start[i] <- start[i - 1] + df3$velocity[i - 1] * 0.1
      } else{
        start[i] <-
          start[i - 1] - (df3$settling_velocity[i - 1] * 0.1 - df3$flow_velocity[i -
                                                                                   1] * 0.1)
      }
    }
    
    df3 <- df3 %>% mutate(dist_traveled = start)
    
    under_40 <- df3$dist_traveled < 40
    if (min(df3$dist_traveled[2:nrow(df3)]) > 0 && df3$time[length_finaldim] < df3$time[length(under_40[under_40 == TRUE])]) {
      truefalse <- TRUE
    }else{
      truefalse <- FALSE
    }
    
    # protting the dataframe
    p <- ggplot(df3) +
      geom_point(mapping = aes(time, dist_traveled, color = overcome)) +
      labs(
        x = "Time (s)",
        y = "Distance traveled upwards in column (cm)",
        title = "Particle position in column",
        fvill = "Particle velocity > settling velocity",
        caption = paste("Is the particle alwarys above 0cm \n and dry before it reaches the top of the column? =", truefalse
      )) +
      
      geom_hline(yintercept = 40, color = "blue") +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = length_finaldim * 0.001, color = "red") +
      annotate(
        "text",
        x = length_finaldim * 0.001,
        y = 20,
        label = "Particle dry",
        angle = 90,
        vjust = -0.5,
        color = "red"
      ) +
      ylim(min(df3$dist_traveled), 50) +
      xlim(0, max(df3$time[df3$dist_traveled < 50])) +
      scale_color_manual(name = "Particle velocity > settling velocity",
                         values = c("yes" = "#00BFC4", "no" = "#F8766D"))
    
    

    
  })
  
  output$distPlot <- renderPlot({
    print(plot_reactive())
    
  })
  
  output$text2 <- renderText({
    print(
      "The overall goal of this plot is to simulate a particles position in a verticle column over time.
      Droplets are produced of a certain size (Initial PM) and start evaporating immidiatly. This evaporation affects the size of
      the relative humidity of the air which in turn affects the evaporation time of the particle. The size of the particle affects
      the time needed for the particle to adjust from the initial gas driven velocity to the constant upwards velocity of the air in the
      column (column m/s). "
    )
  })
  
  
  output$text1 <- renderText({
    paste(
      "Wet particle diameter = ",
      input$initial_pm,
      " \n",
      "Dry particle diameter=",
      input$drydim,
      " \n",
      "needed concentration = ",
      needed_concentration(wetdim = input$initial_pm, drydim = input$drydim),
      " \n",
      "needed ml/min = ",
      round(
        needed_mlmin(
          wanted_ppl = input$ppl,
          flowrate = input$flow_lpm1+input$n2_lpm,
          initial_particle_diameter = input$initial_pm
        ),
        digits = 4
      ),
      " \n",
      "Flow rate = ",
      input$flow_lpm1+input$n2_lpm,
      " \n",
      "Initial velocity =",
      round(lpm2ms(input$n2_lpm, diameter = 0.45), digits = 4),
      " \n",
      "Column m/s = ",
      round(lpm2ms(input$flow_lpm1+input$n2_lpm, diameter = 10), digits = 3),
      " \n",
      "RH = ",
      round(
        rh_function(
          mass_water = needed_mlmin(
            wanted_ppl = input$ppl,
            flowrate = input$flow_lpm1+input$n2_lpm,
            initial_particle_diameter = input$initial_pm
          ) - needed_mlmin(
            input$ppl,
            flowrate = input$flow_lpm1+input$n2_lpm,
            initial_particle_diameter = input$drydim
          ),
          volume_of_space = (input$rhflow * 0.001),
          temp = 293.15
        ),
        digits = 4
      ),
      " \n",
      "Evaporation time to 0 = ",
      round(
        evap_time(
          initial_particle_diameter = input$initial_pm,
          RH = rh_function(
            mass_water = needed_mlmin(
              wanted_ppl = input$ppl,
              flowrate = input$flow_lpm1+input$n2_lpm,
              initial_particle_diameter = input$initial_pm
            ),
            volume_of_space = (input$rhflow * 0.001),
            temp = 293.15
          ) / 100
        ),
        digits = 4
      ),
      " \n",
      "PPL at sampling point = ",
      input$ppl * (input$sample_dilution / 500),
      sep = ""
    )
  })
  
  
  output$data <- renderTable({
    mlmin1 <-
      needed_mlmin(
        wanted_ppl = input$ppl,
        flowrate = input$flow_lpm1+input$n2_lpm,
        initial_particle_diameter = input$initial_pm
      )
    rh1 = rh_function(
      mass_water = mlmin1 - needed_mlmin(
        input$ppl,
        flowrate = input$flow_lpm1+input$n2_lpm,
        initial_particle_diameter = input$drydim
      ),
      volume_of_space = (input$rhflow * 0.001),
      temp = 293.15
    ) / 100
    
    df <- tibble(time = seq(from = 0, to = 15, by = 0.001)) %>%
      #mutate("pm" = initial_pm) %>%
      mutate("pm" = input$initial_pm) %>%
      mutate("velocity" = velocity_as_f_of_time(pdim1 = pm, v0 = initial_velocity, time)) %>%
      mutate("flow_lpm" = input$flow_lpm1+input$n2_lpm) %>%
      mutate("flow_velocity" = lpm2ms(flow_lpm, diameter = 10)) %>%
      mutate("mlmin" = mlmin1) %>%
      add_column("RH" = rh1)
    
    df <-
      df %>% mutate("um_per_second" = pm / evap_time(initial_particle_diameter =  pm, RH = rh1)) %>% mutate("um_per_millisec" = um_per_second *
                                                                                                              0.001)
    
    newpm = max(df$pm)
    for (i in 2:length(rev(df$pm))) {
      newpm[i] <- newpm[i - 1] - df$um_per_millisec[i - 1]
    }
    
    length_finaldim <- length(newpm[newpm > input$drydim])
    newpm[newpm < input$drydim] = input$drydim
    
    df <- df %>% mutate("pm" = newpm) %>%
      mutate("settling_velocity" = settling_velocity(pm) * 0.1)
    
    
    # make two seperate dataframes where one encapsulates the time before particle velocity matches the flow velocity
    # the second is the time when flow velocity is  greater than the particle velocity
    # Assume that the particle continues at the flow velocity
    df1 <- df %>% filter(velocity >= flow_velocity)
    df2 <-
      df %>% filter(velocity <= flow_velocity) %>% mutate("velocity" = flow_velocity)
    df3 <- bind_rows(df1, df2)
    df3 <- df3 %>%
      mutate("overcome" = ifelse(velocity > settling_velocity, "yes", "no")) %>%
      mutate(dist_traveled = c(0, rep(NA, times = nrow(df) - 1)))
    
    
    start <- df3$dist_traveled[1]
    for (i in 2:nrow(df3)) {
      if (df3$overcome[i - 1] == "yes") {
        start[i] <- start[i - 1] + df3$velocity[i - 1] * 0.1
      } else{
        start[i] <- start[i - 1] - df3$settling_velocity[i - 1] * 0.1
      }
    }
    
    df3 <- df3 %>% mutate(dist_traveled = start)
    df3
  })
}

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("Particle position in spredningskolonne"),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        "initial_pm",
        "Wet particle diameter:",
        min = 10,
        max = 60,
        value = 55
      ),
      sliderInput(
        "drydim",
        "Dry particle diameter",
        min = 0,
        max = 50,
        value = 20
      ),
      sliderInput(
        "flow_lpm1",
        "Column flowrate",
        min = 0,
        max = 100,
        value = 25
      ),
      sliderInput(
        "n2_lpm",
        "N2 flowrate",
        min = 0,
        max = 30,
        value = 5
      ),
      numericInput(
        "ppl",
        "Wanted ppL",
        min = 0,
        max = 500000,
        value = 170000
      ),
      sliderInput(
        "rhflow",
        "Humidity flowrate (Should be same as column flowrate + N2 flowrate)",
        min = 0,
        max = 100,
        value = 30
      ),
      sliderInput(
        "sample_dilution",
        "LPM of sample air that ends up in horizontal column (500 LPM)",
        min = 0,
        max = 100,
        value = 5
      )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(tabsetPanel(
      tabPanel("plot", plotOutput("distPlot"), verbatimTextOutput("text1")),
      tabPanel("data", tableOutput("data"))
      # ,
      # tabPanel("Explanation", verbatimTextOutput("text2"))
      
    ))
  )
)

# Run the application
shinyApp(
  ui = ui,
  server = server,
  options = list(display.mode = "showcase")
)
