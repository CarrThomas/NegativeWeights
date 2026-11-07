################################# FRONT MATTER #################################

# Script to build figures and examples for the presentation

rm(list = ls())
gc()

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

########################## WEAK MONOTONICITY EXAMPLE ###########################

# Figure illustrating the difference between weak and strong monotonicity

# Two covariates (X1, X2), plot regions where D(1) > D(0) and D(0) > D(1)
# For strong monotonicity, there is a single region...

compliers <- data.frame(x = c(0.261, 0.673 , 1,  1),
                        y = c(1, 0, 0, 1))
defiers <- data.frame(x = c(0, 0.673, 0.261, 0),
                      y = c(0, 0, 1, 1))
boundary <- data.frame(x = c(0.261, 0.673),
                       y = c(1, 0))
ggplot() + 
  geom_polygon(aes(x = x, y = y, col= "defiers", fill = "defiers"), 
               data = defiers) + 
  geom_polygon(aes(x = x, y = y, col = "compliers", fill = "compliers"),
               data = compliers) +
  geom_line(aes(x = x, y = y, col = "boundary"), 
            data = boundary, linewidth = 1.25) + 
  annotate(geom = "text", x = 0.25, y = 0.45, label = "Compliers") + 
  annotate(geom = "text", x = 0.7, y = 0.55, label = "Defiers") +
  scale_fill_manual(values = c("defiers" = "#F8766D", 
                               "compliers" = "#00BFC4")) + 
  scale_color_manual(values = c("defiers" = "#F8766D",
                                "compliers" = "#00BFC4",
                                "boundary" = "black")) + 
  theme_classic() + 
  labs(x = expression("x"[1]), y = expression("x"[2])) + 
  theme(legend.position = "null")  

  
  
