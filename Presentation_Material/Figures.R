################################# FRONT MATTER #################################

# Script to build figures and examples for the presentation

rm(list = ls())
gc()

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

path <- getwd()
save_path <- paste0(path,"/Presentation_Material/")

#################### WEAK AND STRONG MONOTONICITY EXAMPLES #####################

# Two covariates (X1, X2), plot regions where D(1) > D(0) and D(0) > D(1)

# Weak Monotonicity
compliers <- data.frame(x = c(0.261, 0.673 , 1,  1),
                        y = c(1, 0, 0, 1))
defiers <- data.frame(x = c(0, 0.673, 0.261, 0),
                      y = c(0, 0, 1, 1))
e <- 0.005 # controls the width of the boundary line
boundary <- data.frame(x = c(0.261 - e, 0.261 + e, 0.673 + e, 0.673 - e),
                       y = c(1, 1, 0, 0))

pdf(paste0(save_path, "Weak_Montonicity_Example.pdf"), width = 4.43, height = 3.92)
print(ggplot() + 
        geom_polygon(aes(x = x, y = y, col= "defiers", fill = "defiers"), 
                     data = defiers) + 
        geom_polygon(aes(x = x, y = y, col = "compliers", fill = "compliers"),
                     data = compliers) +
        geom_polygon(aes(x = x, y = y, col = "boundary", fill = "boundary"), 
                     data = boundary) + 
        annotate(geom = "text", x = 0.25, y = 0.45, label = "(Defiers)") + 
        annotate(geom = "text", x = 0.25, y = 0.51, label = expression("D(0)">="D(1)")) + 
        annotate(geom = "text", x = 0.7, y = 0.55, label = "(Compliers)") +
        annotate(geom = "text", x = 0.7, y = 0.61, label = expression("D(1)">="D(0)")) + 
        scale_fill_manual(values = c("defiers" = "#F8766D", 
                                     "compliers" = "#00BFC4",
                                     "boundary" = "black")) + 
        scale_color_manual(values = c("defiers" = "#F8766D",
                                      "compliers" = "#00BFC4",
                                      "boundary" = "black")) + 
        theme_classic() + 
        labs(x = expression("x"[1]), y = expression("x"[2])) + 
        theme(legend.position = "null",
              # axis.ticks.x = element_blank(),
              # axis.ticks.y = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank()))  
dev.off()
  
# Strong Monotonicity
compliers <- data.frame(x = c(0, 0, 1, 1),
                        y = c(1, 0, 0, 1))

pdf(paste0(save_path, "Strong_Montonicity_Example.pdf"), width = 4.43, height = 3.92)
print(ggplot() + 
        geom_polygon(aes(x = x, y = y, col = "compliers", fill = "compliers"),
                     data = compliers) +
        annotate(geom = "text", x = 0.5, y = 0.48, label = "(Compliers)") +
        annotate(geom = "text", x = 0.5, y = 0.54, label = expression("D(1)">="D(0)")) + 
        scale_fill_manual(values = c("compliers" = "#00BFC4")) + 
        scale_color_manual(values = c("compliers" = "#00BFC4")) + 
        theme_classic() + 
        labs(x = expression("x"[1]), y = expression("x"[2])) + 
        theme(legend.position = "null",
              # axis.ticks.x = element_blank(),
              # axis.ticks.y = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank())) 
dev.off()


############################# REORDERED IV EXAMPLE #############################

# Redo the weak monotonicity plot with a change in colours and labels

compliers <- data.frame(x = c(0.261, 0.673 , 1,  1),
                        y = c(1, 0, 0, 1))
defiers <- data.frame(x = c(0, 0.673, 0.261, 0),
                      y = c(0, 0, 1, 1))
e <- 0.005 # controls the width of the boundary line
boundary <- data.frame(x = c(0.261 - e, 0.261 + e, 0.673 + e, 0.673 - e),
                       y = c(1, 1, 0, 0))

pdf(paste0(save_path, "Reordered_IV_Example.pdf"), width = 4.43, height = 3.92)
print(ggplot() + 
        geom_polygon(aes(x = x, y = y, col= "defiers", fill = "defiers"), 
                     data = defiers) + 
        geom_polygon(aes(x = x, y = y, col = "compliers", fill = "compliers"),
                     data = compliers) +
        geom_polygon(aes(x = x, y = y, col = "boundary", fill = "boundary"), 
                     data = boundary) + 
        annotate(geom = "text", x = 0.25, y = 0.45, label = "(Compliers for -Z)") + 
        annotate(geom = "text", x = 0.25, y = 0.51, label = expression("D(0)">="D(-1)")) + 
        annotate(geom = "text", x = 0.7, y = 0.55, label = "(Compliers for Z)") +
        annotate(geom = "text", x = 0.7, y = 0.61, label = expression("D(1)">="D(0)")) + 
        scale_fill_manual(values = c("defiers" = "#00BFC4", 
                                     "compliers" = "#00BFC4",
                                     "boundary" = "black")) + 
        scale_color_manual(values = c("defiers" = "#00BFC4",
                                      "compliers" = "#00BFC4",
                                      "boundary" = "black")) + 
        theme_classic() + 
        labs(x = expression("x"[1]), y = expression("x"[2])) + 
        theme(legend.position = "null",
              # axis.ticks.x = element_blank(),
              # axis.ticks.y = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank()))  
dev.off()

########################### VC SUBGRAPH ILLUSTRATION ###########################

# Illustration with a very complicated classification boundary to show what VC
# subgraph rules out
set.seed(1882)
e <- 0.005
compliers <- data.frame(x = c(-e, -e, 1 + e, 1 + e),
                        y = c(1 + e, -e, -e, 1 + e))
n <- 10001
defiers <- data.frame(x = runif(n, 0, 1), y = runif(n, 0, 1))

pdf(paste0(save_path, "VC_Example.pdf"), width = 4.43, height = 3.92)
print(ggplot() + 
        geom_polygon(aes(x = x, y = y, col = "compliers", fill = "compliers"),
                     data = compliers) +
        geom_point(aes(x = x, y = y, col= "defiers", fill = "defiers"), 
                   data = defiers, size = 0.5) + 
        scale_fill_manual(values = c("defiers" = "#F8766D", 
                                     "compliers" = "#00BFC4")) + 
        scale_color_manual(values = c("defiers" = "#F8766D",
                                      "compliers" = "#00BFC4")) + 
        theme_classic() + 
        labs(x = expression("x"[1]), y = expression("x"[2])) + 
        theme(legend.position = "null",
              # axis.ticks.x = element_blank(),
              # axis.ticks.y = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank())
      ) 
dev.off()

# A simple (but not linear) boundary to show what you can have

compliers <- data.frame(x = c(0, 0, 1, 1),
                        y = c(1, 0, 0, 1))

n <- 1001
x <- seq(0, 1, length.out = n)
poly_order <- 10
a <- runif(poly_order, -1, 1)
a <- a / sqrt(sum(a ^ 2))
y <- sapply(1:poly_order, function(x, k){return(x ^ k)}, x = x) %*% 
  matrix(a, nrow = poly_order, ncol = 1)

defiers <- data.frame(x = x, y = y + 0.9)
compliers_2 <- data.frame(x = x, y = y + 0.3)

pdf(paste0(save_path, "VC_Holds_Example.pdf"), width = 4.43, height = 3.92)
print(ggplot() + 
        geom_polygon(aes(x = x, y = y, col = "compliers", fill = "compliers"),
                     data = compliers) +
        geom_area(aes(x = x, y = y, col = "defiers", fill = "defiers"), 
                    data = defiers) + 
        geom_area(aes(x = x, y = y, col = "compliers2", fill = "compliers"),
                  data = compliers_2) + 
        scale_fill_manual(values = c("defiers" = "#F8766D", 
                                     "compliers" = "#00BFC4",
                                     "compliers2" = "#00BFC4")) + 
        scale_color_manual(values = c("defiers" = "black",
                                      "compliers" = "white",
                                      "compliers2" = "black")) + 
        annotate(geom = "text", x = 0.45, y = 0.6, label = "Defiers") +
        annotate(geom = "text", x = 0.35, y = 0.15, label = "Compliers") +
        annotate(geom = "text", x = 0.85, y = 0.8, label = "Compliers") + 
        ylim(0, 1) + 
        theme_classic() + 
        labs(x = expression("x"[1]), y = expression("x"[2])) + 
        theme(legend.position = "null",
              legend.title = element_blank(),
              # axis.ticks.x = element_blank(),
              # axis.ticks.y = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank())
) 
dev.off()

####################### UNIFORM CONVERGENCE ILLUSTRATION #######################

# Show an example where uniform convergence happens and where it doesn't
n <- 1001
x <- seq(0, 1, length.out = n)
B <- runif(1, -1, 1)
B0 <- runif(1, -1, 1)
B_hat <- B + rnorm(1, sd = 0.1)
B0_hat <- B0 + rnorm(1, sd = 0.1)

PD1 <- pnorm(B0 + x * B)
P_hat <- pnorm(B0_hat + x * B_hat)

e <- rnorm(n, sd = 0.1)
e <- e * rbinom(n, 1, 0.01)
rho <- 0.7
v <- rep(0, n)
#v[1] <- rnorm(1, mean = 0, sd = (1 / (1 - rho ^ 2)))
for (i in 2:n){
  v[i] <- rho * v[i - 1] + e[i]
}
P_hat_2 <- P_hat + v

plot_data <- data.frame(x = x,
                        P = PD1,
                        P_hat = P_hat,
                        P_hat_2 = P_hat_2)

pdf(paste0(save_path, "Uniform_Convergence.pdf"), width = 4.43, height = 3.92)
plot_data %>% ggplot() + 
  geom_line(aes(x = x, y = P, col = "true"), size = 1.5, linetype = "dashed") +
  geom_line(aes(x = x, y = P_hat, col = "estimate"), 
            size = 1.5) + 
  theme_bw() + 
  labs(x = "x", y = "") + 
  theme(legend.position = c(0.2, 0.9),
        axis.text.y = element_blank(),
        legend.title = element_blank()) +
  scale_color_discrete(labels = c(expression(hat("p")("x,z")),
                                  expression("p(x,z)")))
dev.off()


pdf(paste0(save_path, "Not_Uniform_Convergence.pdf"), width = 4.43, height = 3.9)
plot_data %>% ggplot() + 
  geom_line(aes(x = x, y = P, col = "true"), size = 1.5, linetype = "dashed") +
  geom_line(aes(x = x, y = P_hat_2, col = "estimate"),
            size = 1.5) + 
  theme_bw() + 
  labs(x = "x", y = "") + 
  theme(legend.position = c(0.2, 0.9),
        axis.text.y = element_blank(),
        legend.title = element_blank()) +
  scale_color_discrete(labels = c(expression(hat("p")("x,z")),
                                  expression("p(x,z)")))
dev.off()


######################### MARGIN CONDITION ILLUSTRATION ########################

# Illustration to show what the margin condition rules out
# Example process suggested by Claude
# Distribution: Pr(|X| \leq t) = 1 / (1 + ln(1 / t))
# Density: f(|X| = t) = 1 / t * (1 + ln(1 / t)) ^ 2 

# Illustration where the margin condition holds
# X is is normally distributed
# Distribution: Pr(|X| \leq t) = pnorm(t) - pnorm(-t)
# Density; f(|X| = t) = dnorm(-t) + dnorm(t)
# The margin condition is then given by dnorm(0) * t

n <- 1001
t <- seq(0, 0.1, length.out = n)
margin_fails <- data.frame(t = t,
                           dens = 1 / (t * (1 + log(1 / t)) ^ 2),
                           dist = 1 / (1 + log(1 / t)))
sd_Dp <- 0.1
margin_holds <- data.frame(t = t,
                           dens = dnorm(-t, sd = sd_Dp) + dnorm(t, sd = sd_Dp),
                           dist = pnorm(t, sd = sd_Dp) - pnorm(-t, sd = sd_Dp),
                           margin = dnorm(0, sd = sd_Dp) * 2 * t,
                           margin_dens = rep(2 * dnorm(0, sd = sd_Dp), n))



pdf(paste0(save_path, "Margin_Fails_Dens.pdf"), width = 4.43, height = 3.92)
margin_fails %>% ggplot() + 
  geom_line(aes(x = t, y = dens, col = "dens"), size = 1.5)  + 
  theme_bw() + 
  labs(x = "t", y = "") + 
  theme(legend.position = c(0.8, 0.9),
        axis.text.y = element_blank(),
        legend.title = element_blank()) +
  scale_color_discrete(labels = c("f(|p(x,1) - p(x,0)|)=t"))
dev.off()
 
pdf(paste0(save_path, "Margin_Fails_Dist.pdf"), width = 4.43, height = 3.92)
margin_fails %>% ggplot() + 
  geom_line(aes(x = t, y = dist, col = "dist"), size = 1.5) + 
  theme_bw() + 
  labs(x = "t", y = "") + 
  theme(legend.position = c(0.25, 0.9),
        axis.text.y = element_blank(),
        legend.title = element_blank()) +
  scale_color_discrete(labels = c(expression("Pr(|p(x,1) - p(x,0)|)"<="t")))
dev.off()

pdf(paste0(save_path, "Margin_Holds_Dens.pdf"), width = 4.43, height = 3.92)
margin_holds %>% ggplot() + 
  geom_line(aes(x = t, y = dens, col = "dens"), size = 1.5)  + 
  geom_line(aes(x = t, y = margin_dens, col = "margin_dens"), 
            size = 1.5, linetype = "dashed") + 
  theme_bw() + 
  labs(x = "t", y = "") + 
  theme(legend.position = c(0.25, 0.15),
        axis.text.y = element_blank(),
        legend.title = element_blank()) +
  scale_color_discrete(labels = c(expression("f(|p(x,1) - p(x,0)|)=t"), 
                                  expression("C")))
dev.off()

pdf(paste0(save_path, "Margin_Holds_Dist.pdf"), width = 4.43, height = 3.92)
margin_holds %>% ggplot() + 
  geom_line(aes(x = t, y = dist, col = "dist"), size = 1.5) + 
  geom_line(aes(x = t, y = margin, col = "margin"), 
            size = 1.5, linetype = "dashed") + 
  theme_bw() + 
  labs(x = "t", y = "") + 
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(0.25, 0.9)) +
  scale_color_discrete(labels = c(expression("Pr(|p(x,1) - p(x,0)|)"<="t"), 
                                  expression("Ct"^zeta)))
dev.off()

########################## SIMPLE EXAMPLE ILLUSTRATION #########################

# Plot out p(x,1) - p(x,0) as a function of x

n <- 101
x <- seq(-5, 5, length.out = n)
D_p <- x * 0.05
x_vals <- c(-1 * (1:4), 1:4)

plot_data <- data.frame(x = x, D_p = D_p)
lines <- data.frame(x = x_vals, 
                    ymin = pmin(x_vals * 0.05, 0), 
                    ymax = pmax(x_vals * 0.05, 0))
points <- data.frame(x = rep(x_vals, 2),
                     y = c(rep(0, 8), x_vals * 0.05))

pdf(paste0(save_path, "Simple_Example.pdf"), width = 4.43, height = 3.92)
plot_data %>% ggplot() + 
  geom_line(aes(x = x, y = D_p, col = "D_p")) +
  geom_hline(aes(yintercept = 0, col = "x_axis")) + 
  geom_linerange(aes(x = x, ymin = ymin, ymax = ymax, col = "lines"), 
                 data = lines, linetype = "dashed") +
  geom_point(aes(x = x, y = y, col = "points"), data = points) + 
  annotate(geom = "text", 
           x = x_vals, 
           y = c(rep(0.02, 4), rep(-0.02, 4)),
           label = x_vals) + 
  theme_bw() + 
  scale_color_manual(values = c("D_p" = "#00BFC4",
                                "x_axis" = "black",
                                "lines" = "#F8766D",
                                "points" = "#F8766D")) +
  labs(x = "x", y = "p(x,1) - p(x,0)") +
  theme(legend.position = "null")
dev.off()  
