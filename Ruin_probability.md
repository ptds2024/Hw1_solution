Ruin probability for Brownian motion
================

## Initialization of the parameters

``` r
T <- 1 # upper bound of the considered time frame
c <- 1 # drift rate
sigma <- 0.1 # votality rate
tau <- 0.001 # length of small segments
N_t <- ceiling(T/tau)+1 # number of small segments + 1
N <- 100000 # number of simulation
S <- matrix(nrow = N, ncol = N_t) # Brownian motion simulations
u <- 2.5 # threshold 
prob <- 0 # empirical ruin probability
step <- 1/N # contribution to empirical ruin probability
```

## Simulating trajectories

To save time for rerunning the same parts of the code we use
“cache=TRUE”. However, to prevent potential mistakes after changing
parameters we indicate that by “dependson =”parameters”.

``` r
# setup for generating random sequences

set.seed(1) 

for (j in 1:N){
  
# bool is a logical variable, which determines whether a j-th trajectory exceeds the threshold u:
  
bool <-FALSE 

# S[j,1] corresponds to the value of a geometric Brownian motion at t=0, which is equal to 1.

S[j, 1]=1 

# to generate a j-th trajectory, we need to create an auxiliary vector x with i.i.d normal components:

x <- rnorm(n = N_t-1, mean = 0, sd = sqrt(tau))

for (i in 2:N_t){
  
  # calculating the current simulation value
  
  S[j, i]=S[j, i-1]+c*S[j, i-1]*tau+sigma*S[j, i-1]*x[i-1] 
  
  # Below we check if a j-th trajectory exceeds a threshold u at   moment i:
  
  if ((S[j,i]>u) & (bool == FALSE)){ 
    
   # we remember that the j-th trajectory is crossing the threshold u
    
    bool = TRUE 
    
    # updating value of the empirical ruin probability
    
    prob = prob + step 
  }
}
}
print(prob) # printing the resulting value of prob
```

    ## [1] 0.79393

## Random choice of 10 simulations for plotting line graphs

``` r
# to create a nice plot we use ggplot2 library

library("ggplot2")
```

    ## Warning: package 'ggplot2' was built under R version 4.3.2

``` r
# Randomly select 10 trajectory indices and sorting them in the correct order(optional)

sam <- sample(1:N, size = 10, replace = FALSE) |> sort() 

# Next we extract the selected trajectories

S_sampled <- S[sam, ]  # This is a 10 x N_t matrix

# We will Create a data frame suitable for ggplot with columns: time, value, group

df_list <- list()  # Initialize a list to store individual data frames

for (i in 1:length(sam)) {
  traj_index <- sam[i]                # The index of the trajectory in S
  traj_values <- S_sampled[i, ]       # The values of the trajectory
  traj_time <- 1:N_t                  # Time steps
  traj_time_rescaled <- (traj_time - 1) / (N_t - 1)
  traj_group <- rep(traj_index, N_t)  # The group identifier for ggplot
  
  # Create a data frame for the current trajectory
  df_traj <- data.frame(
    time = traj_time_rescaled,
    value = traj_values,
    path_ID = as.factor(traj_index)
  )
  
  df_list[[i]] <- df_traj  # Add the data frame to the list
}

# Combine all individual data frames into one

df <- do.call(rbind, df_list)

# Plot using ggplot2

ggplot <- ggplot(df, aes(x = time, y = value, group = path_ID, color = path_ID)) +
  geom_line() +
  labs(title = "Randomly Selected Brownian Motion Trajectories",
       x = "Time Step",
       y = "Value") +
  geom_hline(yintercept = u)+
  theme_minimal()
ggplot + scale_x_continuous(limits = c(0, T))
```

![](Solution_files/figure-gfm/pressure-1.png)<!-- -->

From the presented plot, we can see that the ruin more likely happens
near the time moment T. In fact, it happens, since the variance of the
process is maximized at t=T.

## Calculating theoretical ruin probability

``` r
mu<- c-sigma^2/2

#prob_theor corresponds to theoretical value of the ruin probability

prob_theor <- pnorm(-(log(u)/(sigma*sqrt(T))-mu*sqrt(T)/sigma), mean = 0, sd = 1)+exp(-2*mu*log(u)/sigma^2)*pnorm(-(log(u)/(sigma*sqrt(T))+mu*sqrt(T)/sigma), mean = 0, sd = 1)

print(prob_theor) # printing the resulting value
```

    ## [1] 0.7843862

The values are relatively close to each other, but are slightly
different.
