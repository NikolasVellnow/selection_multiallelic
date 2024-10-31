rm(list=ls())
library(ggplot2)

columns <- c("freq", "coverage", "prob_detection")

freqs <- c(0.005, 0.001, 0.01, 0.05, 0.1)

num_freqs <- length(freqs)
num_x_values <- 500

#Create a Empty DataFrame with 0 rows and n columns
df = data.frame(matrix(nrow = num_freqs*num_x_values, ncol = length(columns))) 

# Assign column names
colnames(df) = columns

counter <- 0

for (freq in freqs){
  print(freq)
  for (x in seq(1:num_x_values)){
    counter <- counter +1
    df$freq[counter] <- freq
    print(x)
    df$coverage[counter] <- x
    prob <- 1- (1-freq)^x
    print(prob)
    df$prob_detection[counter] <- prob
  }
}

df$freq <- as.factor(df$freq)

p <- ggplot(df, aes(x = coverage, y = prob_detection)) +
  geom_line(aes(color = freq))

p




  