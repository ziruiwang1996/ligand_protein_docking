library('ggpubr')

kendallcor <- function(file){
  data <- read.csv(file)
  cleaned_data <- data[complete.cases(data),]
  test <- cor.test(cleaned_data$EC50..µM., cleaned_data$Computational.Binding.Energy,
                   method='kendall')
  test
  ggscatter(cleaned_data, x='EC50..µM.', y='Computational.Binding.Energy', add = 'reg.line',
            conf.int = TRUE, cor.coef = TRUE, cor.method ='kendall',
            xlab = 'EC50 (µM)', ylab='Computational Free Energy')
}

spearmancor <- function(file){
  data <- read.csv(file)
  cleaned_data <- data[complete.cases(data),]
  test <- cor.test(cleaned_data$EC50..µM., cleaned_data$Computational.Binding.Energy, method='spearman')
  test
}

pearsoncor <- function(file){
  data <- read.csv(file)
  cleaned_data <- data[complete.cases(data),]
  test <- cor.test(cleaned_data$EC50..µM., cleaned_data$Computational.Binding.Energy, method='pearson')
  test
}

