
# Description:
# This R script includes functions for generating gamma-sar samples, calculating bias and mean squared error (MSE)
# of estimators, and generating plots to visualize the results. The code is designed for statistical analysis
# of several estimators under different sample sizes and replication scenarios.

# Functions:
# - generate_samples: Generates replicated gamma-sar samples with specified parameters.
# - calculate_bias_mse: Computes bias and MSE for estimators across various sample sizes and replications.
# - generate_plot: Generates bias and MSE plots for different mu values and estimators, organized in a grid layout.

# Usage:
# - Adjust parameters as needed and call the functions 

# Dependencies:
# - Ensure that the required packages (e.g., ggplot2) are installed for proper execution.

# Note:
# This code may produce warnings related to deprecated ggplot2 syntax, which can be safely ignored.

# source("../MainFunctions/gamma_sar_sample.r")
# source("../MainFunctions/entropy_gamma_sar.r")



generate_samples_gi0 <- function(sample_size, replication, mu, alpha, L) {
  samples <- vector("list", replication)
  for (r in 1:replication) {
    samples[[r]] <- gi0_sample(mu, alpha, L, sample_size)
  }
  return(samples)
}

calculate_entropy_and_test <- function(sample_sizes, R, B, mu, alpha, L, estimators) {
  true_entropy <- entropy_gamma_sar(L, mu)
  output <- data.frame(
    SampleSize = integer(0),
    Estimator = character(0),
    MeanEntropy = numeric(0),
    ZStatistic = numeric(0),
    PValue = numeric(0)
  )
  
  for (ssize in sample_sizes) {
    samples <- generate_samples_gi0(ssize, R, mu, alpha, L)
    
    for (estimator_name in names(estimators)) {
      estimator <- estimators[[estimator_name]]
      v.entropy <- numeric(R)
      
      for (r in 1:R) {
        sample <- samples[[r]]
        
        if (grepl(" ", estimator_name)) {
          v.entropy[r] <- estimator(sample, B = B)
        } else {
          v.entropy[r] <- estimator(sample)
        }
      }
      
      mean_entropy <- mean(v.entropy)
      
      z_statistic <- sqrt(R) * (mean_entropy - true_entropy) / sd(v.entropy)
      p_value <- 2 * (1 - pnorm(abs(z_statistic)))

      output <- rbind(
        output,
        data.frame(
          SampleSize = ssize,
          Estimator = estimator_name,
          MeanEntropy = round(mean_entropy, 5),
          ZStatistic = round(z_statistic, 5),
          PValue = round(p_value, 5)
        )
      )
    }
  }
  
  colnames(output) <- c("$n$", "Estimator", "Mean Entropy", "$Z$ Statistic", "$p$ Value")
  
  return(output)
}



# Function to calculate entropy and perform hypothesis testing
# calculate_entropy_and_test <- function(sample_sizes, R, B, mu, alpha, L, estimators) {
#   true_entropy <- entropy_gamma_sar(L, mu)
#   output <- data.frame(
#     SampleSize = integer(0),
#     Estimator = character(0),
#     MeanEntropy = numeric(0),
#     ZStatistic = numeric(0),
#     PValue = numeric(0)
#   )
#   
#   for (ssize in sample_sizes) {
#     samples <- generate_samples_gi0(ssize, R, mu, alpha, L)
#     
#     for (estimator_name in names(estimators)) {
#       estimator <- estimators[[estimator_name]]
#       v.entropy <- numeric(R)
#       
#       for (r in 1:R) {
#         sample <- samples[[r]]
#         
#         if (grepl(" ", estimator_name)) {
#           v.entropy[r] <- estimator(sample, B = B)
#         } else {
#           v.entropy[r] <- estimator(sample)
#         }
#       }
#       
#       mean_entropy <- mean(v.entropy)
#       
#       z_statistic <- sqrt(R) * (mean_entropy - true_entropy) / sd(v.entropy)
#       p_value <- 2 * (1 - pnorm(abs(z_statistic)))
#       # alpha_t = 0.05
#       # reject_null <- abs(z_statistic) > qnorm(1 - alpha_t / 2)
#       # 
#       # output <- rbind(
#       #   output,
#       #   data.frame(
#       #     SampleSize = ssize,
#       #     Estimator = estimator_name,
#       #     MeanEntropy = round(mean_entropy, 5),
#       #     ZStatistic = round(z_statistic, 5),
#       #     PValue = round(p_value, 5),
#       #     RejectNull = reject_null  
#       #   )
#       # )
#       # 
#       
#       output <- rbind(
#         output,
#         data.frame(
#           SampleSize = ssize,
#           Estimator = estimator_name,
#           MeanEntropy = round(mean_entropy, 5),
#           ZStatistic = round(z_statistic, 5),
#           PValue = round(p_value, 5)
#         )
#       )
#     }
#    }
#   
#   colnames(output) <- c("$n$", "Estimator", "Mean Entropy", "$Z$ Statistic", "$p$ Value")
#   
#   return(output)
# }
# 



calculate_bias_mse_gi0 <- function(sample_sizes, R,B, mu, alpha, L, estimators) {
  true_entropy <- entropy_gI0(mu, alpha, L)
  
  output <- data.frame(n = integer(0), Estimator = character(0), Bias = numeric(0), MSE = numeric(0))
  
  for (ssize in sample_sizes) {
    # Generate samples outside the loop
    samples <- generate_samples_gi0(ssize, R, mu, alpha, L)
    
    for (estimator_name in names(estimators)) {
      estimator <- estimators[[estimator_name]]
      v.entropy <- numeric(R)
      
      for (r in 1:R) {
        sample <- samples[[r]]
        
        if (grepl("Bootstrap", estimator_name)) {
          v.entropy[r] <- estimator(sample, B = B)
        } else {
          v.entropy[r] <- estimator(sample)
        }
      }
      
      mse <- mean((v.entropy - true_entropy)^2)
      bias <- mean(v.entropy) - true_entropy
      
      output <- rbind(output, data.frame(n = ssize, Estimator = estimator_name, Bias = round(bias, 5), MSE = round(mse, 5)))
    }
  }
  
  return(output)
}


# generate_plot_gi0_esp <- function(results_gi0, mu_values, selected_estimators, ncol = 2, nrow = 2) {
#  
#   plot_list <- list()
#   
#   for (mu_val in mu_values) {
#     
#     df <- results_gi0[[as.character(mu_val)]]
#     
#     
#     df_filtered <- df[df$Estimator %in% names(selected_estimators), ]
#     
#     
#     df_filtered$Estimator <- selected_estimators[df_filtered$Estimator]
#     
#     
#     df_filtered$Estimator <- as.character(df_filtered$Estimator)
#     
#     plot_bias <- ggplot(df_filtered, aes(x = n, y = Bias, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", linewidth = 0.5) +
#       labs(y = "Bias", x = expression(italic(n))) + #x = expression("Sample size"~(italic(n))))
#       #labs(y = "Bias", x = expression("Sample size"~"\\big("italic(n)"\\big)")) +
#        # scale_color_discrete(labels = TeX(df_filtered$Estimator)) +
#       scale_color_manual(values = pal_jama()(7)[1:6], labels = TeX(df_filtered$Estimator)) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
#       coord_cartesian(clip = 'off')+
#       theme(legend.position = "bottom", legend.box = "horizontal", legend.box.margin = margin(0, 0, 0, 0)) +
#       guides(color = guide_legend(nrow = 1))
#     
#     
#     # # 
#     # if (mu_val != mu_values[1]) {
#     #   plot_bias = plot_bias + theme(legend.position = "top")
#     # }
#     
#     plot_mse <- ggplot(df_filtered, aes(x = n, y = MSE, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", linewidth = 0.5) +
#       labs(y = "MSE", x = expression(italic(n))) +
#       #scale_color_discrete(labels = TeX(df_filtered$Estimator)) +
#       scale_color_manual(values = pal_jama()(7)[1:6], labels = TeX(df_filtered$Estimator)) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
#       coord_cartesian(clip = 'off')+
#       theme(legend.position = "bottom", legend.box = "horizontal", legend.box.margin = margin(0, 0, 0, 0)) +
#       guides(color = guide_legend(nrow = 1))
#     
#     
#     
#     # if (mu_val != mu_values[1]) {
#     #   plot_mse = plot_mse + theme(legend.position = "top")
#     # }
#     
#     
#     plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
#   }
#   
#   
#   # combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
#   #   plot_layout(guides = "collect")
#   combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
#     plot_layout(guides = "collect")
#  
#   return(combined_plot)
# }

# el anterior codigo genera escala automatica en el eje x

generate_plot_gi0_esp <- function(results_gi0, mu_values, selected_estimators, ncol = 2, nrow = 2) {
  
  plot_list <- list()
  
  for (mu_val in mu_values) {
    
    df <- results_gi0[[as.character(mu_val)]]
    
    df_filtered <- df[df$Estimator %in% names(selected_estimators), ]
    df_filtered$Estimator <- selected_estimators[df_filtered$Estimator]
    df_filtered$Estimator <- as.character(df_filtered$Estimator)
    
    plot_bias <- ggplot(df_filtered, aes(x = n, y = Bias, color = Estimator)) +
      geom_hline(yintercept = 0) +
      geom_point(size = 2) +
      geom_line(linetype = "solid", linewidth = 0.5) +
      labs(y = "Bias", x = expression(italic(n))) +
      scale_x_continuous(breaks = c(9, 25, 49, 81, 121)) +  # Añadir esta línea para escala personalizada
      scale_color_manual(values = pal_jama()(7)[1:6], labels = TeX(df_filtered$Estimator)) + # Eliminar labels y TeX
      annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
      coord_cartesian(clip = 'off') +
      theme(legend.position = "bottom", legend.box = "horizontal", legend.box.margin = margin(0, 0, 0, 0)) +
      guides(color = guide_legend(nrow = 1))
    
    plot_mse <- ggplot(df_filtered, aes(x = n, y = MSE, color = Estimator)) +
      geom_hline(yintercept = 0) +
      geom_point(size = 2) +
      geom_line(linetype = "solid", linewidth = 0.5) +
      labs(y = "MSE", x = expression(italic(n))) +
      scale_x_continuous(breaks = c(9, 25, 49, 81, 121)) +  # Añadir esta línea
      scale_color_manual(values = pal_jama()(7)[1:6], labels = TeX(df_filtered$Estimator)) + # Eliminar labels y TeX
      annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
      coord_cartesian(clip = 'off') +
      theme(legend.position = "bottom", legend.box = "horizontal", legend.box.margin = margin(0, 0, 0, 0)) +
      guides(color = guide_legend(nrow = 1))
    
    plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
  }
  
  combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
    plot_layout(guides = "collect")
  
  return(combined_plot)
}


# generate_plot_gi0_esp <- function(results_gi0, mu_values, selected_estimators, ncol = 2, nrow = 2) {
#   # Lista para almacenar los gráficos
#   plot_list <- list()
# 
#   for (mu_val in mu_values) {
#     # Obtener resultados para el valor actual de mu desde la lista
#     df <- results_gi0[[as.character(mu_val)]]
# 
#     # Filtrar el dataframe para incluir solo los estimadores seleccionados
#     df_filtered <- df[df$Estimator %in% names(selected_estimators), ]
# 
#     # Actualizar los nombres de los estimadores en notación LaTeX
#     df_filtered$Estimator <- selected_estimators[df_filtered$Estimator]
# 
#     plot_bias <- ggplot(df_filtered, aes(x = n, y = Bias, color = as.character(Estimator))) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", linewidth = 0.5) +
#       labs(y = "Bias", x = expression("Sample size"~(n))) +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
#       coord_cartesian(clip = 'off')
# 
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_bias = plot_bias + theme(legend.position = "top")
#     }
# 
#     plot_mse <- ggplot(df_filtered, aes(x = n, y = MSE, color = as.character(Estimator))) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", linewidth = 0.5) +
#       labs(y = "MSE", x = expression("Sample size"~(n))) +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
#       coord_cartesian(clip = 'off')
# 
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_mse = plot_mse + theme(legend.position = "top")
#     }
# 
#     # Agregar los gráficos a la lista
#     plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
#   }
# 
#   # Organizar los gráficos en una sola figura y añadir leyenda común en la parte inferior
#   combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
#     plot_layout(guides = "collect")
# 
#   # No mostrar la figura aquí, devolver el objeto combined_plot
#   return(combined_plot)
# }


# generate_plot_gi0_esp <- function(results_gi0, mu_values, selected_estimators, ncol = 2, nrow = 2) {
#   # Lista para almacenar los gráficos
#   plot_list <- list()
#   
#   for (mu_val in mu_values) {
#     # Obtener resultados para el valor actual de mu desde la lista
#     df <- results_gi0[[as.character(mu_val)]]
#     
#     # Filtrar el dataframe para incluir solo los estimadores seleccionados
#     df_filtered <- df[df$Estimator %in% selected_estimators, ]
#     
#     plot_bias <- ggplot(df_filtered, aes(x = n, y = Bias, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", linewidth = 0.5) +
#       labs(y = "Bias", x =  expression("Sample size"~(n))) +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
#       coord_cartesian(clip = 'off')
#     
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_bias = plot_bias + theme(legend.position = "top")
#     }
#     
#     plot_mse <- ggplot(df_filtered, aes(x = n, y = MSE, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", linewidth = 0.5) +
#       labs(y = "MSE", x =  expression("Sample size"~(n))) +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
#       coord_cartesian(clip = 'off')
#     
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_mse = plot_mse + theme(legend.position = "top")
#     }
#     
#     # Agregar los gráficos a la lista
#     plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
#   }
#   
#   # Organizar los gráficos en una sola figura y añadir leyenda común en la parte inferior
#   combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
#     plot_layout(guides = "collect")
#   
#   # No mostrar la figura aquí, devolver el objeto combined_plot
#   return(combined_plot)
# }

generate_plot_gi0_new <- function(results_gi0, mu_values, ncol = 2, nrow = 2) {
  # Lista para almacenar los gráficos
  plot_list <- list()
  
  for (mu_val in mu_values) {
    # Obtener resultados para el valor actual de mu desde la lista
    df <- results_gi0[[as.character(mu_val)]]
    
    plot_bias <- ggplot(df, aes(x = n, y = Bias, color = Estimator)) +
      geom_hline(yintercept = 0) +
      geom_point(size = 2) +
      geom_line(linetype = "solid", linewidth = 0.5) +
      labs(y = "Bias", x =  expression("Sample size"~(n))) +
      guides(color = guide_legend(title = "Estimator")) +
      annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
      coord_cartesian(clip = 'off')
    
    # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
    if (mu_val != mu_values[1]) {
      plot_bias = plot_bias + theme(legend.position = "top")
    }
    
    plot_mse <- ggplot(df, aes(x = n, y = MSE, color = Estimator)) +
      geom_hline(yintercept = 0) +
      geom_point(size = 2) +
      geom_line(linetype = "solid", linewidth = 0.5) +
      labs(y = "MSE", x =  expression("Sample size"~(n))) +
      guides(color = guide_legend(title = "Estimator")) +
      annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
      coord_cartesian(clip = 'off')
    
    # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
    if (mu_val != mu_values[1]) {
      plot_mse = plot_mse + theme(legend.position = "top")
    }
    
    # Agregar los gráficos a la lista
    plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
  }
  
  # Organizar los gráficos en una sola figura y añadir leyenda común en la parte inferior
  combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
    plot_layout(guides = "collect")
  
  # No mostrar la figura aquí, devolver el objeto combined_plot
  return(combined_plot)
}


generate_plot_gi0 <- function(sample_sizes, R, B, mu_values, alpha, L, estimators, ncol = 2, nrow = 2) {
  # Lista para almacenar los gráficos
  plot_list <- list()
  
  #
  for (mu_val in mu_values) {
    # Calcular resultados para el valor actual de mu
    results <- calculate_bias_mse_gi0 (sample_sizes, R, B, mu_val, alpha, L, estimators)
    #cat("Resultados para mu =", mu_val, "\n")
    
    
    
    df <- as.data.frame(results)
    
    
    plot_bias <- ggplot(df, aes(x = n, y = Bias, color = Estimator)) +
      geom_hline(yintercept = 0) +
      geom_point(size = 2) +
      geom_line(linetype = "solid", linewidth = 0.5) +
      labs(y = "Bias", x =  expression("Sample size"~(n))) +
      guides(color = guide_legend(title = "Estimator")) +
      annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
      coord_cartesian(clip = 'off')
    # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
    if (mu_val != mu_values[1]) {
      plot_bias = plot_bias + theme(legend.position = "top")
    }
    
    
    plot_mse <- ggplot(df, aes(x = n, y = MSE, color = Estimator)) +
      geom_hline(yintercept = 0) +
      geom_point(size = 2) +
      geom_line(linetype = "solid", linewidth = 0.5) +
      labs(y = "MSE", x =  expression("Sample size"~(n))) +
      guides(color = guide_legend(title = "Estimator")) +
      annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
      coord_cartesian(clip = 'off')
    
    
    # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
    if (mu_val != mu_values[1]) {
      plot_mse = plot_mse + theme(legend.position = "top")
    }
    
    # Agregar los gráficos a la lista
    plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
  }
  
  # Organizar los gráficos en una sola figura y añadir leyenda común en la parte inferior
  combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
    plot_layout(guides = "collect")
  
  # No mostrar la figura aquí, devolver el objeto combined_plot
  return(combined_plot)
}



#Gamma sar

generate_samples <- function(sample_size, replication, mu, L) {
  samples <- vector("list", replication)
  for (r in 1:replication) {
    samples[[r]] <- gamma_sar_sample(L, mu, sample_size)
  }
  return(samples)
}


calculate_variance <- function(sample_sizes, R, B, mu, L, estimators) {
  true_entropy <- entropy_gamma_sar(L, mu)
  
  
  output <- data.frame(SampleSize = integer(0), Estimator = character(0), Variance = numeric(0))
  
  for (ssize in sample_sizes) {
    samples <- generate_samples(ssize, R, mu, L)
    
    for (estimator_name in names(estimators)) {
      estimator <- estimators[[estimator_name]]
      v.entropy <- numeric(R)
      
      for (r in 1:R) {
        sample <- samples[[r]]
        
        if (grepl("Bootstrap", estimator_name)) {
          v.entropy[r] <- estimator(sample, B = B)
        } else {
          v.entropy[r] <- estimator(sample)
        }
      }
      
      variance <- var(v.entropy)
      
      output <- rbind(output, data.frame(SampleSize = ssize, Estimator = estimator_name, Variance = round(variance, 5)))
    }
  }
  
  return(output)
}

calculate_bias_mse <- function(sample_sizes, R, B, mu, L, estimators) {
  true_entropy <- entropy_gamma_sar(L, mu)
  
  output <- data.frame(n = integer(0), Estimator = character(0), Bias = numeric(0), MSE = numeric(0))
  
  for (ssize in sample_sizes) {
    # Generate samples outside the loop
    samples <- generate_samples(ssize, R, mu, L)
    
    for (estimator_name in names(estimators)) {
      estimator <- estimators[[estimator_name]]
      v.entropy <- numeric(R)
      
      for (r in 1:R) {
        sample <- samples[[r]]
        
        if (grepl("Bootstrap", estimator_name)) {
          v.entropy[r] <- estimator(sample, B = B)
        } else {
          v.entropy[r] <- estimator(sample)
        }
      }
      
      mse <- mean((v.entropy - true_entropy)^2)
      bias <- mean(v.entropy) - true_entropy
      
      output <- rbind(output, data.frame(n = ssize, Estimator = estimator_name, Bias = round(bias, 5), MSE = round(mse, 5)))
    }
  }
  
  return(output)
}



# 
generate_plot <- function(sample_sizes, R, B, mu_values, L, estimators, ncol = 2, nrow = 2) {
  # Lista para almacenar los gráficos
  plot_list <- list()

  #
  for (mu_val in mu_values) {
    # Calcular resultados para el valor actual de mu
    results <- calculate_bias_mse(sample_sizes, R, B, mu_val, L, estimators)
    #cat("Resultados para mu =", mu_val, "\n")



    df <- as.data.frame(results)


    plot_bias <- ggplot(df, aes(x = n, y = Bias, color = Estimator)) +
      geom_hline(yintercept = 0) +
      geom_point(size = 2) +
      geom_line(linetype = "solid", linewidth = 0.5) +
      labs(y = "Bias", x =  expression("Sample size"~(n))) +
      guides(color = guide_legend(title = "Estimator")) +
      annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
      coord_cartesian(clip = 'off')
    # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
    if (mu_val != mu_values[1]) {
      plot_bias = plot_bias + theme(legend.position = "top")
    }


    plot_mse <- ggplot(df, aes(x = n, y = MSE, color = Estimator)) +
      geom_hline(yintercept = 0) +
      geom_point(size = 2) +
      geom_line(linetype = "solid", linewidth = 0.5) +
      labs(y = "MSE", x =  expression("Sample size"~(n))) +
      guides(color = guide_legend(title = "Estimator")) +
      annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
      coord_cartesian(clip = 'off')


    # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
    if (mu_val != mu_values[1]) {
      plot_mse = plot_mse + theme(legend.position = "top")
    }

    # Agregar los gráficos a la lista
    plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
  }

  # Organizar los gráficos en una sola figura y añadir leyenda común en la parte inferior
  combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
    plot_layout(guides = "collect")

  # No mostrar la figura aquí, devolver el objeto combined_plot
  return(combined_plot)
}



# 
# generate_plot <- function(sample_sizes, R, B, mu_values, L, estimators, ncol = 2, nrow = 2) {
#   # Lista para almacenar los gráficos
#   plot_list <- list()
# 
#   for (mu_val in mu_values) {
#     # Calcular resultados para el valor actual de mu
#     results <- calculate_bias_mse(sample_sizes, R, B, mu_val, L, estimators)
#     df <- as.data.frame(results)
# 
#     # Dividir estimadores en dos grupos (sólidos y punteados)
#     solid_estimators <- names(estimators)[1:(length(estimators) / 2)]
#     dashed_estimators <- names(estimators)[-c(1:(length(estimators) / 2))]
# 
#     # Filtrar dataframes para cada tipo de línea
#     solid_data <- df[df$Estimator %in% solid_estimators, ]
#     dashed_data <- df[df$Estimator %in% dashed_estimators, ]
# 
#     # Plot Bias y MSE
#     plot_bias <- ggplot() +
#       geom_hline(yintercept = 0, linetype = "solid") +
#       geom_point(data = solid_data, aes(x = n, y = Bias, color = Estimator), size = 2) +
#       geom_line(data = solid_data, aes(x = n, y = Bias, color = Estimator), linetype = "solid", linewidth = 0.5) +
#       geom_point(data = dashed_data, aes(x = n, y = Bias, color = Estimator), size = 2) +
#       geom_line(data = dashed_data, aes(x = n, y = Bias, color = Estimator), linetype = "dashed", linewidth = 0.5) +
#       labs(y = "Bias", x = "Sample size") +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
#       coord_cartesian(clip = 'off') +
#       theme(legend.position = "top")
# 
#     plot_mse <- ggplot() +
#       geom_hline(yintercept = 0, linetype = "solid") +
#       geom_point(data = solid_data, aes(x = n, y = MSE, color = Estimator), size = 2) +
#       geom_line(data = solid_data, aes(x = n, y = MSE, color = Estimator), linetype = "solid", linewidth = 0.5) +
#       geom_point(data = dashed_data, aes(x = n, y = MSE, color = Estimator), size = 2) +
#       geom_line(data = dashed_data, aes(x = n, y = MSE, color = Estimator), linetype = "dashed", linewidth = 0.5) +
#       labs(y = "MSE", x = "Sample size") +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
#       coord_cartesian(clip = 'off') +
#       theme(legend.position = "top")
# 
#     plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
#   }
# 
#   combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
#     plot_layout(guides = "collect")
# 
#   return(combined_plot)
# }
# 




# 
generate_table <- function(sample_sizes, R, B, mu_values, L, estimators) {
  # Lista para almacenar las tablas
  table_list <- list()
  
  #
  for (mu_val in mu_values) {
    # Calcular resultados para el valor actual de mu
    results <- calculate_bias_mse(sample_sizes, R, B, mu_val, L, estimators)
    #cat("Resultados para mu =", mu_val, "\n")
    
    # Imprimir la tabla de resultados
    table_list[[as.character(mu_val)]] <- kable(
      results,
      caption = paste("Table for mu =", mu_val),
      format = "latex", 
      booktabs = TRUE
    ) %>%
      kable_styling(latex_options = c("striped", "hold_position"), font_size = 11)
  }
  
  # Return the list of tables
  return(table_list)
}

generate_table2 <- function(sample_sizes, R, B, mu_values, L, estimators) {
  # Crear una tabla para almacenar los resultados de todos los valores de mu
  combined_results <- NULL
  
  # Bucle para iterar sobre cada valor de mu
  for (mu_val in mu_values) {
    # Calcular resultados para el valor actual de mu
    results <- calculate_bias_mse(sample_sizes, R, B, mu_val, L, estimators)
    
    # Agregar una columna con el valor actual de mu en la primera posición
    results <- cbind(mu = mu_val, results)
    
    # Combina los resultados con la tabla principal
    if (is.null(combined_results)) {
      combined_results <- results
    } else {
      combined_results <- rbind(combined_results, results)
    }
  }
  
  # Return the combined table
  return(combined_results)
}



# # Function to generate tables
# generate_table <- function(sample_sizes, R, B, mu_values, L, estimators) {
#   # Lista para almacenar las tablas
#   table_list <- list()
#   
#   # Bucle para iterar sobre cada valor de mu
#   for (mu_val in mu_values) {
#     # Calcular resultados para el valor actual de mu
#     results <- calculate_bias_mse(sample_sizes, R, B, mu_val, L, estimators)
#     #cat("Resultados para mu =", mu_val, "\n")
#     
#     # Imprimir la tabla de resultados
#     table_list[[as.character(mu_val)]] <- results
#   }
#   
#   # Return the list of tables
#   return(table_list)
# }
# 
# # generate_plot <- function(sample_sizes, R, mu_values, L, estimators, ncol = 2, nrow = 2) {
#   # Lista para almacenar los gráficos
#   plot_list <- list()
# 
#   # Bucle para iterar sobre cada valor de mu
#   for (mu_val in mu_values) {
#     # Calcular resultados para el valor actual de mu
#     results <- calculate_bias_mse(sample_sizes, R, mu_val, L, estimators)
#     cat("Resultados para mu =", mu_val, "\n")
# 
#     # Imprimir la tabla de resultados
#     print(results)
# 
#     df <- as.data.frame(results)
# 
#     # Plot Bias y MSE en subgráficos separados
#     plot_bias <- ggplot(df, aes(x = SampleSize, y = Bias, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", linewidth = 0.5) +
#       labs(y = "Bias", x = "Sample size") +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
# 
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_bias = plot_bias + theme(legend.position = "top")
#     }
# 
#     plot_mse <- ggplot(df, aes(x = SampleSize, y = MSE, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", linewidth = 0.5) +
#       labs(y = "MSE", x = "Sample size") +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
# 
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_mse = plot_mse + theme(legend.position = "top")
#     }
# 
#     # Agregar los gráficos a la lista
#     plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
#   }
# 
#   # Organizar los gráficos en una sola figura y añadir leyenda común en la parte inferior
#   combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
#     plot_layout(guides = "collect")
# 
#   # No mostrar la figura aquí, devolver el objeto combined_plot
#   return(combined_plot)
# }
# 



# generate_plot <- function(sample_sizes, R, mu_values, L, estimators, ncol = 2, nrow = 2) {
#   # Lista para almacenar los gráficos
#   plot_list <- list()
#   
#   # Bucle para iterar sobre cada valor de mu
#   for (mu_val in mu_values) {
#     # Calcular resultados para el valor actual de mu
#     results <- calculate_bias_mse(sample_sizes, R, mu_val, L, estimators)
#     cat("Resultados para mu =", mu_val, "\n")
#     
#     # Imprimir la tabla de resultados
#     print(results)
#     
#     df <- as.data.frame(results)
#     
#     # Plot Bias y MSE en subgráficos separados
#     plot_bias <- ggplot(df, aes(x = SampleSize, y = Bias, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", linewidth = 0.5) +
#       labs(y = "Bias", x = "Sample size") +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = bquote(mu == .(mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
#     
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_bias = plot_bias + theme(legend.position = "top")
#     }
#     
#     plot_mse <- ggplot(df, aes(x = SampleSize, y = MSE, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", linewidth = 0.5) +
#       labs(y = "MSE", x = "Sample size") +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = bquote(mu == .(mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
#     
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_mse = plot_mse + theme(legend.position = "top")
#     }
#     
#     # Agregar los gráficos a la lista
#     plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
#   }
#   
#   # Organizar los gráficos en una sola figura y añadir leyenda común en la parte inferior
#   combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
#     plot_layout(guides = "collect")
#   
#   # No mostrar la figura aquí, devolver el objeto combined_plot
#   return(combined_plot)
# }






# generate_plot <- function(sample_sizes, R, mu_values, L, estimators) {
#   # Lista para almacenar los gráficos
#   plot_list <- list()
#   
#   # Bucle para iterar sobre cada valor de mu
#   for (mu_val in mu_values) {
#     # Calcular resultados para el valor actual de mu
#     results <- calculate_bias_mse(sample_sizes, R, mu_val, L, estimators)
#     cat("Resultados para mu =", mu_val, "\n")
#     
#     # Imprimir la tabla de resultados
#     print(results)
#     
#     df <- as.data.frame(results)
#     
#     # Plot Bias y MSE en subgráficos separados
#     plot_bias <- ggplot(df, aes(x = SampleSize, y = Bias, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", size = 0.5) +
#       labs(y = "Bias", x = "Sample size") +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
#     
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_bias = plot_bias + theme(legend.position = "top")
#     }
#     
#     plot_mse <- ggplot(df, aes(x = SampleSize, y = MSE, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", size = 0.5) +
#       labs(y = "MSE", x = "Sample size") +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
#     
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_mse = plot_mse + theme(legend.position = "top")
#     }
#     
#     # Agregar los gráficos a la lista
#     plot_list[[as.character(mu_val)]] <- list(plot_bias, plot_mse)
#   }
#   
#   # Devolver la lista de gráficos
#   return(plot_list)
# }

# generate_plot <- function(sample_sizes, R, mu_values, L, estimators) {
#   # Lista para almacenar los gráficos
#   plot_list <- list()
# 
#   # Bucle para iterar sobre cada valor de mu
#   for (mu_val in mu_values) {
#     # Calcular resultados para el valor actual de mu
#     results <- calculate_bias_mse(sample_sizes, R, mu_val, L, estimators)
#     cat("Resultados para mu =", mu_val, "\n")
# 
#     # Imprimir la tabla de resultados
#     print(results)
# 
#     df <- as.data.frame(results)
# 
#     # Plot Bias y MSE en subgráficos separados
#     plot_bias <- ggplot(df, aes(x = SampleSize, y = Bias, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", size = 0.5) +
#       labs(y = "Bias", x = "Sample size") +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
# 
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_bias = plot_bias + theme(legend.position = "top")
#     }
# 
#     plot_mse <- ggplot(df, aes(x = SampleSize, y = MSE, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", size = 0.5) +
#       labs(y = "MSE", x = "Sample size") +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
# 
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_mse = plot_mse + theme(legend.position = "top")
#     }
# 
#     # Agregar los gráficos a la lista
#     plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
#   }
# 
#   # Organizar los gráficos en una sola figura y añadir leyenda común en la parte inferior
#   combined_plot <- wrap_plots(plot_list, ncol = 2, nrow = 2) +
#     plot_layout(guides = "collect")
# 
#   # Mostrar la figura
#   print(combined_plot)
# }

# ...

# ...

# generate_plots <- function(results, mu_val) {
#   df <- as.data.frame(results)
#   
#   plot_bias <- ggplot(df, aes(x = SampleSize, y = Bias, color = Estimator)) +
#     geom_hline(yintercept = 0) +
#     geom_point(size = 2) +
#     geom_line(linetype = "solid", size = 0.5) +
#     labs(y = "Bias", x = "Sample size") +
#     guides(color = guide_legend(title = "Estimator")) +
#     annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
#   
#   if (mu_val != mu_values[1]) {
#     plot_bias = plot_bias + theme(legend.position = "top")
#   }
#   
#   plot_mse <- ggplot(df, aes(x = SampleSize, y = MSE, color = Estimator)) +
#     geom_hline(yintercept = 0) +
#     geom_point(size = 2) +
#     geom_line(linetype = "solid", size = 0.5) +
#     labs(y = "MSE", x = "Sample size") +
#     guides(color = guide_legend(title = "Estimator")) +
#     annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
#   
#   if (mu_val != mu_values[1]) {
#     plot_mse = plot_mse + theme(legend.position = "top")
#   }
#   
#   return(plot_bias + plot_mse)
# }
# 
# 
# calculate_bias_mse <- function(sample_sizes, R, mu, L, estimators) {
#   true_entropy <- entropy_gamma_sar(L, mu)
#   
#   output <- data.frame(SampleSize = integer(0), Estimator = character(0), Bias = numeric(0), MSE = numeric(0))
#   
#   for (ssize in sample_sizes) {
#     # Generate samples outside the loop
#     samples <- generate_samples(ssize, R, mu, L)
#     
#     for (estimator_name in names(estimators)) {
#       estimator <- estimators[[estimator_name]]
#       v.entropy <- numeric(R)
#       
#       for (r in 1:R) {
#         sample <- samples[[r]]
#         
#         if (grepl("Bootstrap", estimator_name)) {
#           v.entropy[r] <- estimator(sample, B = 10)
#         } else {
#           v.entropy[r] <- estimator(sample)
#         }
#       }
#       
#       mse <- mean((v.entropy - true_entropy)^2)
#       bias <- mean(v.entropy) - true_entropy
#       
#       output <- rbind(output, data.frame(SampleSize = ssize, Estimator = estimator_name, Bias = round(bias, 5), MSE = round(mse, 5)))
#     }
#   }
#   
#   return(output)
# }

