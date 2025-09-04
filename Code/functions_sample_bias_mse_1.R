
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



generate_plot_gi0_esp <- function(results_gi0, mu_values, selected_estimators, ncol = 2, nrow = 2) {
  
  plot_list <- list()
  
  for (mu_val in mu_values) {
    
    df <- results_gi0[[as.character(mu_val)]]
    
    df_filtered <- df[df$Estimator %in% names(selected_estimators), ]
    df_filtered$Estimator <- selected_estimators[df_filtered$Estimator]
    df_filtered$Estimator <- as.character(df_filtered$Estimator)
    
    plot_bias <- ggplot(df_filtered, aes(x = n, y = Bias, color = Estimator)) +
      geom_hline(yintercept = 0) +
      geom_point(size = 2.4) +
      geom_line(linetype = "solid", linewidth = 1.2, alpha=.7) +
      labs(y = "Bias", x = expression(italic(n))) +
      scale_x_continuous(breaks = c(9, 25, 49, 81, 121)) +  #  escala personalizada
      scale_color_manual(values = pal_jama()(7)[1:6], labels = TeX(df_filtered$Estimator)) + # Eliminar labels y TeX
      annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 4.0, vjust = -0.1, size = 3) +
      coord_cartesian(clip = 'off') +
      theme(legend.position = "bottom", legend.box = "horizontal", legend.box.margin = margin(0, 0, 0, 0)) +
      guides(color = guide_legend(nrow = 1))
    
    plot_mse <- ggplot(df_filtered, aes(x = n, y = MSE, color = Estimator)) +
      geom_hline(yintercept = 0) +
      geom_point(size = 2.4) +
      geom_line(linetype = "solid", linewidth = 1.2, alpha=.7) +
      labs(y = "MSE", x = expression(italic(n))) +
      scale_x_continuous(breaks = c(9, 25, 49, 81, 121)) +  # 
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





generate_plot_gi0 <- function(sample_sizes, R, B, mu_values, alpha, L, estimators, ncol = 2, nrow = 2) {
  
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
    # Eliminar la leyenda para todos los grC!ficos excepto el primero de cada fila
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
    
    
    # Eliminar la leyenda para todos los grC!ficos excepto el primero de cada fila
    if (mu_val != mu_values[1]) {
      plot_mse = plot_mse + theme(legend.position = "top")
    }
    
    # Agregar los grC!ficos a la lista
    plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
  }
  
  # Organizar los grC!ficos en una sola figura y aC1adir leyenda comC:n en la parte inferior
  combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
    plot_layout(guides = "collect")
  
  # No mostrar la figura aquC-, devolver el objeto combined_plot
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


calculate_bias_mse_r <- function(sample_sizes, R, B, mu, L, alpha_values, estimators) {
  output <- data.frame(n = integer(0), Estimator = character(0), Alpha = numeric(0), Bias = numeric(0), MSE = numeric(0))
  
  for (alpha in alpha_values) {
    true_entropy <- entropy_renyi_gamma_sar(L, mu, alpha)
    
    for (ssize in sample_sizes) {
      samples <- generate_samples(ssize, R, mu, L)
      
      for (estimator_name in names(estimators)) {
        estimator <- estimators[[estimator_name]]
        v.entropy <- numeric(R)
        
        for (r in 1:R) {
          sample <- samples[[r]]
          
          if (grepl("Bootstrap", estimator_name)) {
            v.entropy[r] <- estimator(sample, B = B, a = alpha)
          } else {
            v.entropy[r] <- estimator(sample, a = alpha)
          }
        }
        
        mse <- mean((v.entropy - true_entropy)^2)
        bias <- mean(v.entropy) - true_entropy
        
        output <- rbind(output, data.frame(n = ssize, Estimator = estimator_name, Alpha = alpha, Bias = round(bias, 5), MSE = round(mse, 5)))
      }
    }
  }
  
  return(output)
}


generate_plot_renyi <- function(results_df, alpha_values, selected_estimators, ncol = 1, nrow = 6) {
 # library(ggplot2)
 # library(patchwork)
  
  plot_list <- list()
  
  # for (alpha_val in alpha_values) {
  #   df <- results_df[results_df$Alpha == alpha_val, ]
  #   df_filtered <- df[df$Estimator %in% names(selected_estimators), ]
  #   df_filtered$Estimator <- selected_estimators[df_filtered$Estimator]
  #   df_filtered$Estimator <- as.character(df_filtered$Estimator)
  #   
  #   plot_bias <- ggplot(df_filtered, aes(x = n, y = Bias, color = Estimator)) +
  #     geom_hline(yintercept = 0, color = "gray50", linetype = "dashed",linewidth = 0.5) +
  #     geom_point(size = 2.4) +
  #     geom_line(linetype = "solid", linewidth = 1.2, alpha = 0.7) +
  #     labs(y = "Bias", x = expression(italic(n))) +
  #     scale_x_continuous(breaks = unique(df_filtered$n)) +
  #     scale_color_manual(values = c("blue", "red"), labels = TeX(df_filtered$Estimator)) +#TeX(df_filtered$Estimator)
  #     ggtitle(bquote(lambda == .(alpha_val))) +
  #     theme_minimal() +
  #     theme(text = element_text(family = "serif"),
  #           axis.text = element_text(size = 16),
  #           axis.title = element_text(size = 18),
  #           plot.title = element_text(size = 20, hjust = 0),  # Aumenta el tamaño del título de lambda
  #           legend.position = "bottom",
  #           legend.box = "horizontal",
  #           legend.text = element_text(size = 16),  # Aumenta el tamaño del texto de la leyenda
  #           legend.title = element_text(size = 18)) +  # Aumenta el tamaño del título de la leyenda
  #     guides(color = guide_legend(nrow = 1))
  #   
  #   plot_mse <- ggplot(df_filtered, aes(x = n, y = MSE, color = Estimator)) +
  #     geom_hline(yintercept = 0,color = "gray50", linetype = "dashed",linewidth = 0.5) +
  #     geom_point(size = 2.4) +
  #     geom_line(linetype = "solid", linewidth = 1.2, alpha = 0.7) +
  #     labs(y = "MSE", x = expression(italic(n))) +
  #     scale_x_continuous(breaks = unique(df_filtered$n)) +
  #     scale_color_manual(values = c("blue", "red"), labels = TeX(df_filtered$Estimator)) +
  #     theme_minimal() +
  #     theme(text = element_text(family = "serif"),
  #           axis.text = element_text(size = 16),
  #           axis.title = element_text(size = 18),
  #           plot.title = element_text(size = 20, hjust = 0),  # Aumenta el tamaño del título de lambda
  #           legend.position = "bottom",
  #           legend.box = "horizontal",
  #           legend.text = element_text(size = 16),  # Aumenta el tamaño del texto de la leyenda
  #           legend.title = element_text(size = 18)) +  # Aumenta el tamaño del título de la leyenda
  #     guides(color = guide_legend(nrow = 1))
  #   
  #   # Combine bias and MSE plots into one
  #   combined_plot <- plot_bias / plot_mse
  #   
  #   plot_list[[as.character(alpha_val)]] <- combined_plot
  # }
  # 
  # # Arrange all plots
  # combined_plot_all <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
  #   plot_layout(guides = "collect")
  
  library(ggplot2)
  library(patchwork)
  
  for (alpha_val in alpha_values) {
    df <- results_df[results_df$Alpha == alpha_val, ]
    df_filtered <- df[df$Estimator %in% names(selected_estimators), ]
    df_filtered$Estimator <- selected_estimators[df_filtered$Estimator]
    df_filtered$Estimator <- as.character(df_filtered$Estimator)
    
    plot_bias <- ggplot(df_filtered, aes(x = n, y = Bias, color = Estimator)) +
      geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.5) +
      geom_point(size = 3.2) +  # Aumenta el tamaño de los puntos
      geom_line(linetype = "solid", linewidth = 1.5, alpha = 0.8) +  # Aumenta el grosor de la línea
      labs(y = "Bias", x = expression(italic(n))) +
      #scale_x_continuous(breaks = unique(df_filtered$n)) +
      scale_x_continuous(breaks = unique(df_filtered$n), minor_breaks = NULL)+
      scale_y_continuous(minor_breaks = NULL)+
      scale_color_manual(values = c("blue", "red"), labels = TeX(df_filtered$Estimator)) +
      ggtitle(bquote(lambda == .(alpha_val))) +
      theme_minimal() +
      theme(text = element_text(family = "serif"),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            plot.title = element_text(size = 20,  hjust = 0),  # Aumenta el tamaño del título de lambda
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.text = element_text(size = 16),  # Aumenta el tamaño del texto de la leyenda
            legend.title = element_text(size = 18)) +  # Aumenta el tamaño del título de la leyenda
      guides(color = guide_legend(nrow = 1))
    
    plot_mse <- ggplot(df_filtered, aes(x = n, y = MSE, color = Estimator)) +
      geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.5) +
      geom_point(size = 3.2) +
      geom_line(linetype = "solid", linewidth = 1.5, alpha = 0.8) +
      labs(y = "MSE", x = expression(italic(n))) +
      #scale_x_continuous(breaks = unique(df_filtered$n)) +
      scale_x_continuous(breaks = unique(df_filtered$n), minor_breaks = NULL)+
      scale_y_continuous(minor_breaks = NULL)+
      scale_color_manual(values = c("blue", "red"), labels = TeX(df_filtered$Estimator)) +
      theme_minimal() +
      theme(text = element_text(family = "serif"),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18)) +
      guides(color = guide_legend(nrow = 1))
    
    # Combinar Bias y MSE en formato vertical
    combined_plot <- plot_bias / plot_mse
    
    plot_list[[as.character(alpha_val)]] <- combined_plot
  }
  
  # Organizar todos los gráficos
  combined_plot_all <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
    plot_layout(guides = "collect")
  
  
  return(combined_plot_all)
}

# calculate_bias_mse <- function(sample_sizes, R, B, mu, L, estimators) {
#   true_entropy <- entropy_gamma_sar(L, mu)
#   
#   output <- data.frame(n = integer(0), Estimator = character(0), Bias = numeric(0), MSE = numeric(0))
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
#           v.entropy[r] <- estimator(sample, B = B)
#         } else {
#           v.entropy[r] <- estimator(sample)
#         }
#       }
#       
#       mse <- mean((v.entropy - true_entropy)^2)
#       bias <- mean(v.entropy) - true_entropy
#       
#       output <- rbind(output, data.frame(n = ssize, Estimator = estimator_name, Bias = round(bias, 5), MSE = round(mse, 5)))
#     }
#   }
#   
#   return(output)
# }
# 


# 
generate_plot <- function(sample_sizes, R, B, mu_values, L, estimators, ncol = 2, nrow = 2) {
  # Lista para almacenar los grC!ficos
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
    # Eliminar la leyenda para todos los grC!ficos excepto el primero de cada fila
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


    # Eliminar la leyenda para todos los grC!ficos excepto el primero de cada fila
    if (mu_val != mu_values[1]) {
      plot_mse = plot_mse + theme(legend.position = "top")
    }

    # Agregar los grC!ficos a la lista
    plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
  }

  # Organizar los grC!ficos en una sola figura y aC1adir leyenda comC:n en la parte inferior
  combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
    plot_layout(guides = "collect")

  # No mostrar la figura aquC-, devolver el objeto combined_plot
  return(combined_plot)
}



# 
# generate_plot <- function(sample_sizes, R, B, mu_values, L, estimators, ncol = 2, nrow = 2) {
#   # Lista para almacenar los grC!ficos
#   plot_list <- list()
# 
#   for (mu_val in mu_values) {
#     # Calcular resultados para el valor actual de mu
#     results <- calculate_bias_mse(sample_sizes, R, B, mu_val, L, estimators)
#     df <- as.data.frame(results)
# 
#     # Dividir estimadores en dos grupos (sC3lidos y punteados)
#     solid_estimators <- names(estimators)[1:(length(estimators) / 2)]
#     dashed_estimators <- names(estimators)[-c(1:(length(estimators) / 2))]
# 
#     # Filtrar dataframes para cada tipo de lC-nea
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
    
    # Agregar una columna con el valor actual de mu en la primera posiciC3n
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



