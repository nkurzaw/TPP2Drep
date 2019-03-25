#' @importFrom stats lm
#' @importFrom stats coef
getStartParameters <- function(df, unique_temp, len_temp,
                               slopEC50 = FALSE){
  # function to generate vector of start parameters
  # for h1 model optimization
  log_conc <- log2_value <- temperature <- NULL
  
  start_par <- c(mean(unique(df$log_conc)[
    which(is.finite(unique(df$log_conc)))]))
    if(slopEC50){
      start_par <- c(start_par, 0)
    }
    start_par <- c(
      start_par,
      coef(lm(log2_value ~ log_conc, 
              filter(df, is.finite(log_conc))))[[2]],
      max(unlist(lapply(unique_temp, function(temp){
        max(filter(df, temperature == temp)$log2_value) -
          min(filter(df, temperature == temp)$log2_value)
        }))),
      unlist(lapply(unique_temp, function(temp){
        mean(filter(df, temperature == temp)$log2_value)
        })),
      rep(0, len_temp))
    return(start_par)
}

#' @importFrom stats lm
#' @importFrom stats coef
getStartParametersRep <- function(
    df, unique_temp, len_temp, slopEC50 = FALSE){
    # function to generate vector of start parameters
    # for h1 model optimization
    log_conc <- log2_value <- temperature <- NULL
    
    start_par <- c(mean(unique(df$log_conc)[
        which(is.finite(unique(df$log_conc)))]))
    if(slopEC50){
        start_par <- c(start_par, 0)
    }
    start_par <- c(
        start_par,
        coef(lm(log2_value ~ log_conc, 
                filter(df, is.finite(log_conc))))[[2]],
        max(unlist(lapply(unique_temp, function(temp){
            max(filter(df, temperature == temp)$log2_value) -
                min(filter(df, temperature == temp)$log2_value)
        }))),
        (df %>% group_by(replicate, temperature) %>% 
             summarize(mlv = mean(log2_value)))$mlv,
        rep(0, len_temp))
    return(start_par)
}

paste_rmNA <- function(x, sep = "|"){
  # function that pastes non-dedundant and na-filtered
  # vector elements into a string
  return(paste(x[!is.na(x)], collapse = sep))
}

getOptLimits <- function(ec50Limits, len_temp, 
                         len_rep_temp = NULL,
                         slopEC50 = FALSE){
  # function to generate list of vectors defining lower 
  # and upper limits of optimization parameters
  if(!slopEC50){
      if(!is.null(len_rep_temp)){
          out_list <- list(lower = c(ec50Limits[1],
                                     rep(-Inf, 2 + len_rep_temp), 
                                     rep(0, len_temp)),
                           upper = c(ec50Limits[2],
                                     rep(Inf, 2 + len_rep_temp), 
                                     rep(1, len_temp)))
      }else{
          out_list <- list(lower = c(ec50Limits[1],
                                     rep(-Inf, 2 + len_temp), 
                                     rep(0, len_temp)),
                           upper = c(ec50Limits[2],
                                     rep(Inf, 2 + len_temp), 
                                     rep(1, len_temp)))
      }
    
  }else{
      if(!is.null(len_rep_temp)){
          out_list <- list(lower = c(ec50Limits[1],
                                     rep(-Inf, 3 + len_rep_temp), 
                                     rep(0, len_temp)),
                           upper = c(ec50Limits[2],
                                     rep(Inf, 3 + len_rep_temp), 
                                     rep(1, len_temp)))
      }else{
          out_list <- list(lower = c(ec50Limits[1],
                                     rep(-Inf, 3 + len_temp), 
                                     rep(0, len_temp)),
                           upper = c(ec50Limits[2],
                                     rep(Inf, 3 + len_temp), 
                                     rep(1, len_temp)))
      }
  }
  return(out_list)
}

checkDfColumns <- function(df){
  # internal function to check all needed columns
  # are presented in supplied input data frame
  req_coln <- c("representative",
                "clustername",
                "temperature",
                "log_conc",
                "log2_value")
  coln_df <- colnames(df)
  if(!all(req_coln %in% coln_df)){
    stop(paste(c("Input data frame, requires at least the", 
                 "following columns: representative,", 
                 "clustername temperature, log_conc,",
                 "log2_value!"), collapse = " "))
  }
}

dense_rank_2var <- function(x, y, all_possible = TRUE){
    if(all_possible){
        exp_grid <- expand.grid(unique(x),
                                unique(y)) %>% 
            tbl_df %>% 
            dplyr::select(x = Var1, y = Var2) %>% 
            arrange(x, y) %>% 
            mutate(z = 1:n())
        
        df <- data_frame(x, y) %>% 
            left_join(exp_grid, by = c("x", "y"))
        
        return(df$z)
    }else{
        ndup_df <- data_frame(x, y) %>% 
            filter(!duplicated(paste0(as.character(x),
                                      as.character(y)))) %>% 
            arrange(x, y) %>% 
            mutate(z = 1:n())
        
        df <- data_frame(x, y) %>% 
            left_join(ndup_df, by = c("x", "y"))
        
        return(df$z)
    }
}