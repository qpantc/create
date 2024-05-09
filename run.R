# trait_ = Input1
# Input_path = Input2
# Output_path = Input3


library(phytools)
library(tidyverse)
climate_list <- c('MAT','MAP','PAW')

## FILES in HPC 
MATRIX_PATH <- '/scratch/gent/vo/000/gvo00074/vsc44253/phd/paper2/HPC/Distance.matrix/sma.distance/Outputs/'
trait_combinations <- read.csv('/scratch/gent/vo/000/gvo00074/vsc44253/phd/paper2/HPC/Distance.matrix/sma.distance/Outputs/trait.pairs.csv')


## FILES upload
response_distance_path <- paste0(Input_path,'response_distance/')
Trait_response <- read.csv(paste0(Input_path,'Trait.environment.response.csv'))
species_df <- read.csv(paste0(Input_path,'analyzed_species.csv'))
Trait_response <- left_join(Trait_response,species_df[c('accepted_bin','cluster')])


for (climate_ in climate_list) {
    all_results <- data.frame()
    for (i in c(1:nrow(trait_combinations))) {
        trait_combination_ <- trait_combinations$combination[i]
        load(paste0(MATRIX_PATH,trait_combination_,'.RData')) # file save as distance_matrix
        ### lambda
        one_test <- subset(Trait_response, trait == trait_ & climate == climate_)
        distance_matrix[is.infinite(distance_matrix)]=0
        distance_matrix <- round(distance_matrix,3)
        tree_matrix <- distance_matrix

        phy.tree = as.phylo(hclust(as.dist(tree_matrix), method = "ward.D"))
        
        Pagels_lambda <- list(lambda = NA, P = NA)
        # Blomberg_K <- list(K = NA, P = NA)
        # try(Blomberg_K <- phylosig(phy.tree,setNames(one_test$Estimate,one_test$cluster),test=TRUE,method="K"))
        try(Pagels_lambda <- phylosig(phy.tree,setNames(one_test$Estimate,one_test$cluster),test=TRUE,method = "lambda"))
        
        ## cor.test
        distance_matrix[lower.tri(distance_matrix)] <- NA
        
        cor_distance <- distance_matrix
        
        load(paste0(response_distance_path,climate_,'_',trait_,'.RData')) # file save as response_distance_matrix
        response_distance_matrix[lower.tri(response_distance_matrix)] <- NA
        regression_distance <- response_distance_matrix
        diag(regression_distance) <- NA
        regression_distance_list <- regression_distance[!is.na(regression_distance)]
        # regression_distance_list <- (regression_distance_list-min(regression_distance_list))/max(regression_distance_list)
        cor_distance_list <- cor_distance[!is.na(regression_distance)]
        
        cor_test_summary <- cor.test(regression_distance_list,cor_distance_list,method = 'spearman')
        
        one_result <- data.frame(
        trait = trait_,
        climate = climate_,
        trait_combination = trait_combination_,
        cor = cor_test_summary$estimate,
        p_value = cor_test_summary$p.value,
        # Blomberg_K = Blomberg_K$K,
        # Blomberg_K_p = Blomberg_K$P,
        Pagels_lambda = Pagels_lambda$lambda,
        Pagels_lambda_p = Pagels_lambda$P#,lambda_phy_p1 = mod_phy_lambda[3]
        )
        if (one_result$Pagels_lambda_p < 0.05 & one_result$Pagels_lambda<1) {
            all_results <- rbind(all_results,one_result)
        }
        if (nrow(all_results)>1000) {
            all_results <- all_results[order(all_results$Pagels_lambda,decreasing = TRUE),]
            all_results <- all_results[1:1000,]
        }
        if (i %% 1000 == 0) {
            write.csv(all_results,
                        paste0(Output_path,climate_,'_',trait_,'_regression.similarity.cor.csv'),
                        row.names = FALSE)
            print(paste(i,Sys.time(),trait_combination_))
        }
    }
}
