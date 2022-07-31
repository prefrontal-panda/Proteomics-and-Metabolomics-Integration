#This is the workflow for MOFA - based on the examples provided by the developers.

#We start with formating the data for MOFA analysis.
#Setting working directory
setwd("C:/Users/dcho0009/Desktop/multiomics_trial")
#Libraries
library(tidyr)
library(dplyr)

#Read in files
prot <- read.csv('prot_SCx_log_scale_OutReplaceRow.csv', stringsAsFactors = F)
colnames(prot)[1] <- "strain"

met <- read.csv('met_SCx_log_OutReplaceRow.csv', stringsAsFactors = F, fileEncoding = 'UTF-8-BOM')
colnames(met)[1] <- "strain"

#Editing the proteomic dataset
prot_long <- prot %>%
  pivot_longer(!strain, #all columns except strain
               names_to = "features", #proteins to go into column; 'features'
               values_to =  "value") %>% #values to go into coumn; 'value'
  arrange(features) %>% #arranging by proteins (ascending order) 
  mutate(omic = "Proteomic", .after = features) #add 'view' column (or in this case, omic)
#Editing metabolomics dataset
met_long <- met %>%
  pivot_longer(!strain,
               names_to = "features",
               values_to = "value") %>%
  arrange(features) %>%
  mutate(omic = "Metabolomic", .after = features)

#Concatenating both datasets
join <- rbind(prot_long, met_long)

#Saving file
write.csv(join, 'MOFA_input.csv')

#Now we run the MOFA analysis
#library
library(MOFA2)
library(data.table)
library(ggplot2)
library(GGally) #for scatter matrix
library(ggpubr)

#Training model
#Loading in data file
input <- read.csv('MOFA_input.csv', stringsAsFactors = F)
#Checking structure
dim(input)
input[1:5,]
input <- input[,-1]
colnames(input) <- c("sample", "feature", "view", "value") #columns have to be named this
#otherwise: Error in create_mofa_from_df(input) : 
#all(c("sample", "feature", "value") %in% colnames(df)) is not TRUE

#Creating MOFA object
model_obj <- create_mofa_from_df(input)
print(model_obj)

#Visualising data structure
plot_data_overview(model_obj)

#Selecting options
data_opts <- get_default_data_options(model_obj)
head(data_opts)
data_opts$scale_views <- TRUE
head(data_opts)
#If something requires changing
#data_opts$
#head(data_opts)

model_opts <- get_default_model_options(model_obj)
head(model_opts)
#model_opts$likelihoods <- c('bernoulli', 'bernoulli')
#model_opts[["likelihoods"]][["Metabolomic"]] <- 'bernoulli' #data not binary - canot use this.
#model_opts[["likelihoods"]][["Proteomic"]] <- 'bernoulli'
model_opts$num_factors <- 3
head(model_opts)

train_opts <- get_default_training_options(model_obj)
head(train_opts)
#Need to do a few runs to find the model with the highest ELBO
train_opts$maxiter <- 3000
train_opts$convergence_mode <- "slow"
#head(train_opts)

#Training model
model_train <- prepare_mofa(
  object = model_obj,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts)

file = file.path(getwd(), "MOFA_model_trained.hdf5")
trained_model <- run_mofa(model_train, file, save_data = T, use_basilisk = T)

#Loading trained model
filepath <- file.path(getwd(), "MOFA_model_trained.hdf5")
print(filepath)
trained_model <- load_model(filepath)

plot_factor_cor(trained_model)

#Downstream analysis
#Overview of data
plot_data_overview(trained_model)

#Loading metadata
#The metadata is stored as a data.frame object in model@samples_metadata, and it requires at 
#least the column sample. The column group is required only if you are doing multi-group inference.
#The number of rows must match the total number of samples in the model (sum(model@dimensions$N)).
metadata <- read.csv("metadata.csv")
head(metadata)

samples_metadata(trained_model) <- metadata
head(trained_model@samples_metadata)

#Variance decomposition
head(trained_model@cache$variance_explained$r2_total) #total variance
head(trained_model@cache$variance_explained$r2_per_factor) #per factor variance
#plotting variance
plot_variance_explained(trained_model, x="view", y="factor") #individual factor plot
plot_variance_explained(trained_model, x="group", y="factor", plot_total = T)

#Visualisation of single factors
plot_factor(trained_model,
            factors = 1:3,
            color_by = "strain",
            legend = T)
#color_by = "sample",
#shape_by = 'strain',

#Visualisation of combination of factors
plot_factors(trained_model,
             factors = 1:3,
             color_by = "strain")#,
# shape_by = "sample",
# alpha = 0.8)

#Visualisation of feature weights
plot_weights(trained_model,
             view = "Metabolomic",
             factors = 1,
             nfeatures = 10,
             text_size = 4)

plot_top_weights(trained_model,
                 view = "Metabolomic",
                 factors = 1,
                 nfeatures = 10)

#Visualisation of patterns
#Heatmap
plot_data_heatmap(trained_model,
                  factor = 1,
                  features = 25,
                  view = "Proteomic",
                  culster_rows = T)

plot_data_heatmap(trained_model,
                  factor = 1,
                  features = 25,
                  view = "Metabolomic",
                  cluster_rows = T)

#Scatterplot
plot_data_scatter(trained_model,
                  factor = 1,
                  features = 6,
                  view = "Metabolomic",
                  color_by = "strain",
                  add_lm = T)