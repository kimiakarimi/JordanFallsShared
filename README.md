# JordanFallsShare
The files here include data sets and codes for the following dissertaion chapter: Assessing Controls on Nutrient Loading at the Watershed Scale through Data-Driven Modeling.

Items in this repository:

data_base.rds- This is a R data file that includes all the input data sets for the base model. It needs to be imported into R with the appropriate command.

data_buffer.rds- This is a R data file that includes all the input data sets for the buffer model. It needs to be imported into R with the appropriate command.

data_scm.rds- This is a R data file that includes all the input data sets for the SCM model. It needs to be imported into R with the appropriate command.

data_buffer_scm.rds- This is a R data file that includes all the input data sets for the Buffer+SCM model. It needs to be imported into R with the appropriate command.

RStan_codes.R- RStan code needed to run the models. The "rstan" (https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) and "rstudioapi" packages need to be installed prior to running the models. 
# Example
Step 1: Install "rstan" and "rstudioapi" packages. 

Step 2: Compile the "stanmodelcode_TN_base" model (line #5 in the "RStan_codes.R" code) .

Step 3: Read the base model input data set (data_base.rds), if necessary (line #736)

Step 4: Run the stan function (line #738). Suggest using the following parameters: iter=4000, warmup=2000, thin=5, chains=3,cores=3,adapt_delta =0.99 ,max_treedepth =25. The estimated runtime on desktop is about one day. Note that more iterations and warmup steps may improve model convergence.
