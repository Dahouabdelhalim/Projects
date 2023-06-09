############################## Config Settings ##############################

#General settings
options(stringsAsFactors = FALSE)

#Indicate the directory with original data
original_directory = "/Original_data"

#Indicate the directory for input files
input_directory = "/Preprocessed_data"

#Indicate the output directory
output_directory = "/Generated_results"

#List excel files with original data
HC_excel_files = list("HC_combinedfiles_network_2days.xlsx", 
                      "HC_combinedfiles_network_10days.xlsx", 
                      "HC_combinedfiles_network_8weeks.xlsx")
PHC_excel_files = list("PHC_combinedfiles_network_2days.xlsx", 
                       "PHC_combinedfiles_network_10days.xlsx", 
                       "PHC_combinedfiles_network_8weeks.xlsx")

time_label = c("2 days", "10 days", "8 weeks")
