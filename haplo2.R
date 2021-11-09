library('blupADC', lib.loc = "/usr/home/qgg/egill/R/x86_64-pc-linux-gnu-library/4.0")

data_path = "/usr/home/qgg/egill/Project4/MoBPS/20Test_05_11_2021/H1/"
data_name = "haplo"
data_type = "Haplotype"

kinship_result=cal_kinship(
                input_data_path=data_path,       # input data path 
				input_data_name=data_name,       # input data name,  for vcf data, you don't need to add the suffix  
				input_data_type=data_type,           # input data type
				phased_genotype=TRUE,                 #whether the vcf data has been phased
				haplotype_window_nSNP=5,         # options 
                kinship_type=c("G_A","G_Ainv"),           #type of  kinship matrix
				bigmemory_cal=TRUE,             # format conversion via bigmemory object
				bigmemory_data_path=getwd(),    # path of bigmemory data 
				bigmemory_data_name="test_blupADC", #name of bigmemory data 				
                return_result=TRUE)              #return result              
G_Amatrix=kinship_result$G_A$A
G_Ainvmatrix=kinship_result$G_A$Ainv

write.table(G_Amatrix,"testA",row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = "\t")
write.table(G_Ainvmatrix,"testA",row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = "\t")

