# GPCR - do BLAST on protein sequences:
#
# description:  GPCR data 
#               Comparing PrEST amino acid sequence for off-target binding
# required_packages: tidyverse, reshape2, scales
# creator: Annika Bendes
######################################################################
# empty global environment
rm(list = ls(all.names = TRUE))

# load packages
library(tidyverse) #Innehåller dplyr + ggplot2
library(reshape2)
library(scales)
library(gridExtra)
require(ggrepel)
library(RColorBrewer)
require(dplyr)
library(readxl)
library(Biostrings)
library("xlsx")
########  Install rBLAST
#install.packages("devtools")
#devtools::install_github("mhahsler/rBLAST")
library(rBLAST)
library(taxonomizr)
library(ggplot2)


## Load antibody tables:
reactivity_information <- read_excel("antibody_result_table.xlsx")

# Load antibody epitope sequences
PrEST_information <- read_excel("antibody_epitope_sequence_information.xlsx")

# Load protein sequences
protein_sequences <- read_excel("protein_sequences.xlsx")


#### Take out all the cross-reactive antibodies:

crossreacrive_abs <- reactivity_information[grep("TRUE", reactivity_information$crossreact),  ]


####      Loop through every ab to compare the PrESTs aa:

selected_binders <- crossreacrive_abs$unique_binder_name

crossreactivity_information_table <- NULL

for(i in selected_binders){
  
  ab_position <- grep(i, crossreacrive_abs$unique_binder_name)
  
  crossreactive_vector_i <- unlist(strsplit(crossreacrive_abs$crossreact_names_unique[grep(i, crossreacrive_abs$unique_binder_name)], ","))
  
  crossreactivity_information_i  <- cbind(crossreacrive_abs[     ab_position, 1:9], crossreactive_vector_i )
  
  crossreactivity_information_table <- rbind(crossreactivity_information_table, crossreactivity_information_i)
  
}


##################        Remove all antibodies that are not HPA antibodies (because they don't have PrESTs)
crossreactivity_information_table_all_abs <- crossreactivity_information_table
crossreactivity_information_table <- crossreactivity_information_table[grep("HPA", crossreactivity_information_table$antibody),      ]

######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##
##                                                  Boxplot where we compare sequence length for the different 
##
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####                                                Is the risk of off-target binding dependent or independent on PrEST/Protein length?

#                                 Take out number of amino acids for the PrEST for on target (Group 1) and off-target (Group 2)
all_hpa_abs_table <- reactivity_information[grep("HPA", reactivity_information$antibody), ]

amino_acid_length_table <-NULL

for (ab in 1:nrow(all_hpa_abs_table)) {
  antibody <- all_hpa_abs_table$unique_binder_name[ab]
  hpa_id <- all_hpa_abs_table$antibody[ab]
  target <- all_hpa_abs_table$gene_name[ab]
  PREST_name <- all_hpa_abs_table$prest[ab]
  prest_seq <- PrEST_information$`PrEST seq (aa)`[grep(antibody, PrEST_information$unique_binder_name)]
  prest_length <- nchar(prest_seq)
  cross_reactivity <- all_hpa_abs_table$crossreact[ab]
  
  amino_acid_lenght_i <- c(antibody, target, hpa_id, cross_reactivity, PREST_name, prest_seq, prest_length)
  amino_acid_length_table <- rbind(amino_acid_length_table, amino_acid_lenght_i)
  
  
}

amino_acid_length_table <- data.frame(amino_acid_length_table)
rownames(amino_acid_length_table) <- amino_acid_length_table[,1]
colnames(amino_acid_length_table) <- c("unique_binder_name","gene_name", "antibody", "crossreact",  "prest", "prest_seq", "prest_length" )
amino_acid_length_table$prest_length <- as.numeric(as.character(amino_acid_length_table$prest_length))

amino_acid_length_table_wo_ctrl <- amino_acid_length_table[-grep("ctrl", amino_acid_length_table$unique_binder_name),]

### Add only one of the HPA008070 antibodies (positive control ab)
amino_acid_length_table_with_one_ctrl <- rbind(amino_acid_length_table_wo_ctrl, amino_acid_length_table[219,]) # Same ab, select one because they have same PrEST length

########## Add now information if the ab capture anything at all

amino_acid_length_table_with_one_ctrl$ab.status <- amino_acid_length_table_with_one_ctrl$crossreact

# Take out the antibodies who did not capture anything:
no_target_abs_table <- reactivity_information %>%
  filter(detect_target == "FALSE", crossreact == "FALSE" )

no_target_abs_table$no_target <- "YES"

amino_acid_length_table_with_one_ctrl_no_target_info <- merge(amino_acid_length_table_with_one_ctrl, no_target_abs_table[, c(1, 35)], by = "unique_binder_name", all = TRUE)

### Remove the commercial ones that are not HPA:
amino_acid_length_table_with_one_ctrl_no_target_info <- amino_acid_length_table_with_one_ctrl_no_target_info[-which(is.na(amino_acid_length_table_with_one_ctrl_no_target_info$gene_name)), ]


amino_acid_length_table_with_one_ctrl_no_target_info$ab.status <- as.character(amino_acid_length_table_with_one_ctrl_no_target_info$ab.status )

amino_acid_length_table_with_one_ctrl_no_target_info <- amino_acid_length_table_with_one_ctrl_no_target_info %>% mutate(ab.status = case_when(no_target == "YES" ~ "no_target",
                                                                                      T ~ ab.status))

amino_acid_length_table_with_one_ctrl_no_target_info <- amino_acid_length_table_with_one_ctrl_no_target_info %>% mutate(ab.status = case_when(ab.status == "TRUE" ~ "off-target",
                                                                                                                                              T ~ ab.status))
amino_acid_length_table_with_one_ctrl_no_target_info <- amino_acid_length_table_with_one_ctrl_no_target_info %>% mutate(ab.status = case_when(ab.status == "FALSE" ~ "on-target",
                                                                                                                                              T ~ ab.status))


amino_acid_length_table_with_one_ctrl_no_target_info$ab.status <- as.factor(amino_acid_length_table_with_one_ctrl_no_target_info$ab.status )


###### Boxplots

my_comparisons <- list( c("no_target", "off-target"), c("off-target", "on-target"), c("no_target", "on-target") )

prest_length_boxplot <- ggplot(amino_acid_length_table_with_one_ctrl_no_target_info, aes(x = ab.status,
                                                                    y = as.numeric(prest_length),
                                                                    fill = ab.status
)) +
  theme_light() +
  labs(title = paste0("PrEST length"),
       subtitle = paste0("(Only HPA) Number of antibodies: ", nrow(amino_acid_length_table_with_one_ctrl)),
       tag=format(Sys.time(),"%Y-%m-%d")) +
  geom_boxplot( ) +
  theme(plot.subtitle = element_text(size=7),
        plot.caption = element_text(size=5)) +
  ylab(paste0("Length of antigen [N]")) +
  stat_compare_means(label.y = 190) +  # add global KW comparison
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  scale_fill_brewer(palette="Dark2") +
  geom_beeswarm( cex = 0.7) +
  theme( axis.line = element_line(size = 1) , axis.text = element_text(size = 10) ,
         axis.title.y = element_text(size = 23, face = "bold"), plot.title = element_text(size = 45, face = "bold"), legend.title = element_text(size = 16), 
         axis.ticks.length.x = unit(0.3, "cm"), axis.ticks.length.y = unit(0.3, "cm"), legend.position="none" ) 



#pdf(paste0("PrEST length comparison on off and no target abs_color_", format(Sys.time(),"%Y-%m-%d_%H%M%S"),
#           "Kruskal-wallis.pdf"),
#    width = 7, height = 10, useDingbats=F, onefile=T)
print(prest_length_boxplot)
#dev.off()


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                                         Alluvial plots
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sankey_data_matrix <- crossreactivity_information_table_all_abs

sankey_data_matrix$crossreactive_gpcr <- NA
sankey_data_matrix$crossreactive_vector_i <- gsub(" ", "", sankey_data_matrix$crossreactive_vector_i, fixed = TRUE)

for ( r in (1:nrow(sankey_data_matrix))){
  gpcr_cross <- as.character(sankey_data_matrix$crossreactive_vector_i[r])
  
  sankey_data_matrix$crossreactive_gpcr[r] <- gsub("\\..*","", gpcr_cross)
  
}


############### Prepare a pivot table with binder name and crossreactive_gpcr
#install.packages("pivottabler")
library(pivottabler)

pivot <- qpvt(sankey_data_matrix, "unique_binder_name", "crossreactive_gpcr", "n()") # TOC = Train Operating Company 


#install.packages("reshape2")             # Install & load reshape2
library("reshape2")

pivot <- dcast(sankey_data_matrix,                # Create pivot table
               unique_binder_name ~ crossreactive_gpcr)


data_long <- melt(pivot)

data_long <- data_long[-grep(0, data_long$value) , ]


colnames(data_long) <- c("source", "target", "value")
data_long$target <- paste(data_long$target, " ", sep="")

######### Sankey using alluvial package

library(ggalluvial)

sba_name <- unique(crossreactivity_information_table_all_abs$sba)

#library(dplyr)
q <- reactivity_information %>%
  filter(crossreact_names_unique != "") %>%
  select(gene_name_id, crossreact_names_unique, sba) %>%
  separate_rows(crossreact_names_unique, sep = ", ") %>%
  separate(crossreact_names_unique, into = c("other_gpcr", "other_interactor"), sep = "\\.")




#pdf("alluvial plot crossreacting antibodies vs crossreacting gpcr_221018_per sba_calrcl excluded_lines.pdf", width = 8, height = 8, useDingbats=F)
#layout(matrix(ncol=2, data=1:4, byrow = T))


#par(cex = 0.7, mar = c(0, 0, 0, 0))

for (sba in sba_name){
  
  SBA <- sba
  alluvial_plot <- ggplot(q %>% filter(sba == SBA), aes(axis1 = gene_name_id, axis2 = other_gpcr)) +
    geom_alluvium(aes(fill = gene_name_id), colour = "black", show.legend = F) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_void()+ 
    ggtitle(paste0(SBA))
  
  
  print(alluvial_plot)
  
  
}

#dev.off()

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


##########################                                                    BLAST                                           ###########################



##############          Prepare fasta files, one for the PrEST and one for the target + cross-reactive proteins:


library(msa)
library(seqinr)
library(Biostrings)
library(stringr)


binder_names <- unique(crossreactivity_information_table$unique_binder_name)
cross_reactive_targets_no_aa <- NULL
blast_result_big_table <- NULL

for(ab in binder_names[c(1:length(binder_names))]){
  
  #### Extract first PrEST
  
  
  
  PrEST_Antibody <- PrEST_information[  grep( ab, PrEST_information$unique_binder_name    ), 22]
  
  ### Extract crossreactive target:
  
  # Take out the sequences that are for the cross-reactive antibody
  
  ### Take out also the sequence for the target of interest:
  
  target_of_interest <- unique(crossreactivity_information_table$gene_name[grep( ab , crossreactivity_information_table$unique_binder_name)])
  
  cross_info_table_ab <- crossreactivity_information_table[grep( ab , crossreactivity_information_table$unique_binder_name), ]
  
  ##### Add PrEST to fasta file:
  sequence_list <- c(PrEST_Antibody )
  seq_names <- c("PrEST")
  
  write.fasta(sequences = list(sequence_list), names = seq_names, file.out = "PrEST_fastafile", open = "w", nbchar = 60, as.string = FALSE )
  
  
  ###### Add target of interest to a new fasta file:
  
  sequence_list <- as.character(protein_sequences[grep(paste0("^",target_of_interest, "$"), protein_sequences$`Receptor Name` ), 7])
  seq_names <- as.character(target_of_interest)
  
  write.fasta(sequences = list(sequence_list), names = seq_names, file.out = "proteins_fastafile", open = "w", nbchar = 60, as.string = FALSE )
  
  used_targets <- NULL
  
  
  for(i in 1:nrow(cross_info_table_ab)){
    
    
    crossreactive_target_i <- as.character(cross_info_table_ab$crossreactive_vector_i[i])
    crossreactive_target_i <- strsplit(crossreactive_target_i, ".", fixed = TRUE)
    crossreactive_target_i <- str_trim( unlist(crossreactive_target_i)[1], "left")
    
    aa_protein_construct <- as.character(protein_sequences[grep(    paste0("^",crossreactive_target_i, "$")     , protein_sequences$`Receptor Name` ), 7])
    
    
    if( nchar(aa_protein_construct) > 12  && crossreactive_target_i %in% used_targets == FALSE ) {
      
      
      print("Processing...")

      sequence_list <- aa_protein_construct

      seq_names <- strsplit(as.character(cross_info_table_ab$crossreactive_vector_i[i]), ".", fixed = TRUE)
      seq_names <- str_trim( unlist(seq_names)[1], "left")
      
      
      write.fasta(sequences = list(sequence_list), names = seq_names, file.out = "proteins_fastafile", open = "a", nbchar = 60, as.string = FALSE )
      
      used_targets <- c(used_targets, crossreactive_target_i)
      
      
    }
    
  }
  
  
  if( nchar(aa_protein_construct) > 12 ) {
    
    prest_file <- "PrEST_fastafile"
    protein_file <- "proteins_fastafile"
    
    prest_aa <- readAAStringSet(prest_file, format="fasta",
                          nrec=-1L, skip=0L, seek.first.rec=FALSE,
                          use.names=TRUE, with.qualities=FALSE)
    
    makeblastdb(protein_file, dbtype = "prot")
    
    bl <- blast(db = "proteins_fastafile", type = "blastp")  # Protein + crosscaptured # Add path to the file
    
    #Run BLAST query:
    blast_result <- predict(bl, prest_aa) 
    
    
    
    
    if(nrow(blast_result) >=  1) {
      
      #Add Binder info as QueryID
      blast_result$QueryID <- ab
      
      ### Add name for the binder PrEST:
      blast_result$binder_prest <- as.character(PrEST_information[grep(ab, PrEST_information$unique_binder_name), 9])
      
      ### Add a list counting how many targets were not shown in the blast_results (to unsimilar)
      
      
      blast_result_big_table <- rbind(blast_result_big_table, blast_result)
    
    }

    
  } else  {
    cross_reactive_targets_no_aa <- c(cross_reactive_targets_no_aa , crossreactive_target_i)
    
  }
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###                                                                                                  Load expression level table:
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load GPCR expression level information
seperator <- ","
expression_level <- read.csv("gpcr_expression_level_file.csv", sep=seperator, dec=",", stringsAsFactors=F, header=T)

#Remove CALCRL that is the positive control since this one is excluded from cross reactivity analysis
expression_level <- expression_level[-grep("^CALCRL.posctrl$", expression_level$GPCR),]

## Adjust names of GPCRs so they become common between lists:
reactivity_information$gene_name[ grep("NPY4R;NPY4R2" , reactivity_information$gene_name)] <- "NPY4R"
reactivity_information$gene_name[ grep("PPAN-P2RY11;P2RY11" , reactivity_information$gene_name)] <- "P2RY11"
reactivity_information$gene_name[ grep("HCAR2;HCAR3" , reactivity_information$gene_name)] <- "HCAR2"


#################===================================================================================================================================================================
##                                  Calculate the expression ratio of all GPCRs within each subfamilyto get a representative distribution of expresison ratio
#################===================================================================================================================================================================
sba_names <- c("alpha", "beta", "gamma", "delta", "gsaf", "other")

gpcr_expression_table <- NULL
for (sba in sba_names){
  
  #### Extract all the GPCRs that are within the SBA:
  gpcr_protein_name <- unique(reactivity_information$gene_name[ grep(sba, reactivity_information$sba)])
  
  #########           Calculate the expression ratio between all proteins:
  
  for (gpcr in gpcr_protein_name){
    #Remove the gpcr we are using now so we dont get protein1/protein1
    gpcr_protein_name_i <- gpcr_protein_name[-grep(gpcr, gpcr_protein_name )]
    
    
    for ( other_gpcr in gpcr_protein_name_i){
      
      average_expression_ratio <-  as.numeric(expression_level[ grep( paste0( "^",  gpcr, "$")  , expression_level$GPCR) ,  2 ]) / as.numeric(expression_level[ grep( paste0("^", other_gpcr, "$") , expression_level$GPCR) ,  2 ])
      protein_ratio_name <- paste0(gpcr, "_", other_gpcr)
      
      input_vector <- cbind(sba, protein_ratio_name, average_expression_ratio)
      
      gpcr_expression_table <- rbind(gpcr_expression_table, input_vector)
      
      
    }
  }
}

gpcr_expression_table <- data.frame(gpcr_expression_table)

gpcr_expression_table$average_expression_ratio <- as.numeric(as.character(gpcr_expression_table$average_expression_ratio))


###################### Keep only rows with expresison ratio > 1 

gpcr_expression_table_above1 <- subset(gpcr_expression_table, average_expression_ratio > 1) 

###### Make one row per cross-reactive sample
selected_binders <- crossreacrive_abs$unique_binder_name

crossreactivity_information_table <- NULL

for(i in selected_binders){
  
  crossreacrive_abs$crossreact_names_unique[grep(i, crossreacrive_abs$unique_binder_name)]
  
  ab_position <- grep(i, crossreacrive_abs$unique_binder_name)
  
  crossreactive_vector_i <- unlist(strsplit(crossreacrive_abs$crossreact_names_unique[grep(i, crossreacrive_abs$unique_binder_name)], ","))
  
  crossreactivity_information_i  <- cbind(crossreacrive_abs[     ab_position, 1:9], crossreactive_vector_i )
  
  crossreactivity_information_table <- rbind(crossreactivity_information_table, crossreactivity_information_i)
  
  
}

########### Remove sample information so that only GPCR name is left, then count how many times it appears in the list:

crossreactivity_information_table$crossreactive_gpcr <- NA
crossreactivity_information_table$crossreactive_vector_i <- gsub(" ", "", crossreactivity_information_table$crossreactive_vector_i, fixed = TRUE)

for ( r in (1:nrow(crossreactivity_information_table))){
  gpcr_cross <- as.character(crossreactivity_information_table$crossreactive_vector_i[r])
  
  crossreactivity_information_table$crossreactive_gpcr[r] <- gsub("\\..*","", gpcr_cross)
  
}
###################         Calculate ratio between expression level of target / expression level of crossreactive protein:



crossreactivity_information_table$expression_ratio <- NA

## Change the ones with too long gene name to the same as in the expression table:
crossreactivity_information_table$gene_name[grep("NPY4R;NPY4R2", crossreactivity_information_table$gene_name)] <- "NPY4R"
crossreactivity_information_table$gene_name[grep("PPAN-P2RY11;P2RY11", crossreactivity_information_table$gene_name)] <- "P2RY11"

## Change name of all cross-reactive CALCRL.posctrl CALCRL:
expression_level$GPCR[grep("CALCRL", expression_level$GPCR)] <- "CALCRL"


for (i in 1:nrow(crossreactivity_information_table)){
  correct_target <- as.character(crossreactivity_information_table$gene_name[i])
  
  cross_reactive_protein <- as.character(crossreactivity_information_table$crossreactive_gpcr[i])

  ### TAKE off target/ target
  
  average_expression_ratio <-  as.numeric(expression_level[ grep( paste0( "^",  cross_reactive_protein, "$")  , expression_level$GPCR) ,  2 ]) / as.numeric(expression_level[ grep( paste0("^", correct_target, "$") , expression_level$GPCR) ,  2 ])
  
  average_expression_ratio_otherway <- as.numeric(expression_level[ grep( paste0("^", correct_target, "$") , expression_level$GPCR) ,  2 ]) / as.numeric(expression_level[ grep( paste0( "^",  cross_reactive_protein, "$")  , expression_level$GPCR) ,  2 ])
  
  
  crossreactivity_information_table$expression_ratio[i] <- average_expression_ratio
  
  crossreactivity_information_table$expression_ratio_otherway[i] <- average_expression_ratio_otherway
}

crossreactivity_information_table$expression_ratio_log2 <- log2(crossreactivity_information_table$expression_ratio)
crossreactivity_information_table$expression_ratio_otherway_log2 <- log2(crossreactivity_information_table$expression_ratio_otherway)

##################################### Repeat the density plot but include the off target expression as it's own

crossreactivity_information_table$protein_ratio_name <- paste0(crossreactivity_information_table$crossreactive_gpcr, "_", crossreactivity_information_table$gene_name)
crossreactivity_information_table_original_sba <- crossreactivity_information_table
crossreactivity_information_table$sba <- "crossreactive"

gpcr_crossreactivity_short <- crossreactivity_information_table[,c(4, 16, 12)]

## Only keep unique ones:
gpcr_crossreactivity_short <- distinct(gpcr_crossreactivity_short) 

colnames(gpcr_crossreactivity_short) <- colnames(gpcr_expression_table_above1)
expression_ratio_table_all_proteins <- rbind(gpcr_expression_table_above1, gpcr_crossreactivity_short)

gpcr_expression_table_everyone <- rbind(gpcr_expression_table, gpcr_crossreactivity_short)


########## Prepare a histogram for how many in each sba crossreacts

table(crossreacrive_abs$sba)

all_abs_table_wo_controls <- reactivity_information[-grep("ctrl", reactivity_information$sba), ]

df_melt_all_abs <- melt(table(all_abs_table_wo_controls$sba))

df_melt <- melt(table(crossreacrive_abs$sba))

df_percentage <- cbind(df_melt_all_abs, df_melt)
colnames(df_percentage) <- c("sba", "all", "sba2", "cross")
df_percentage$procent <- df_percentage$cross/df_percentage$all*100
df_percentage$cross_text <- paste0("N=", df_percentage$cross)

theme_set(theme_classic(base_size = 12))

p<-ggplot(data=df_percentage, aes(x=sba, y=procent)) +
  geom_bar(stat="identity", fill = "gray")+
  ggtitle("% crossreactive Abs per SBA") +
  xlab("SBA") + ylab("%") +
  ylim(0, 100)+
  geom_text(aes(label=cross_text), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(axis.line = element_line(size = 1) , axis.text.y = element_text(size = 20, color = "black"), axis.text.x = element_text(size = 15, color = "black") , axis.title = element_text(size = 16, face = "bold"), plot.title = element_text(size = 20, face = "bold"), legend.title = element_text(size = 16), axis.ticks.length.x = unit(0.3, "cm"), axis.ticks.length.y = unit(0.3, "cm") ) 


#pdf("percentage crossreactive antibodies per sba_221021_2.pdf", width = 5, height = 5, useDingbats=F)
print(p)
#dev.off()

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                     Density plots with log2 expression ratio
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### Add crossreactive proteins the other way around: (target / crossreactive protein)

gpcr_cross_otherway <- crossreactivity_information_table[,c(4, 16, 13, 15)]

### Keep only unique ones:
#gpcr_cross_otherway <- distinct(gpcr_cross_otherway) 

colnames(gpcr_cross_otherway) <- colnames(gpcr_expression_table_everyone)

gpcr_expression_table_everyone_other_way_too <- rbind(gpcr_expression_table_everyone, gpcr_cross_otherway)

expression_ratio_table_all_proteins$expression_ratio_log2 <- log2(expression_ratio_table_all_proteins$average_expression_ratio)
gpcr_expression_table_everyone$expression_ratio_log2 <- log2(gpcr_expression_table_everyone$average_expression_ratio)

gpcr_expression_table_everyone_other_way_too$expression_ratio_log2 <- log2(gpcr_expression_table_everyone_other_way_too$average_expression_ratio)

gpcr_expression_table_everyone_wo_cross <- gpcr_expression_table_everyone[-grep("crossreactive", gpcr_expression_table_everyone$sba), ]


for (combination in 1:nrow(gpcr_cross_otherway)){
  combination_i <- gpcr_cross_otherway$protein_ratio_name[combination]
  
  if (combination_i %in% gpcr_expression_table_everyone_wo_cross$protein_ratio_name){
    
    print (combination_i)
    gpcr_expression_table_everyone_wo_cross <- gpcr_expression_table_everyone_wo_cross[!grepl(combination_i, gpcr_expression_table_everyone_wo_cross$protein_ratio_name),]
    
  }
}


gpcr_cross_otherway$expression_ratio_log2 <- log2(gpcr_cross_otherway$average_expression_ratio)

gpcr_expression_table_everyone_other_way_too_only_cross_in_cross <- rbind(gpcr_expression_table_everyone_wo_cross, gpcr_cross_otherway, gpcr_expression_table_everyone[grep("crossreactive", gpcr_expression_table_everyone$sba), ])



##                                                   Prepare density plots within each subfamily, use cross reactives (both ways) and only in cross reactive group

plot_list_sba_density <- list()
for ( sba in sba_names){
  gpcr_expression_table_everyone_sba_i <- gpcr_expression_table_everyone_other_way_too_only_cross_in_cross[grep(sba, gpcr_expression_table_everyone_other_way_too_only_cross_in_cross$sba), ]
  gpcr_expression_table_everyone_sba_i$crossreactive <- "False"
  
  crossreactivity_sba_i <- crossreactivity_information_table_original_sba[grep(sba, crossreactivity_information_table_original_sba$sba), ]
  crossreactivity_sba_i$crossreactive <- "True"
  crossreactivity_sba_i_otherway <- crossreactivity_sba_i[, c(4, 16, 13, 15, 17)]
  crossreactivity_sba_i <- crossreactivity_sba_i[, c(4, 16, 12, 14, 17)]
  crossreactivity_sba_i_otherway <- crossreactivity_sba_i_otherway[!duplicated(crossreactivity_sba_i_otherway$protein_ratio_name),]
  crossreactivity_sba_i <- crossreactivity_sba_i[!duplicated(crossreactivity_sba_i$protein_ratio_name),]
  
  colnames(crossreactivity_sba_i) <- colnames(gpcr_expression_table_everyone_sba_i)
  colnames(crossreactivity_sba_i_otherway) <- colnames(gpcr_expression_table_everyone_sba_i)
  
  gpcr_expression_ratio_table_all_abs_cross_reactive_sba_i<- rbind( gpcr_expression_table_everyone_sba_i,  crossreactivity_sba_i,crossreactivity_sba_i_otherway )
  
  
  mu <- ddply(gpcr_expression_ratio_table_all_abs_cross_reactive_sba_i, "crossreactive", summarise, grp.mean=mean(expression_ratio_log2))
  head(mu)
  p_i <- ggplot(gpcr_expression_ratio_table_all_abs_cross_reactive_sba_i, aes(x=expression_ratio_log2, color=crossreactive)) +
    geom_density( )+
    ggtitle(paste0("Expression ratio density - ", sba )) +
    geom_vline(data=mu, aes(xintercept=grp.mean, color=crossreactive),
               linetype="dashed") +
    scale_color_manual(values=c(colorpalette_dark2[1:5], "black"))
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1), axis.line = element_line(size = 1) , axis.text = element_text(size = 10) , axis.title = element_text(size = 23, face = "bold"), plot.title = element_text(size = 45, face = "bold"), legend.title = element_text(size = 16), axis.ticks.length.x = unit(0.3, "cm"), axis.ticks.length.y = unit(0.3, "cm") ) 
  
  plot_list_sba_density[[sba]] <- p_i
  
}



spl <- split(plot_list_sba_density, (seq_along(plot_list_sba_density)-1) %/% 6)

ppl <- lapply(spl, function(g)
  marrangeGrob(grobs = g, layout_matrix = matrix(1:6, ncol=3, byrow=T), top=NULL))


#pdf(paste0("Density plots per subfamily", format(Sys.time(),"%Y-%m-%d_%H%M%S"),
#           ".pdf"),
#    width = 15, height = 9, useDingbats=F, onefile=T)

print(ppl)

#dev.off()



crossreactivity_information_table$log_expression_ratio <- log(crossreactivity_information_table$expression_ratio)
crossreactivity_information_table$log2_expression_ratio <- log2(crossreactivity_information_table$expression_ratio)
# 
# 
expression_ratio_table <- crossreactivity_information_table[ , c(1, 6, 7, 9:14)  ]
expression_ratio_table$target_crossreactiveprotein <- NULL
expression_ratio_table$antibody_crossreactiveprotein <- NULL
for ( i in 1:nrow(expression_ratio_table)){
  expression_ratio_table$target_crossreactiveprotein[i] <- paste0(expression_ratio_table$gene_name[i], ".", expression_ratio_table$crossreactive_gpcr[i])
  expression_ratio_table$antibody_crossreactiveprotein[i] <- paste0(expression_ratio_table$antibody[i],  ".", expression_ratio_table$gene_name[i], ".", expression_ratio_table$crossreactive_gpcr[i])
  
}
# 
# ### Remove replicated expression ratios:
expression_ratio_table_unique <- expression_ratio_table[!duplicated(expression_ratio_table$target_crossreactiveprotein),]
# 

E_value_table_PrEST <- blast_result_big_table

### ###Add column with binder crossreactive protein to both E_value_table_PrEST and expression_ratio_table:

E_value_table_PrEST$binder_crossprotein <- NA

E_value_table_PrEST$binder_crossprotein <- paste0(E_value_table_PrEST$QueryID,".", E_value_table_PrEST$SubjectID)

expression_ratio_table$E_value_PrESTvsProtein <- NA

for (row in 1:nrow(expression_ratio_table)){
  search_for <- paste0(expression_ratio_table$unique_binder_name[row],".", expression_ratio_table$crossreactive_gpcr[row])
  if (length(grep(search_for,E_value_table_PrEST$binder_crossprotein)) > 0 ){
    expression_ratio_table$E_value_PrESTvsProtein[row] <- E_value_table_PrEST$E[grep(search_for,E_value_table_PrEST$binder_crossprotein )]
  } else {
    expression_ratio_table$E_value_PrESTvsProtein[row] <-NA
  }
  
}

##### Add also the e-value to crossreactivity_information_table

crossreactivity_information_table$E_value_PrESTvsProtein <- NA

for (row in 1:nrow(crossreactivity_information_table)){
  search_for <- paste0(crossreactivity_information_table$unique_binder_name[row],".", crossreactivity_information_table$crossreactive_gpcr[row])
  if (length(grep(search_for,E_value_table_PrEST$binder_crossprotein)) > 0 ){
    crossreactivity_information_table$E_value_PrESTvsProtein[row] <- E_value_table_PrEST$E[grep(search_for,E_value_table_PrEST$binder_crossprotein )]
  } else {
    crossreactivity_information_table$E_value_PrESTvsProtein[row] <-NA
  }
  
}


crossreactivity_information_table$expression_ratio_log2 <- log2(crossreactivity_information_table$expression_ratio)
crossreactivity_information_table$target_cross <- paste0(crossreactivity_information_table$gene_name, ".", crossreactivity_information_table$crossreactive_gpcr)


##############################       Add E-value binary level column to crossreactivity_information_table
crossreactivity_information_table$prest_e_value_binary <- crossreactivity_information_table$E_value_PrESTvsProtein
crossreactivity_information_table$prest_e_value_binary[is.na(crossreactivity_information_table$prest_e_value_binary)] <- 0

for (row in 1:nrow(crossreactivity_information_table)){
  if (crossreactivity_information_table$prest_e_value_binary[row] > 0) {
    crossreactivity_information_table$prest_e_value_binary[row] <- 1
  } else {
    crossreactivity_information_table$prest_e_value_binary[row] <- 0
  }
}

crossreactivity_information_table_unique_crosstarget <- crossreactivity_information_table[!duplicated(crossreactivity_information_table$protein_ratio_name),]

crossreactivity_information_table_unique_crosstarget$prest_e_value_binary <- as.factor(crossreactivity_information_table_unique_crosstarget$prest_e_value_binary)


################################################################################################################################################################################################################

## Add manually the E-value for the CAB that was made from web
e_value_table_CAB <- read_excel("e_value_for_cab.xlsx")

## Target = GPR25, cross protien = PRLHR
crossreactivity_information_table_unique_crosstarget$E_value_PrESTvsProtein[grep("DF4968", crossreactivity_information_table_unique_crosstarget$antibody)] <- e_value_table_CAB$E[grep("PRLHR",e_value_table_CAB$SubjectID )]

crossreactivity_information_table_unique_crosstarget$E_value_PrESTvsProtein_dotsize <- NA

rownames(crossreactivity_information_table_unique_crosstarget) <- 1:nrow(crossreactivity_information_table_unique_crosstarget)

for (row in 1:nrow(crossreactivity_information_table_unique_crosstarget)){
  e_value_i <- crossreactivity_information_table_unique_crosstarget$E_value_PrESTvsProtein[row]
  
  if (is.na(e_value_i)){
    crossreactivity_information_table_unique_crosstarget$E_value_PrESTvsProtein_dotsize[row] <- 10
    
  } else if ( e_value_i < 1 && e_value_i >  0.1  ) {
    
    crossreactivity_information_table_unique_crosstarget$E_value_PrESTvsProtein_dotsize[row] <- 30
    
  } else if ( e_value_i < 10 && e_value_i >  1  ) {
    
    crossreactivity_information_table_unique_crosstarget$E_value_PrESTvsProtein_dotsize[row] <- 20
    
  } else if (e_value_i < 0.1 ) {
    crossreactivity_information_table_unique_crosstarget$E_value_PrESTvsProtein_dotsize[row] <-40
  }
  
}

crossreactivity_information_table_unique_crosstarget$E_value_PrESTvsProtein_dotsize <- as.character(crossreactivity_information_table_unique_crosstarget$E_value_PrESTvsProtein_dotsize)

################    Add a specific color column to expression_ratio_table so we can have a different colors for the GPCRs that are captured in more than one sample:

expression_ratio_table$colorcolumn <- 1

expression_ratio_table$colorcolumn[duplicated( expression_ratio_table$target_crossreactiveprotein)] <- 2
expression_ratio_table$colorcolumn <- as.factor(expression_ratio_table$colorcolumn)

crossreactivity_information_table_unique_crosstarget$target_crossreactiveprotein <- paste0(crossreactivity_information_table_unique_crosstarget$gene_name, "." , crossreactivity_information_table_unique_crosstarget$crossreactive_gpcr)

crossreactivity_information_table_unique_crosstarget$colorcolumn <- 1
for (row in 1:nrow(crossreactivity_information_table_unique_crosstarget)){
  if (crossreactivity_information_table_unique_crosstarget$target_crossreactiveprotein[row] %in% expression_ratio_table$target_crossreactiveprotein[duplicated(expression_ratio_table$target_crossreactiveprotein)] ){
    crossreactivity_information_table_unique_crosstarget$colorcolumn[row] <- 2
    
  }
}


crossreactivity_information_table_unique_crosstarget$colorcolumn <- as.factor(crossreactivity_information_table_unique_crosstarget$colorcolumn)

expression_ratio_evalue_plot <- ggplot(data = crossreactivity_information_table_unique_crosstarget, aes( x = reorder(target_crossreactiveprotein, - expression_ratio_log2), y = expression_ratio_log2, size=E_value_PrESTvsProtein_dotsize)) +
  geom_point(aes(size=E_value_PrESTvsProtein_dotsize, fill = colorcolumn) ,stat = "identity", pch = 21, alpha = 0.75) + 
  scale_fill_manual(values = c("#F8766D", "#00BFC4"), labels = c( "1" = "1", "2" = "> 1"), name = "Captured") +
  ggtitle(paste0("Expression ratio of cross reactive GPCR to target GPCR")) +
  labs(y = "log2(Expression ratio)", x = "Target.Crossreactive protein") +
  geom_hline(yintercept=1, linetype="dashed", color = "black", lwd = 1) + 
  geom_hline(yintercept=0, linetype="dotted", color = "black", lwd = 1) + 
  geom_hline(yintercept=-1, linetype="dashed", color = "black", lwd = 1) + 
  scale_y_continuous(n.breaks = 6) +
  scale_size_manual(name   = "E-value" ,
                    values = c("10" = 4 , "20" = 8 , "30" = 12, "40" = 16),
                    labels = c("10"= "E > 10", "20" = "1 < E < 10", "30" ="0.1 < E < 1", "40" =" E < 0.1" ) )  + 
  guides(fill = guide_legend(override.aes = list(size = c(6, 6)))) + 
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1), axis.line = element_line(size = 1) , axis.text = element_text(size = 10) , axis.title = element_text(size = 15, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold"), legend.title = element_text(size = 15), axis.ticks.length.x = unit(0.3, "cm"), axis.ticks.length.y = unit(0.3, "cm"),
        axis.text.y = element_text(size = 20)) 

#pdf("scatter plot Expression ratio and Evalue_PrESTvsProtein_20221024.pdf", width = 10, height = 7, useDingbats=F)
print(expression_ratio_evalue_plot)
#dev.off()








