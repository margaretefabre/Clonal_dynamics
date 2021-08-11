# Adapted from Mike Spencer Chapman's script

#This script is designed to be interactive & quick to allow rapid exploring of filtering parameters

# on farm5:
# module load mpboot
# module load vagrent

library(stringr)
library(ape)
library(seqinr)
library(tidyverse)
library(data.table)
library(phangorn)


my_working_directory="/lustre/scratch119/casm/team163gv/mf13/COLONIES/PD34493/PD34493_1/"	#Edit

treemut_dir="/nfs/casm/team163gv/mf13/COLONIES/treemut/treemut/"	#Edit
setwd(my_working_directory)

R_function_files = list.files("/nfs/casm/team163gv/mf13/COLONIES/my_functions/my_functions/",pattern=".R",full.names=TRUE)	#Edit
sapply(R_function_files,source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

#Set the file paths for saved files
Run_ID = "PD34493"	#Edit
filtering_ID = "standard"	#Edit
mats_and_params_file = paste0("mats_and_params_", Run_ID)
filtered_muts_file = paste0("filtered_muts/filtered_muts_",Run_ID,"_",filtering_ID)
dna_string_file = paste0("mpboot_files/Filter_", Run_ID,"_", filtering_ID,".fa")
mpboot_tree_file = paste0(dna_string_file,".treefile")
tree_file_path = paste0("tree_files/tree_", Run_ID,"_",filtering_ID, ".tree")
vcf_header_path = "mutation_vcfs/VCF_header_for_VaGrent.txt"
vcf_path = paste0("mutation_vcfs/mutations_all_", Run_ID,"_",filtering_ID,".vcf")
shared_vcf_path = paste0("mutation_vcfs/mutations_shared_", Run_ID,"_",filtering_ID,".vcf")
vagrent_input_path = paste0("mutation_vcfs/mutations_all_", Run_ID,"_",filtering_ID,"_header.vcf")
vagrent_output_path = paste0(vagrent_input_path,".annot")
file_annot = paste0("annotated_muts/annotated_mut_set_", Run_ID,"_",filtering_ID)	#Edit
CHgenes_path = "/nfs/casm/team163gv/mf13/COLONIES/files/CHgenes"	#Includes 56 genes in CH targeted sequencing bait set
CancerCensusGenes_path = "/nfs/casm/team163gv/mf13/COLONIES/files/CancerCensusGenes"	#Includes tier 1 genes in Cancer Gene Census, n=577

sensitivity_analysis_path <- "sensitivity_analysis_PD34493" #Output table from the sensitivity analysis	#Edit
#####


#If have previously run analysis - load up files, or can rerun the analysis
previously_run = F	#Edit
if(previously_run) {
  load(file_annot); tree <- read.tree(tree_file_path)
} else {
  load(mats_and_params_file )
  
  
#Can specify coverage cut-off here, to exclude low coverage samples from analysis - samples too difficult to even place on tree
  min_sample_mean_cov = 4	#Edit


#Remove mixed samples  
mixed_samples_already_excluded=FALSE	#Edit
if(mixed_samples_already_excluded) {
  mixed_samples=NULL
} else {
  colnames(COMB_mats$NV) <- gsub(pattern = "_MTR", replacement = "",x = colnames(COMB_mats$NV))
  sample_peak_vaf = sapply(colnames(COMB_mats$NV), check_peak_vaf, COMB_mats=COMB_mats, filter_params = filter_params)
  
  #Exclude samples with a peak VAF < 0.4 #suggests mixed colony 
  mixed_samples = names(sample_peak_vaf[sample_peak_vaf <0.4])	#Edit
  pass_samples = names(sample_peak_vaf[sample_peak_vaf >=0.4])	#Edit
}


#Visualize the density plots of the mixed and pass samples for sanity check .
pdf(paste0(my_working_directory, "plots/mixedSample_FAIL.pdf", Run_ID, ".pdf"),width=15,height=7)	#Edit
par(mfrow=c(3,2))
sapply(mixed_samples, vaf_density_plot, COMB_mats, filter_params)
dev.off()
pdf(paste0(my_working_directory, "plots/mixedSample_PASS.pdf", Run_ID, ".pdf"),width=15,height=7)	#Edit
par(mfrow=c(3,2))
sapply(pass_samples, vaf_density_plot, COMB_mats, filter_params)
dev.off()


#Remove (i) low coverage samples and (ii) mixed samples, and their private mutations
if(min_sample_mean_cov > 0|!is.null(mixed_samples)) {
  output = remove_low_coverage_samples(COMB_mats = COMB_mats,
                                       filter_params = filter_params,
                                       min_sample_mean_cov = min_sample_mean_cov,
                                       other_samples_to_remove = mixed_samples,
                                       min_variant_reads_auto = 3, #these parameters are to remove mutations from the matrix that are no longer positive in any samples
                                       min_variant_reads_xy = 2)
  COMB_mats= output$COMB_mats
  filter_params = output$filter_params
}

#Removing 3 samples:
#	0 low depth 
#	3 mixed samples (PD34493k_42, PD34493k_49, PD34493k_75)
#	0 duplicate 
#	3672 mutations removed as no positives in any remaining samples
  
  
#REVIEW MEAN DEPTH HISTOGRAMS TO DECIDE MANUAL CUT-OFFS FOR EXCLUDING OUTLIERS
#hist(filter_params$mean_depth, breaks = 100, xlim = c(0,60))
#hist(filter_params$mean_depth[COMB_mats$mat$Chrom %in% c("X","Y")], breaks = 200, xlim = c(0,30))
#hist(filter_params$mean_depth[!COMB_mats$mat$Chrom %in% c("X","Y")], breaks = 200, xlim = c(0,30))
XY_low_depth_cutoff = 3; XY_high_depth_cutoff = 11; AUTO_low_depth_cutoff = 9; AUTO_high_depth_cutoff = 22 	#Edit


#Get the filtered mutation set - returns object with the filtered mat/ NV/NR matrices, as well as a full genotype matrix, and matrix of shared mutations only
filtered_muts = get_filtered_mut_set(input_set_ID = Run_ID,  #the Run_ID of the unfiltered mutation set used as input - just gets stored with output as a record
                                     COMB_mats = COMB_mats,  #the main full mutation matrix as outputed by the "HSC_filtering_treebuild_table.R" script
                                     filter_params = filter_params,  #the filter_params matrix as outputed by the "HSC_filtering_treebuild_table.R" script
                                     gender = COMB_mats$gender, #patients gender. Must be "male" or "female". This is determined in the "HSC_filtering_treebuild_table.R" script
                                     retain_muts = NA,  #the mut_refs (i.e. Chr-Pos-Ref-Alt) of any mutations that should be manually retained, despite not meeting filtering criteria, NULL by default
                                     
                                     #PARAMETERS FOR FILTERING THE MUTATION SET. If parameter is set to null, filter will not be applied.
                                     germline_pval = -10,  #the log10 p-value cutoff for mutations coming from an expected germline distribution
                                     rho = 0.1,  #rho cutoff for the beta-binomial filter
                                     mean_depth = c(AUTO_low_depth_cutoff,AUTO_high_depth_cutoff, XY_low_depth_cutoff, XY_high_depth_cutoff),   #Numeric vector of length 4 defining mean depth at mutation site cut-offs. This is in the order (1) lower threshold for autosomes, (2) upper threshold for autosomes, (3) lower threshold for XY, (4) upper threshold for XY. This removes mis-mapping/ low reliability loci.
                                     min_depth = c(6,4), #Numeric vector of length 2 defining minimum depths that at least one positive sample must have for mutation to be retained (AUTO and XY)
                                     min_vaf = NA, #Numeric vector of length 2 defining minimum vaf in at least one sample for mutation to be retained (AUTO and XY). Mutations with value > cutoff are retained.
                                     
                                     #These filters are only appropriate for pure clonal samples i.e. single-cell colonies
                                     pval_dp2 = NA,  #the p-value cut-off if using the "pval within pos" filter, with positive samples defined as having >= 2 reads
                                     pval_dp3 = 0.001,   #the p-value cut-off if using the "pval within pos" filter, with positive samples defined as having >= 3 reads (allows for more index hopping)
                                     min_pval_for_true_somatic = 0.1,   #the minimum p-value that at least one sample must have for the variant:normal read distribution coming from that expected for a true somatic
                                     
                                     #PARAMETERS FOR DECIDING GENOTYPE FOR EACH SAMPLE FOR EACH RETAINED MUTATION - should be equal to, or more relaxed than the above. At least one parameter must be used.
                                     min_variant_reads_SHARED = 2,  #the minimum number of reads for subsequent samples to be assigned a positive genotype
                                     min_pval_for_true_somatic_SHARED = 0.05, #the minimum p-value for coming from "true somatic mutation" read distribution for subsequent samples to be assigned a positive genotype
                                     min_vaf_SHARED = NA #Numeric vector of length 2, defining minimum vaf to be assigned a positive genotype


  
#Decide an ID for this filtered set, depending on approach taken, and save
save(filtered_muts, file = filtered_muts_file)

#Write a fasta file of the dummy DNA strings
write.fasta(filtered_muts$dna_strings, names=names(filtered_muts$dna_strings), dna_string_file)

#BUILD TREE with MPBoot
system(paste0("mpboot -s ", dna_string_file," -bb 1000"))

#Import the tree into R using ape
tree <- read.tree(mpboot_tree_file)
tree <- root(tree, "Ancestral", resolve.root = TRUE)
tree <- multi2di(tree)
tree$edge.length = rep(1, nrow(tree$edge)) #Initially need to assign edge lengths of 1 for the tree_muts package to work
  
#Assign mutations back to the tree using treemut package
df = reconstruct_genotype_summary(tree) #Define df (data frame) for treeshape

mtr = filtered_muts$COMB_mats.tree.build$NV
mtr$Ancestral = 0
mtr = as.matrix(mtr)
depth = filtered_muts$COMB_mats.tree.build$NR
depth$Ancestral = 10
depth = as.matrix(depth)
p.error = c(rep(0.01, ncol(filtered_muts$COMB_mats.tree.build$NR)),1e-6)
res = assign_to_tree(mtr=mtr[,df$samples],dep=depth[,df$samples], df=df, error_rate=p.error)
treefit_pval_cutoff = 1e-3
poor_fit = res$summary$pval < treefit_pval_cutoff  #See how many mutations don't have read counts that fit the tree very well
sum(poor_fit)

#Remove poor_fit muts if they look erroneous
filter_poor_fit = "no"
if(filter_poor_fit == "yes") {
  filtered_muts$COMB_mats.tree.build <- list_subset(filtered_muts$COMB_mats.tree.build,select_vector = !poor_fit)
  mtr = as.matrix(filtered_muts$COMB_mats.tree.build$NV)
  depth = as.matrix(filtered_muts$COMB_mats.tree.build$NR)
  p.error = c(rep(0.01, ncol(filtered_muts$COMB_mats.tree.build$NR)))
  res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
}
  
#Add node and pval information to the filtered_muts object
filtered_muts$COMB_mats.tree.build$mat$node <- tree$edge[res$summary$edge_ml,2]
filtered_muts$COMB_mats.tree.build$mat$pval <- res$summary$pval

#Assign edge lengths to the tree - safest to use this approach to assignment if have filtered any poor fit mutations
tree$edge.length <- sapply(tree$edge[,2],function(node) {sum(filtered_muts$COMB_mats.tree.build$mat$node==node)})
  

#Save the filtered_muts file and the tree file
write.tree(tree, file = tree_file_path)
save(filtered_muts, file = paste0(filtered_muts_file, "_beforeAnnotation"))
  
#----------------------------------------------- 
#ANNOTATING THE MUTATIONS
#-----------------------------------------------
#Write vcf files for VariantCaller analysis - separate out "All mutations" & Shared mutations"
vcf_file = create_vcf_files(filtered_muts$COMB_mats.tree.build$mat)
write.table(vcf_file, sep = "\t", quote = FALSE, file = vcf_path, row.names = FALSE)

#Also save the shared mutations for mutational signature analysis
shared_vcf_file = create_vcf_files(filtered_muts$COMB_mats.tree.build$mat,select_vector = which(!filtered_muts$COMB_mats.tree.build$mat$node%in%1:length(tree$tip.label)))
write.table(shared_vcf_file, sep = "\t", quote = FALSE, file = shared_vcf_path, row.names = FALSE)

#1. paste vcf file to a dummy header file
system(paste0("cat ",vcf_header_path," ",vcf_path," > ", vagrent_input_path))

#2. commands to run vagrent
system(paste0("AnnotateVcf.pl -i ",vagrent_input_path," -o ",vagrent_output_path," -sp Human -as NCBI37 -c /lustre/scratch117/casm/team78pipelines/reference/human/GRCh37d5/vagrent/e75/vagrent.cache.gz"))
   
#3. import vagrent output and add into the filtered_muts object
vagrent_output = fread(vagrent_output_path,skip = "#CHROM")
annot_info = as.data.frame(str_split(vagrent_output$INFO, pattern = ";",simplify = TRUE), stringsAsFactors = FALSE)
colnames(annot_info) <- c("VT","VD","VC","VW")

annot_info$VC <- gsub(x=annot_info$VC, pattern = "VC=", replacement = "")
annot_info$VT <- gsub(x=annot_info$VT, pattern = "VT=", replacement = "")
annot_info$VW <- gsub(x=annot_info$VW, pattern = "VW=", replacement = "")
annot_info$VD <- gsub(x=annot_info$VD, pattern = "VD=", replacement = "")

filtered_muts$COMB_mats.tree.build$mat <- cbind(filtered_muts$COMB_mats.tree.build$mat,split_vagrent_output(df = annot_info,split_col = "VD"))  

#Save the annotated filtered_muts files (post tree filtering)
save(filtered_muts, file = file_annot)
}

#----------------------------------------------- 
#Additional check for mixed colonies
#----------------------------------------------- 

#Maps the VAFs of all the mutations on the tree (branch by branch) of individual colonies.  
#Flagged if inconsistent with being either (1) clonal (with a degree of flexibility - here set to at least 85%)  or (2) absent.  

details=filtered_muts$COMB_mats.tree.build$mat
matrices=list(mtr=filtered_muts$COMB_mats.tree.build$NV,dep=filtered_muts$COMB_mats.tree.build$NR)

#Test branch VAFs for each sample to see if there are any that are inconsistent with a clonal sample
sample_test=lapply(tree$tip.label,function(samples,min.depth=1) {
  print(samples)
  node_test=sapply(tree$edge[,2],function(node) {
    #print(node)
    info=get_edge_info(tree,details,node)
    if(length(info$idx.in.details)>0) {
      df=data.frame(mtr=sum(matrices$mtr[info$idx,samples],na.rm = TRUE),
                    dep=sum(matrices$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
      df=df[which(df$dep>=min.depth),]
      df$vaf=df$mtr/df$dep
      df=df[which(!is.na(df$vaf)),]
      N=dim(df)[1]
      MTR=sum(df$mtr)
      DEP=sum(df$dep)
      min.mean.vaf=0.425
      if(DEP>0) {
        z=binom.test(MTR,DEP,alternative = "less",p=min.mean.vaf)
        z2=binom.test(MTR,DEP,alternative = "greater",p=0.05)
        z$p.value=max(z$p.value,z2$p.value)
        if(z$p.value<0.05){
          if(z$p.value<0.05/dim(tree$edge)[1]){
            return(2)
          }else{
            return(1)
          }
        } else {
          return(0)
        }
      } else {
        return(NA)
      }
    } else {
      return(NA)
    }
  })
  node_test=node_test[!is.na(node_test)]
  df=data.frame(sample=samples,high_prob=sum(node_test==2),low_prob=sum(node_test==1))
  return(df)
})
sample_test_df=Reduce(rbind,sample_test)

#Define which samples have "passed" this test of clonality
sample_test_df$result=ifelse(sample_test_df$high_prob>0|sample_test_df$low_prob>2,"FAIL","PASS")
mixed_samples=as.character(sample_test_df$sample[sample_test_df$result=="FAIL"])
print(paste(mixed_samples,"will be removed as failed the test for clonality"))

#Plot the vaf trees of the failing samples
pdf(paste0(Run_ID,"_",filtering_ID,"_mixed_vaf_trees.pdf"),width=30,height = 10)
lapply(mixed_samples,function(sample) {
  tree=plot_tree(tree,cex.label=0)
  plot_tree_vaf(tree,
                details=filtered_muts$COMB_mats.tree.build$mat,
                matrices=list(mtr=filtered_muts$COMB_mats.tree.build$NV,dep=filtered_muts$COMB_mats.tree.build$NR),
                samples=sample)
})
dev.off()

#If any additional mixed colonies flagged, need to remove and re-build tree.


#---------------------------------------------------------------------------------------------------
#CORRECT THE EDGE LENGTHS BASED ON SAMPLE SENSITIVITY & INVITRO MUTATION PROPORTION
#---------------------------------------------------------------------------------------------------
#NB - for function to work correctly, the sensitivity_df must be set up in exactly the correct format

#Create trees of only SNVs or only INDELs (this function does not do any correction)
tree_SNV = get_subset_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, v.field="Mut_type", value = "SNV")
tree_INDEL = get_subset_tree(tree = tree, details = filtered_muts$COMB_mats.tree.build$mat, v.field="Mut_type", value = "INDEL")

#Import the sensitivity analysis file
sensitivity_df <- read.delim(sensitivity_analysis_path, header = TRUE, stringsAsFactors = FALSE, sep = " ")

#add a line to sensitivity_df for ancestral:
anc <- data.frame("Sample" = "Ancestral", "SNV_sensitivity" = 1, "INDEL_sensitivity" = 1)
sens1 <- bind_rows(sensitivity_df, anc)

#Create corrected SNV tree
tree_SNV_c = get_corrected_tree(tree = tree_SNV, details = filtered_muts$COMB_mats.tree.build$mat, include_indels = FALSE, sensitivity_df = sens1, get_edge_from_tree=F)



#---------------------------------------------------------------------------------------------------
#MAKE ULTRAMETRIC TREES
#---------------------------------------------------------------------------------------------------

make.ultrametric.tree <- function(tree) {
    root.number <- length(tree$tip.label) + 1
    ultra.tree <- length.normalise(tree, tree, root.number, 1)
    return(ultra.tree)
}

length.normalise <- function(orig.tree, new.tree, curr.node, remaining.stick) {
    curr.node.children <- unlist(Descendants(orig.tree, curr.node, "children"))

    for (j in curr.node.children) {
      index <- which(orig.tree$edge[,2] == j & orig.tree$edge[,1] == curr.node, arr.ind = TRUE)

      if (j %in% orig.tree$tip.label) {
        new.tree$edge.length[index] <- remaining.stick
      } else {
        curr.node.tips <- unlist(Descendants(orig.tree, j, "tips"))
        curr.dist <- find.distance(orig.tree, curr.node, j)
        if (curr.dist == 0) {curr.dist <- 0.01} # So that no edge lengths are zero
        desc.lengths <- sapply(curr.node.tips, FUN = find.distance, tree = orig.tree, from = curr.node)
        new.tree$edge.length[index] <- remaining.stick * curr.dist / mean(desc.lengths)
        shorter.stick <- remaining.stick - new.tree$edge.length[index]

        # Call as recursive function
        new.tree <- length.normalise(orig.tree, new.tree, j, shorter.stick)
      }
    }
    return(new.tree)
} 
  
find.distance <- function(tree, from, to) {
  path <- nodepath(tree, from, to)
  res <- 0
  for (j in 2:length(path)) {
    index <- which(tree$edge[,2] == path[j] & tree$edge[,1] == path[j-1], arr.ind = TRUE)
    res <- res + tree$edge.length[index]
  }
  return(res)
}  

tree_SNV_c_ultra <- make.ultrametric.tree(tree_SNV_c)


#---------------------------------------------------------------------------------------------------
#ANNOTATE MUTATIONS CAUSING CODING CHANGES IN SHARED BRANCHES (the most robust calls)
#---------------------------------------------------------------------------------------------------

details<-filtered_muts$COMB_mats.tree.build$mat
NV <- filtered_muts$COMB_mats.tree.build$NV
NR <- filtered_muts$COMB_mats.tree.build$NR

details$variant_ID <- paste(details$Gene, details$Protein, sep = " ") #Add a "variant_ID" that includes that gene & protein coding change
details$Type[details$Type == ""] <- "no_annotation" #Label blank field as "no_annotation"


#Define the coding changes in shared branches
details$shared_coding_change <- ifelse(!details$node %in% 1:length(tree$tip.label) & details$Type %in% c("protein_coding:exon:CDS:substitution:codon_variant:non_synonymous_codon",
                                                                                                         "protein_coding:exon:CDS:insertion:frameshift_variant",
                                                                                                         "protein_coding:exon:CDS:deletion:frameshift_variant",
                                                                                                         "protein_coding:exon:CDS:substitution:codon_variant:stop_gained",
                                                                                                         "protein_coding:exon:CDS:substitution:codon_variant:initiator_codon_change",
                                                                                                         "protein_coding:exon:CDS:deletion:inframe_variant:inframe_codon_loss"), "Shared coding change", "no")

#Define the CH gene coding changes in shared branches
CHgenes <- read.table(CHgenes_path, header=F, stringsAsFactors=F)
details$shared_coding_change_CHgene <- ifelse(!details$node %in% 1:length(tree$tip.label) & details$Gene %in% CHgenes$V1 & details$Type %in% c("protein_coding:exon:CDS:substitution:codon_variant:non_synonymous_codon",
                                                                                                         "protein_coding:exon:CDS:insertion:frameshift_variant",
                                                                                                         "protein_coding:exon:CDS:deletion:frameshift_variant",
                                                                                                         "protein_coding:exon:CDS:substitution:codon_variant:stop_gained",
                                                                                                         "protein_coding:exon:CDS:substitution:codon_variant:initiator_codon_change",
                                                                                                         "protein_coding:exon:CDS:deletion:inframe_variant:inframe_codon_loss"), "Shared coding change", "no")
                                       
                                                                                                      
#Define the cancer gene census (tier 1 genes only n=577) and CH gene coding changes in shared branches
#Add in cancer gene census annotation (tier 1 genes only n=577)
censusGenes <- read.table(CancerCensusGenes_path, header=F, stringsAsFactors=F)
details$shared_coding_change_CHorCensusgene <- ifelse(!details2$node %in% 1:length(tree$tip.label) & (details2$Gene %in% CHgenes$V1 | details2$Gene %in% censusGenes$V1) & details2$Type %in% c("protein_coding:exon:CDS:substitution:codon_variant:non_synonymous_codon",
                                                                                                         "protein_coding:exon:CDS:insertion:frameshift_variant",
                                                                                                         "protein_coding:exon:CDS:deletion:frameshift_variant",
                                                                                                         "protein_coding:exon:CDS:substitution:codon_variant:stop_gained",
                                                                                                         "protein_coding:exon:CDS:substitution:codon_variant:initiator_codon_change",
                                                                                                         "protein_coding:exon:CDS:deletion:inframe_variant:inframe_codon_loss"), "Shared coding change", "no")


---------------------------------------------------------------------------------------------------------------------------------

# Mutation burden stats:
get_mut_burden_stats_ancestral <- function(tree) {
  mut_burden = get_mut_burden_ancestral(tree)
  cat(paste("Mean mutation burden is", round(mean(mut_burden),digits = 1),"\n"))
  cat(paste("Range of mutation burden is", round(range(mut_burden)[1],digits = 1),"to",round(range(mut_burden)[2],digits = 1),"\n"))
  cat(paste("Standard deviation of mutation burden is", round(sd(mut_burden),digits = 1),"\n"))
  }

get_mut_burden_stats_ancestral(tree_SNV_c)

---------------------------------------------------------------------------------------------------------------------------------


save(tree_SNV_c, file=paste0(my_working_directory, "trees/tree"))
save(tree_SNV_c_ultra, file=paste0(my_working_directory, "trees/tree_ultra"))
save(details, file=paste0(my_working_directory, "trees/details")) 


---------------------------------------------------------------------------------------------------------------------------------




























