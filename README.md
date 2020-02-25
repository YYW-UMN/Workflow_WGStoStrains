# Workflow_WGStoStrains
Reproducible workflow in R to process vcf files and characterize SNP patterns from WGS data

## Objective
Process raw variant calling files to characterize SNP patterns from WGS data 

## Data required:
- Variant calling files (genotypes, depth, minor allele read count, quality, etc)
- Meta data such as individual ID, state, collection date, sample type, etc.
- Reference genome in fasta format
- Gene annotation files

## Steps in the workflow:
Vcf files from running the WGS pipeline https://bitbucket.org/jgarbe/gopher-pipelines/src/default/ on https://www.msi.umn.edu.  

### 1. Creat BSgenome packages from reference genome
Follow this tutorial to create the R package for each reference genome
https://www.bioconductor.org/packages//2.7/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf
<p align="center">
<img src="https://github.com/YYW-UMN/Workflow_WGStoStrains/blob/master/Figures/BSgenome%20PackageExample.png" width="500" />
</p>

### 2. Process vcf files: 
  - filter out INDELs and regions prone to sequencing error 
  - filter out singleton SNPs
  - annotate and identify SNPs with functional consequences
  
### 3. Characterize SNP patterns:
  - Looking at all samples indepedently: this provides overall summaries at a single time point:
    - rare and common SNPs
    - SNPs shared among different geographical grouping variables such as State (figure below left)
    - SNPs presented in different type of samples, such as tissue, fecal, or blood samples (figure below right)
<p align="center">
<img src="https://github.com/YYW-UMN/Workflow_WGStoStrains/blob/master/Figures/VennDiagram.png" width="330" />
<img src="https://github.com/YYW-UMN/Workflow_WGStoStrains/blob/master/Figures/SampleType.png" width="400" />
</p>

  - Group samples by individuals: this provides temporal summaries for all time points:
    - SNP patterns within an individual over time (figure below left)
    - SNP patterns between herds or states over time (figure below right)
<p align="center">
<img src="https://github.com/YYW-UMN/Workflow_WGStoStrains/blob/master/Figures/Heatmap_samples_over_time.png" width="400" />
<img src="https://github.com/YYW-UMN/Workflow_WGStoStrains/blob/master/Figures/SNPovertime.png" width="400" />
</p>
    
### 4. Infer strains and ancestral relations: 
  - Phylogenetic tree 
  - Minimum spanning tree (figure below left)
<p align="center">
<img src="https://github.com/YYW-UMN/Workflow_WGStoStrains/blob/master/Figures/MixtureSamples.png" width="300" />
<img src="https://github.com/YYW-UMN/Workflow_WGStoStrains/blob/master/Figures/Strains_Clusters.png" width="300" />
</p>

  - Cluster sample to strains:
    - Identify mixture samples and their proportions 
    - Investigate how strains evolve over time (figure below)
<p align="center">
<img src="https://github.com/YYW-UMN/Workflow_WGStoStrains/blob/master/Figures/StrainsEvolve.png" width="700" />
</p>

