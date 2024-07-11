# Immune-ICC
Title: Multi-Omics-Profiling-Unveils-the-Complexity-and-Dynamics-of-Immune-Infiltrates-in-ICC

# Highlights
&nbsp;&nbsp;&nbsp;&nbsp;1. In the immune response mechanism of intrahepatic cholangiocarcinoma, the activation effect of antigen-presenting cells (such as B cells and macrophages) on T cells is enhanced, while the activation effect on Natural Killer cells is weakened.<br>
&nbsp;&nbsp;&nbsp;&nbsp;2. Infiltrating immune cells within the tumor microenvironment of intrahepatic cholangiocarcinoma exhibit varying immunoreactivity levels.<br>
&nbsp;&nbsp;&nbsp;&nbsp;3. The joint analysis of multi-omics data reveals the existence of four distinct infiltrating immune cell subpopulations, each possessing unique characteristics, within the tumor microenvironment of intrahepatic cholangiocarcinoma.<br>
&nbsp;&nbsp;&nbsp;&nbsp;4. The immune response model of intrahepatic cholangiocarcinoma has been revised, updated, and constructed.<br>
&nbsp;&nbsp;&nbsp;&nbsp;5. Novel methodology is implemented to examine alterations in immune cell activity, and a relatively innovative analytical procedure is devised.

# Innovation
I. From the perspective of bioinformatics
1) In the traditional classical immune response mechanism, when tumor cells stimulate APC (antigen presenting cells), the activated APC will have two effects: a) Activate NK cells to eliminate tumor cells; b) Activate T cells to eliminate tumor cells. However, in the globally public TCGA database (now also called the GDC database), the analysis of intrahepatic cholangiocarcinoma found that in the immune response mechanism, APC has a stronger activation effect on T cells, while its activation effect on NK cells is weakened (the proof method is in the first innovation of the second aspect), which has not been detailed in previous research findings.
2) APC mainly includes macrophages and B cells. From the GEO database, we borrowed the single-cell data of intrahepatic cholangiocarcinoma published by other teams (School of Life Sciences and Bioengineering, Beijing University of Technology and Cancer Prevention and Treatment Center, Sun Yat-sen University), and analyzed B cells, macrophages, NK cells, and T cells (in the second innovation of the second aspect). It is proposed that there are two subpopulations of B cells, one subpopulation of macrophages, and one subpopulation of T cells involved in the immune response mechanism of intrahepatic cholangiocarcinoma. And the impact of these four subpopulations on the survival rate of intrahepatic cholangiocarcinoma is described respectively.

II. From the perspective of analytical methods
1) According to the previous analytical methods, the cell proportion index is used to evaluate cell activity (such as Cibersort, xCell, etc.) for cell activity analysis. The larger the proportion, the stronger the activity is. However, under the guidance of clinical prior knowledge, we found that the cell proportion is easily affected by the cardinality (the denominator of the proportion), and it is mistakenly believed that the proportion decreases and its activity decreases. Because when a certain cell proliferates a certain number, its activity is essentially enhanced. However, because the number of other cells proliferates more, resulting in a greater increase in the cardinality, it will produce the illusion that the proportion result is smaller.  
&nbsp;&nbsp;&nbsp;&nbsp;Therefore, our analysis of cell activity is to observe the change in the slope of the regression line in the regression analysis. That is, the ssGSEA analysis method is used to calculate the enrichment score (which can be considered as the immune activity expression value) of each immune cell in the sample on the 16 immune-related gene sets we collected from TISCH2. Then, regression analysis is performed on the immune activity expression values ​​and proportions of B cells, macrophages, T cells, and NK cells. The slope of the regression line in the control group can be considered as the immune cell activity in a healthy state, and the slope of the regression line in the tumor group can be considered as the immune cell activity in a cancerous state. It can be observed that the slope of all regression lines is k>0. Because the immune activity expression value and proportion are positively correlated. The focus is to look at the slope k1 of the control group and the slope k2 of the tumor group, and observe the state of k1 changing to k2.  
2) We designed a different analysis process.  
&nbsp;&nbsp;&nbsp;&nbsp;a> The immune activity expression values ​​of the samples were clustered according to the sum of square differences to separate the high, medium and low immune activity groups. The tumor purity of the three groups was also verified, and the results were consistent;<br>
&nbsp;&nbsp;&nbsp;&nbsp;b> The high and low immune activity groups were differentially analyzed to screen out significant differential genes, and then WGCNA was used to further screen the differential genes to obtain four modular significant genes;<br>
&nbsp;&nbsp;&nbsp;&nbsp;c> The COX proportional hazard model was used to perform multivariate (modular significant genes) regression analysis on the samples of the high and low immune activity groups, and the survival ROC curve was used to predict the five-year survival time of the module genes, showing good results;<br>
&nbsp;&nbsp;&nbsp;&nbsp;d> PCA was used to analyze the module genes. The cumulative variance contribution of the first three principal components (PC1, PC2, PC3) exceeded 75%, and the high and low immune activity groups were shuffled. The first three principal components can be used to distinguish the samples well, with an error rate of 0. This verifies the previous screening;<br>
&nbsp;&nbsp;&nbsp;&nbsp;e> We used the Scissor tool to jointly analyze single-cell data, batch sample RNA-Seq data, and clinical survival data. Scissor has the ability to measure the Pearson correlation between single-cell and bulk samples to quantify the similarity between single-cell data and bulk sample data. Then, it is combined with clinical phenotypes to generate a correlation matrix. Finally, we set Scissor to select the COX regression model in the correlation matrix and apply sparse penalties and graph regularization to the regression model in order to select single cells with high confidence;<br>
&nbsp;&nbsp;&nbsp;&nbsp;f> At the end of the process, we selected four subsets of immune cells. And observed the significant genes corresponding to the subgroups and the performance of the enriched pathways.

# Data Availability Statement
&nbsp;&nbsp;&nbsp;&nbsp;The authors used publicly available datasets. 
The ICC transcriptome samples and their corresponding clinical phenotypes from GDC data portal (https://portal.gdc.cancer.gov/);<br>
UCSC Xena portal (https://xena.ucsc.edu/);<br>
Single cell samples from NCBI gene expression Omnibus (GEO) database (https://www.ncbi.nlm.nih.gov/; GSE138709, GSE159929);<br>
Gene set LM22 matrix obtained from CIBERSORT website (https://cibersortx.stanford.edu/);<br>
The collection of 16 immune-related gene sets contains 929 genes from TISCH2(http://tisch.comp-genomics.org/).

# Acknowledgments
&nbsp;&nbsp;&nbsp;&nbsp;The authors would like to thank Professor Ying Xu’s team from the School of Medicine of Southern University of Science and Technology and Professor Renchu Guan’s team from the School of Computer Science and Technology of Jilin University for their help and support in this research. They would also like to thank the organizations behind the online public database platforms mentioned in the study.

