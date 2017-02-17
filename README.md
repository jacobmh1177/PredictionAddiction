# Project Title: Predicting Drug Sensitivity from Genomic Factors

# Team Name: Prediction Addiction

# Members: Nicole Ferraro, Jacob Hoffman, David Nolan, Nicholas Rodriguez

# Overview:

* In the age of personalized medicine, identifying factors that can reveal information about a patient’s drug sensitivity is increasingly important as we increase our ability to customize treatment plans. Cancer is a disease area where personalized approaches are becoming increasingly common. The Cancer Cell Line Encyclopedia (https://portals.broadinstitute.org/ccle/home) has profiled 947 human cancer cell lines, providing gene expression data, copy number, and massively parallel sequencing to target mutations and indels in 1651 genes. Additionally, they have assessed the pharmacological profile for 24 anticancer drugs across 479 of these cancer cell lines. This data provides a rich opportunity to understand the molecular factors impacting a patient’s sensitivity to a particular treatment, and to infer how this impact changes depending on the type of cancer studied. Barretina, et al. has leveraged this data to train an elastic net regression model to predict sensitivity, as well as a naive Bayes classifier for discrete sensitivity classes [1]. We propose to replicate their results, and also attempt to improve upon their performance by implementing additional modeling strategies to understand which best represent the data.  We will assess the change in important features across both cell lines and drugs. Our results can be further applied to data from Iorio, et al., who also characterized drug sensitivity across a variety of cancer cell lines and drugs [2], and have made their data available (http://www.cancerrxgene.org/downloads), to assess the model robustness. Additionally, we hope to incorporate natural language processing techniques to mine the existing literature after we determine important feature-drug interactions to find existing work that supports our findings. Our preliminary approach will be to leverage pre-built word vectors found here (http://bio.nlplab.org/). These word vectors are built off of PubMed and PMC texts. We expect that while some of our results will support our the biological relationships revealed by the predictive models, we may also find some relationships outside the scope of these datasets.

# Specific Aims:
* Recapitulate the elastic net regression and Naive Bayes classifier in [1] and compare results, including important features and accuracy across cell lines, to those previously found. [Nicole]
* Implement different machine learning approaches to try and improve performance, such as random forests or neural net classifiers. [Nick]
* Validate results on an external dataset, and compare the results from this test dataset to those found in [2]. [David]
* Incorporate data mining of the existing literature via NLP to better elucidate the biological mechanism behind the connections between molecular features and drug sensitivity for particular drugs. [Jacob]

# Potential Setbacks:
* We may be unable to replicate the findings in the original Cancer Cell Line Encyclopedia paper, and it could be difficult to determine the exact cause of any discrepancies, whether they be genuine errors or differences in data processing.
* The run-time of training a variety of models may limit the amount of strategies we are able to implement, and our ability to optimize the code could be a limiting factor.
* Time management (as always…) 
* Limited experience with NLP tools.

# References:

1. Barretina, J., Caponigro, G., Stransky, N., Venkatesan, K., Margolin, A. A., Kim, S., ... & Reddy, A. (2012). The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. Nature, 483(7391), 603-607.
2. Iorio, F., Knijnenburg, T. A., Vis, D. J., Bignell, G. R., Menden, M. P., Schubert, M., ... & Cokelaer, T. (2016). A landscape of pharmacogenomic interactions in cancer. Cell, 166(3), 740-754.

