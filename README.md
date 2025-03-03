
# LocusCCA.CVRtesting:  A Low-Rank, Sparse Canonical Correlation Analysis and Canonical Variate Regression Testing Framework for Brain Connectivity Analysis

`locusCCA.CVRtesting` is an R package implementing a two-stage statistical framework (Locus-CCA + CVR testing) designed for analyzing associations between brain connectivity and cognitive, behavioral, or clinical outcomes. The approach addresses challenges inherent in brain connectome analyses such as high dimensionality and noise, providing interpretable and powerful insights into neurodevelopmental and neuropsychiatric studies.


-   I. Installation
-   II. Method
-   III. Detailed Descriptions of the Functions
-   IV. A Toy Example

## I. Package Installation
You can easily install `locusCCA.CVRtesting` from GitHub with:

```r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("Emory-CBIS/locusCCA.CVRtesting")
library(locusCCA.CVRtesting)
```


## II. Method

###  Locus-CCA Method
Locus-CCA (**Low-rank Characterization of Brain Connectivity Matrices with Universal Sparsity by Canonical Correlation Analysis**) is a specialized form of Canonical Correlation Analysis designed explicitly for brain connectivity data. The method uses a low-rank decomposition combined with sparsity-inducing penalties, enabling it to identify robust, interpretable connectivity patterns associated with clinical or cognitive variables. 

Formally, Locus-CCA identifies $m$  canonical directions  𝐔 and 𝐕 by solving:

<img src="fig/objective.png" width="650" align="center"/>

where:

- **X** is an $n \times p$ stacked brain connectivity data matrix, where $p = V(V-1)/2$ represents the length of the vectorized upper triangle  of connectivity matrix of size $V\times V$.
- **Y** is an $n \times q$ matrix of clinical or cognitive variables.
- 𝐔 and 𝐕 are canonical direction weight matrices for connectivity and clinical variables, with dimensions $p \times m$ and $q \times m$ respectively.
- **𝐉ⱼ** and **Dⱼ** are low rank parameters, with dimensions $p \times R_j$ and $R_j \times R_j$ respectively. $R_j$ is latent rank of $j$ th canonical direction on brain connectivity.
- The function $\mathcal{L}(\cdot)$ extracts the vectorized upper-triangular portion from symmetric connectivity trait matrices.
- $\rho_1$ is the regularization parameters penalizing the sparisty of 𝐔. 

#### Method Highlights

- **Low-rank factorization**:  
  Locus-CCA employs low-rank decomposition to capture intrinsic, structured patterns within connectivity matrices. This efficiently reduces model complexity and enhances interpretability.

- **Universal Sparsity Regularization**:  
  Sparsity regularization (L1, Hardthreshold or SCAD penalties) is applied element-wise to canonical weights, ensuring robust and interpretable connectivity patterns that represent meaningful neural circuitry associated with clinical or behavioral phenotypes.

- **Canonical Correlation Maximization**:  
  Canonical directions (**XU**, **YV**) derived from Locus-CCA are optimized to achieve maximum canonical correlation, capturing the strongest possible linear relationships between brain connectivity and clinical/behavioral outcomes.

### CVR Testing Method

Canonical Variate Regression (CVR) Testing is a robust statistical procedure designed to evaluate the significance of canonical directions derived from Locus-CCA in predicting clinical or behavioral outcomes. Specifically, CVR testing assesses whether each identified canonical connectivity pattern (**XU**) significantly contributes to explaining variation in the univariate clinical or behavioral measures (**z**).

Formally, for each canonical direction  $j \leq m$, CVR testing calculates a test statistic $T_j$:


$$T_j = \frac{\sqrt{n}\mathbf{S}_j}{\sqrt{\mathbf{I}_j}} \xrightarrow{d} N(0,1)$$

*Comprehensive and complete methodological explanations of the CVR Testing framework are described in detail in our accompanying research paper (in preparation).*


<!--
\mathbf{S}_j = \frac{1}{n \hat{\sigma}^2}\left(\mathbf{z}-\hat{\mathbf{z}}\right)^\top\left(\mathbf{f}_j - \hat{\mathbf{f}}_j\right),\quad
\mathbf{I}_j = \frac{1}{n \hat{\sigma}^2}\left(\|\mathbf{f}_j\|^2 - \mathbf{f}_j^\top\hat{\mathbf{f}}_j\right),
$$

and

- $$\mathbf{z}$$ is the observed $$n$$-dimensional clinical or behavioral outcome vector.
- $$\hat{\mathbf{z}}$$ represents predictions from a Lasso regression of $$\mathbf{z}$$ onto the connectivity matrix $$\mathbf{X}$$.
- $$\hat{\sigma}^2$$ is the estimated residual variance from the Lasso regression.
- $f_j$ is the $j$ th canonical factor computed as the standardized projection of connectivity matrix onto the canonical direction:  $f_j = {XU_{.j}}$
- $$\hat{\mathbf{f}}_j$$ represents predictions of $$\mathbf{f}_j$$ obtained by regressing $$\mathbf{f}_j$$ onto the remaining connectivity features ($$\mathbf{X}$$) via a second Lasso regression to control for collinearity among predictors.
-->


#### Method Highlights

- **Statistical Significance Testing**:  
  The CVR testing procedure provides formal hypothesis tests for assessing the statistical significance of each canonical variate, determining whether connectivity patterns identified by Locus-CCA significantly predict clinical or behavioral outcomes.

- **Account for High Dimensionality and Edge Dependence**:  
  CVR transforms the high-dimensional testing challenge into testing an aggregated scalar statistic, effectively accounting for dependency structures among brain connectivity edges.

- **Reduced Multiple Testing Burden**:  
  By evaluating canonical variates rather than individual connectivity edges, CVR greatly reduces multiple comparison issues, enhancing statistical power and interpretability of findings.

By integrating Locus-CCA with CVR testing, this approach provides a comprehensive, statistically rigorous framework for linking complex brain connectivity data to meaningful clinical and behavioral phenotypes.


### Functions Overview
The structure of the package is as follows, and detailed descriptions of the function arguments are provided in the section below:

-   **Main Function:**
    -   `Locus_CCA`: performs Locus_CCA on brain connectivity and clinical/behavorial variables from the same group of subjects.
    -  `CVR_testing`: Peforms CVR testing based on the results from Locus_CCA.
-   **Tuning Parameter Selection:**
    -   `bic_cal()`: selects the tuning parameters $\rho_1$.
-   **Helper Functions:**
    -   `Ltrinv` and `Ltrans`: transform the brain connectivity to vectorized upper triangle and transform it back.
    - `plot_conn`:   both 
## III. Detailed Descriptions of the Functions

### 1. dynaLOCUS function

```         
dynaLOCUS(Y, q, V, n_subject, preprocess = TRUE, penalty_S = "L1", phi = 2, rho = 0.9, 
penalty_A = TRUE, lambda = 0.1, maxIteration = 100, speed_up = FALSE, espli1 = 0.01, espli2 = 0.05, silent = TRUE, demean = TRUE)
```

-   `Y`: Group-level dynamic connectivity data represented as a matrix of dimension $NT \times p$, where $NT$ denotes the number of subjects multiplied by the number of time points, and $p$ represents the number of edges in the connectivity network. To construct `Y` from dFC matrices, suppose we have $T$ dFC matrices from each of the $N$ subjects, where each matrix is a $V \times V$ symmetric matrix. We use the `Ltrans()` function to extract the upper triangular elements of each dFC matrix and convert them into a row vector of length $p = \frac{(V-1)V}{2}$. We then concatenate these vectors across time and subjects to obtain the group dFC data `Y`, which has dimensions $NT \times p$.
-   `q`: The number of connectivity traits to extract.
-   `V`: The number of nodes. Note that $p$ needs to be equal to $\frac{(V-1)V}{2}$.
-   `n_subject`: The total number of subjects in the study.
-   `preprocess`: If `TRUE`, the concatenated group-level connectivity data will be preprocessed. Defaults to `TRUE`.
-   `penalty_S`: The option for the penalization function for the sparsity regularization for the connectivity traits. Users can choose `"NULL"`, `"L1"`, or `"SCAD"` (smoothly clipped absolute deviation), introduced by [Fan and Li, 2001](https://www.jstor.org/stable/3085904). Defaults to `"L1"`.
-   `phi`: A tuning parameter for the element-wise penalty on `S`. Defaults to 2.
-   `rho`: A proportional tuning parameter ranging from 0 to 1 for the adaptive selection method to determine the number of ranks for modeling each connectivity trait. The value of `rho` represents the closeness of the connectivity traits estimated with and without the low-rank structure. A higher value of `rho` will lead to a higher rank. Defaults to 0.9.
-   `penalty_A`: The option for the temporal smoothness regularization for the trait loadings. If `TRUE`, a temporal smoothness penalty is included in the optimization function to encourage similarity in the trait loadings in adjacent time windows. Defaults to `TRUE`.
-   `lambda`: A numeric tuning parameter for the temporal smooth lasso penalty on `A`. Defaults to 0.1.
-   `maxIteration`: The maximum number of iterations. Defaults to 100.
-   `speed_up`: If `FALSE` (Default), use the Node-rotation algorithm (see Algorithm 1 in Appendix A of the paper) for learning dyna-LOCUS. If `TRUE`, use the Alternative algorithm (see Supplementary Materials Section 7) to further speed up the computation for learning dyna-LOCUS.
-   `espli1`: A number describing the tolerance for change on `A`.
-   `espli2`: A number describing the tolerance for change on `S`.
-   `demean`: If `TRUE`, demean each column of the group-level connectivity matrix `Y`. Defaults to `TRUE`.
-   `silent`: If `FALSE`, print out the penalty added on `A` and `S`. Defaults to `TRUE`.

The `dynaLOCUS` function serves as the primary function in the algorithm, implementing the novel decomposition method for brain network dynamic connectivity matrices using low-rank structure with uniform sparsity and temporal group lasso. Users can provide the group-level concatenated dynamic connectivity data as input, along with specifying parameters such as the number of connectivity traits to extract, the number of nodes under consideration, and the total number of subjects in the study, etc. The output of the `dynaLOCUS` function is a list comprising 3 components:

-   `A`: The mixing matrix ${a_{itl}}$ of dimension $NT \times q$.
-   `S`: The subnetworks extracted, with a dimension of $q \times p$, where $p$ is the number of edges. Each row in S represents a vectorized subnetwork, generated using the `Ltrans()` function in the package.
-   `theta`: A list of length $q$, where `theta[[i]]` contains the symmetric low-rank decomposition of $i$th subnetwork.

### 2. bic_selection function

```         
bic_selection(Y, q, V, n_subject, preprocess = TRUE, penalty = "L1", 
phi_grid_search = seq(1, 2, 0.2), rho_grid_search = c(0.9), maxIteration = 100, espli1 = 0.01, espli2 = 0.05, demean = TRUE, save_output = FALSE)
```

- `Y`: Group-level dynamic connectivity data is represented as a matrix of dimension $NT \times p$, where $NT$ denotes the number of subjects multiplied by the number of time points, and $p$ represents the number of edges in the connectivity network.
- `q`: The number of connectivity traits to extract.
- `V`: The number of nodes. Note that $p$ needs to be equal to $(V-1)V/2$.
- `n_subject`: The total number of subjects in the study.
- `preprocess`: If `TRUE`, the concatenated group-level connectivity data will be preprocessed. Defaults to `TRUE`.
- `penalty`: The option for the penalization function for the sparsity regularization for the connectivity traits. Users can choose "NULL", "L1", or "SCAD" (smoothly clipped absolute deviation), introduced by [Fan and Li, 2021](https://www.jstor.org/stable/3085904). Defaults to "L1".
- `phi_grid_search`: Grid search candidates for the tuning parameter $\phi$ for the penalty on connectivity traits. Defaults to `seq(1, 2, 0.2)`.
- `rho_grid_search`: Grid search candidates for the tuning parameter for the adaptive selection method to determine the number of ranks for modeling each connectivity trait. Defaults to `c(0.9)`.
- `maxIteration`: The maximum number of iterations. Defaults to `100`.
- `espli1`: A number describing the tolerance for change on A.
- `espli2`: A number describing the tolerance for change on S.
- `demean`: If `TRUE`, performs demeaning on the input data. Defaults to `TRUE`.
- `save_output`: If `TRUE`, saves the output. Defaults to `FALSE`.

The `bic_selection` function serves as a valuable guide for tuning the parameters $\phi$ and $\rho`. However, it is worth noting that in certain datasets, the choice may not be straightforward solely based on BIC. Tuning parameters can also be selected based on visual inspection of the extracted connectivity traits to achieve the desired level of sparsity and appealing neuroscience interpretation. The function outputs a list comprising two components:

- `bic_tab`: A dataframe containing BIC values per $\phi$ and $\rho$.
- `results`: A list of dyna-LOCUS output, if `save_output` is `TRUE`.

### 3. lambda_selection function

```         
lambda_selection(Y, q, V, n_subject, preprocess = TRUE, penalty_S = "L1", phi = 2, rho = 0.9, 
lambda_grid_search = c(NA, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1), maxIteration = 100, espli1 = 0.01, espli2 = 0.05, demean = TRUE, save_output = FALSE)
```

-   `Y`: Group-level dynamic connectivity data represented as a matrix of dimension $N \times p$, where $N$ denotes the number of subjects multiplied by the number of time points, and $p$ represents the number of edges in the connectivity network.
-   `q`: The number of connectivity traits to extract.
-   `V`: The number of nodes. Note that $p$ needs to be equal to $\frac{(V-1)V}{2}$.
-   `n_subject`: The total number of subjects in the study.
-   `preprocess`: If `TRUE`, the concatenated group-level connectivity data will be preprocessed. Defaults to `TRUE`.
-   `penalty_S`: The option for the penalization function for the sparsity regularization for the connectivity traits. Users can choose `"NULL"`, `"L1"`, or `"SCAD"`. Defaults to `"L1"`. 
-   `phi`: A tuning parameter for the element-wise penalty on `S`. Defaults to 2.
-   `rho`: A proportional tuning parameter ranging from 0 to 1 for the adaptive selection method to determine the number of ranks for modeling each connectivity trait. The value of `rho` represents the closeness of the connectivity traits estimated with and without the low-rank structure. A higher value of `rho` will lead to a higher rank. Defaults to 0.9.
-   `lambda_grid_search`: Grid search candidates of tuning parameter $\lambda$. Defaults to `c(NA, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1)`, where `NA` represents no temporal smoothness penalty.
-   `maxIteration`: The maximum number of iterations. Defaults to 100.
-   `espli1`: A number describing the tolerance for change on `A`.
-   `espli2`: A number describing the tolerance for change on `S`.
-   `demean`: If `TRUE`, performs demeaning on the input data. Defaults to `TRUE`.
-   `save_output`: If `TRUE`, saves the output. Defaults to `FALSE`.

The `lambda_selection` function serves as a valuable guide for tuning the parameter $\lambda$ with fixed values of $\phi$ and $\rho$. The function outputs a list comprising 3 components:

-   `selected_lambda`: The $\lambda$ value that achieves the minimal goodness-of-fit.
-   `goodness_of_fit`: Goodness-of-fit values for each $\lambda$ with no temporal smoothness penalty.
-   `results`: A list of dyna-LOCUS outputs, if `save_output` is set to `TRUE`.

### 4. reliability_index function

```         
reliability_index(Y, q_choices, V, n_subject, preprocess = TRUE, penalty_S = "L1", 
phi = 2, rho = 0.9, penalty_A = TRUE, lambda = 0.1, maxIteration = 100, espli1 = 0.01, espli2 = 0.05, silent = TRUE, demean = TRUE,seeds = 1:50)
```

-   `Y`: Group-level dynamic connectivity data is represented as a matrix of dimension $NT \times p$, where $NT$ denotes the number of subjects multiplied by the number of time points, and $p$ represents the number of edges in the connectivity network.
-   `q_choices`: A vector of integers representing different choices for the number of connectivity traits to extract.
-   `V`: The number of nodes. Note that $p$ needs to be equal to $(V-1)V/2$.
-   `n_subject`: The total number of subjects in the study.
-   `preprocess`: If `TRUE`, the concatenated group-level connectivity data will be preprocessed. Defaults to `TRUE`.
-   `penalty_S`: The option for the penalization function for the sparsity regularization for the connectivity traits. Users can choose `"NULL"`, `"L1"`, or `"SCAD"`. Defaults to `"L1"`. 
-   `phi`: A tuning parameter for the element-wise penalty on `S`. Defaults to 2.
-   `rho`: A proportional tuning parameter ranging from 0 to 1 for the adaptive selection method to determine the number of ranks for modeling each connectivity trait. The value of `rho` represents the closeness of the connectivity traits estimated with and without the low-rank structure. A higher value of `rho` will lead to a higher rank. Defaults to 0.9.
-   `penalty_A`: The option for the temporal smoothness regularization for the trait loadings. If `TRUE`, temporal smoothness penalty is included in the optimization function to encourage similarity in the trait loadings in adjacent time windows. Defaults to `TRUE`.
-   `lambda`: A numeric tuning parameter for the temporal smooth penalty on `A`. Defaults to 0.1.
-   `maxIteration`: The maximum number of iterations for the dyna-LOCUS algorithm.
-   `espli1`: A number describing the tolerance for change on `A`.
-   `espli2`: A number describing the tolerance for change on `S`.
-   `silent`: Logical, if `FALSE`, prints out the penalty added on `A` and `S`.
-   `demean`: Logical, if `TRUE`, demeans the data before processing.
-   `seeds`: A vector of integers to set seeds for bootstrap sampling, default is 1:50.

The `reliability_index` function calculates the reliability index and reproducibility (RI) for the dyna-LOCUS method based on input choices of $q$. The output of the `dynaLOCUS` function is a list comprising 2 components:

-   `reliability`: A list where each component corresponds to a different choice of $q$. Each component is of length $q$, representing the reliability index for each trait based on that choice of $q$.
-   `RI`: A list where each component corresponds to a different choice of $q$. Each component is of length $q$, representing the reproducibility (RI) for each trait based on that choice of $q$.

### 5. energy_var_cal function

```         
energy_var_cal(A, q, t_length)
```

-   `A`: The loading matrix.
-   `q`: The number of connectivity traits.
-   `t_length`: The length of the time series per subject.

The function calculates the energy and variation for connectivity traits. The output of the function is a list comprising 2 components:

-   `trait_energy`: A vector representing the mean log-transformed energy of each connectivity trait across subjects.
-   `trait_variation`: A vector representing the mean variation of each connectivity trait across subjects.

### 6. ccf_calculation function

```         
ccf_calculation(A, q, n_subject)
```

-   `A`: The trait loading matrix.
-   `q`: The number of connectivity traits.
-   `n_subject`: The number of subjects.

The function calculates cross-correlation function (CCF) for connectivity traits. The output of the function is an array:

-   `ccf_value_array`: An array containing the final CCF values for each pair of connectivity traits.

### 7. dFC_states function

```         
dfC_states(A, S, n_subject, q, n_examplers = 10, gap = 5, seed = 1, n_centers)
```

-   `A`: The loading matrix.
-   `S`: The dyna-LOCUS result connectivity traits.
-   `n_subject`: The total number of subjects in the study.
-   `q`: The number of connectivity traits.
-   `n_examplers`: The number of exemplars to select from each subject for clustering initialization. Default is 10.
-   `gap`: The gap between selected exemplars in time points. Default is 5.
-   `seed`: An integer for setting the random seed to ensure reproducibility. Default is 1.
-   `n_centers`: The number of centers (or clusters) to use in k-means clustering.

This function retrieves whole-brain dynamic functional connectivity (dFC) states based on the results from the dyna-LOCUS method. The output of the function is a matrix:

-   `Y_construct`: A matrix of reconstructed connectivity patterns based on the clustered loading matrix. Each row of the matrix represents the upper triangular part of a dFC state.

## IV. A Toy Example

In this section, we provide a toy example to demonstrate the implementation of the package. The toy example data contains dynamic functional connectivity (dFC) matrices from 20 subjects, each with 16 dFC matrices across time. Each dFC matrix is symmetric with dimensions of $V \times V$, where $V = 264$ is the number of nodes. These dFC matrices are generated based on connectivity traits extracted from real imaging data, using [Power's brain atlas](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3222858/). The input $Y$ matrix would be of dimension $NT \times p$, where $N = 20$ subjects, $T = 16$ time windows, and $p = V(V-1)/2$ edges. The first 16 rows of the $Y$ matrix contain the concatenated dFC matrices from subject 1, and the data in the following rows are arranged similarly for the subsequent subjects.

**How to prepare the input matrix Y from dFC matrices**: Suppose we have $T$ dFC matrices from each of the $N$ subjects, where each matrix is a $V \times V$ symmetric matrix. To generate our input matrix $Y$, we use the `Ltrans()` function to extract the upper triangular elements of each dFC matrix and convert them into a row vector of length $p = \frac{(V-1)V}{2}$. We then concatenate these vectors across time and subjects to obtain the group dFC data $Y$, which has dimensions $NT \times p$.

``` r
# library 
library(dynaLOCUS)

# import the toy example data from github
# note: the data size is 85.1 MB, and loading may take around 22 seconds
load(url("https://raw.githubusercontent.com/Emory-CBIS/dynaLOCUS_tutorial_data/main/Y.RData"))

# check the dimension
dim(Y)
```

We propose to select the number of connectivity traits $q$ to extract based on the reliability of the extracted traits evaluated via the `reliability_index()` function. Note that the function can be time-consuming as it requires two layers of for loops. In this toy example, we illustrate with 3 bootstrap samples for each choice of $q$. In real data analysis, please consider using at least 50 bootstrap samples (i.e., `seeds = 1:50`).

``` r
# calculate reliability
output = reliability_index(Y, q_choices = 2:4, V = 264, n_subject = 20, seeds = 1:3) 

# extract reliability
reliability = output$reliability

# calculate the average reliability among latent connectivity traits for each choice of q
lapply(reliability, mean)
```

The results are

```         
$`q = 2`
[1] 0.4370009

$`q = 3`
[1] 0.6598781

$`q = 4`
[1] 0.7318346
```

In this toy example, we consider possible values for $q$ of 2, 3, and 4. The reliability index increases considerably from $q = 2$ to $q = 3$, but only slightly when further increasing from $q = 3$ to $q = 4$. Furthermore, the reliability index for $q = 3$ already indicates substantial reliability. Therefore, we decide to select $q = 3$.

Next, we proceed to use the BIC-type criterion to select the hyperparameters $\rho$ and $\phi$. In this toy example, we fix $\rho$ at 0.9 and explore various values for $\phi$ to observe their impact on the BIC value. Our experiments with several real imaging datasets suggest the range of 0.5 to 3.5 works well for $\phi$. We recommend initially considering the range $seq(0.5, 3.5, 0.5)$ to evaluate the BIC.

``` r
# bic selection
bic_output = bic_selection(Y, q = 3, V = 264, n_subject = 20,
                           phi_grid_search = seq(0.1, 1, 0.2), rho_grid_search = 0.9) 
# retrive the results
bic_tab = bic_output$bic_tab

# visuliazation
library(ggplot2)
ggplot(bic_tab, aes(x = phi, y = log(bic_value))) +
  geom_line() +
  labs(x = "phi", y = "Log BIC") +
  theme_bw() +
  theme(text = element_text(size = 14))
```


We select $\phi = 0.5$ and $\rho = 0.9$ based on the BIC-type criterion. It is worth noting that the BIC criterion serves as a valuable guide in selecting the tuning the parameters $\phi$ and $\rho$. However, the choice may not always be straightforward solely based on BIC in practice. Therefore, besides the BIC criterion, users can also employ supplementary selection strategies, such as specifying tuning parameters based on the desired sparsity level and the neuroscience interpretations they aim to achieve in the extracted connectivity traits.

Next, we determine the value of $\lambda$ by evaluating the goodness of fit.

``` r
# lambda selection
goodness_of_fit = lambda_selection(Y, q = 3, V = 264, n_subject = 20, phi = 0.5, rho = 0.9,
                                   lambda_grid_search = c(NA, 1e-4, 1e-3, 1e-2, 1e-1, 5e-1))
goodness_of_fit$selected_lambda
# [1] 0.01
```

The results show that the selected $\lambda = 0.01$.

Next, we perform the dyna-LOCUS decomposition using the parameters we have just selected.

``` r
# dyna-LOCUS decomposition
result = dynaLOCUS(Y, q = 3, V = 264, n_subject = 20, 
                   phi = 0.5, rho = 0.9, penalty_A = TRUE, lambda = 0.01)
# [1] "Finished!"
```

We visualize the connectivity traits based on the Power's atlas. Please note that the visualization code is prepared based on the Power's atlas, and please modify as needed if other atlases are used. Please ensure that packages `ggplot2`, `gridExtra`, and `lattice` are installed, as they are required for plotting.

``` r
source(url("https://raw.githubusercontent.com/Emory-CBIS/dynaLOCUS_tutorial_data/main/visualization.R"))
plots = lapply(1:3, function(j) {
  levelplotX(result$S[j,], maint = paste0("Trait ", j))
})
# arrange the plots in a grid layout
g = do.call(arrangeGrob, c(plots, nrow = 1, ncol = 3))
# store the plot
ggsave(filename = "traits.png", g, width = 15, height = 5)
```

<img src="fig/traits.png" width="650" align="center"/>

For studies with large sample sizes and brain atlases involving a large number of nodes, users can use the alternative estimation algorithm to accelerate computation in learning dyna-LOCUS by specifying the `speed_up = TRUE` option. We present the code and result as follows:

``` r
# dyna-LOCUS decomposition with faster computation
result = dynaLOCUS(Y, q = 3, V = 264, n_subject = 20, speed_up = TRUE,
                   phi = 0.5, rho = 0.9, penalty_A = TRUE, lambda = 0.01)
# [1] "Finished!"
plots = lapply(1:3, function(j) {
  levelplotX(result$S[j,], maint = paste0("Trait ", j))
})
# arrange the plots in a grid layout
g = do.call(arrangeGrob, c(plots, nrow = 1, ncol = 3))
ggsave(filename = "traits_speed.png", g, width = 15, height = 5)
```


The functions have been thoroughly tested with PNC data. Due to the variability across different datasets, users are encouraged to try out both the `speed_up = TRUE` and `speed_up = FALSE` options to see which works good for their analysis. If any problems arise, please feel free to reach out for further assistance.

After extracting the connectivity traits, we can evaluate the reproducibility of the extracted traits.

``` r
# calculate reproducibility
output = reliability_index(Y, q_choices = 3, V = 264, n_subject = 20, seeds = 1:5, 
                           phi = 0.5, lambda = 0.01) 
# extract reproducibility
output$RI
# $`q = 3`
# [1] 0.6626315 0.9493696 0.9999014
```

For each connectivity trait, we can evaluate its energy and variation based on its temporal trait loadings:

``` r
# energy and variation
energy_var_cal(A = result$A, q = 3, t_length = 16)
# $trait_energy
# [1] -6.1659538 -0.8114545 -5.4254209
# 
# $trait_variation
# [1] 0.0008743673 0.0096093112 0.0010382041
```

Furthermore, we can use the cross-correlation function (CCF) function to evaluate the synchronization between traits for each subject:

```         
# ccf calculation
ccf_calculation(A = result$A, q = 3, n_subject = 20)
```
