
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

- **X** is an $n \times p$ stacked brain connectivity data matrix, where $p = node(*node-1)/2$ represents the length of the vectorized upper triangle  of connectivity matrix of size $node \times node$.
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
    - `plot_conn`:  plots the canonical weights on brain connectivity in the form of heatmap for adjancency connectivity matrix.
## III. Detailed Descriptions of the Functions

### 1. Locus_CCA function

```         
Locus_CCA(X, Y, node, m, rho, gamma = 2.1,
                         penalt = 'L1', proportion = 0.9
                         , silent = FALSE, tol = 1e-3)
```

-   `X`: Group-level brain connectivity data represented as a matrix of dimension $n \times p$, where $n$ denotes the number of subjects, and $p$ represents the number of edges in the connectivity network. To construct `X` from connectivity matrices, suppose each matrix is a $node \times node$ symmetric matrix. We use the `Ltrans()` function to extract the upper triangular elements of each connectivity matrix and convert them into a row vector of length $p = (node-1)*node/{2}$. We then concatenate these vectors across  subjects to obtain the group connectivity data `X`, which has dimensions $n \times p$.
-   `Y`: Group-level clinical/behavioral subscale scores for the same subjects in `X`, the dimension is $n\times q$.
-   `node`: The number of nodes. Note that $p$ needs to be equal to $(node-1)*node/{2}$.
-   `m`: The number of canonical correlation components to extract.
-   `rho`: A tuning parameter for the element-wise penalty on 𝐔.
-   `gamma`: A tuning parameter only used if  SCAD penalty is used, default to be 2.1.
-   `penalt`: The option for the penalization function for the sparsity regularization for the canonical correlation directions on brain connectivity. Users can choose `"NULL"`,  `"Hardthreshold"`, `"L1"`, or `"SCAD"` (smoothly clipped absolute deviation), introduced by [Fan and Li, 2001](https://www.jstor.org/stable/3085904). Defaults to `"L1"`.
-   `proportion`: A proportional tuning parameter ranging from 0 to 1 to determine the number of ranks for modeling each connectivity trait. The value of `rho` represents the closeness of the connectivity traits estimated with and without the low-rank structure. A higher value of `rho` will lead to a higher rank. Defaults to 0.9.
-  `silent`: If `FALSE`, print out the training progress. Defaults to `FALSE`.
-   `tol`: A number describing the tolerance for change on parameters 𝐔 and 𝐕  .

The `Locus_CCA` function serves as the primary function in the algorithm, implementing the novel CCA method for investigating the association between brain network connectivity matrices and clinical/behavioral subscale scores using low-rank structure with uniform sparsity. Users can provide the group-level concatenated connectivity data and clinical/behavioral subscale as input, along with specifying parameters such as the number of connectivity traits to extract, the number of nodes under consideration, etc. The output of the `Locus_CCA` function is a list comprising 4 components:

-   `U`: The canonical correlation directions on brain connectivity with  dimension $p \times m$.
-   `V`: The canonical correlation directions on subscale scores with  dimension $q \times m$.
-   `CC`: A m by m matrix. Canonical correlations between the each corresponding projection, i.e the pairwise correlations between columns in $XU$ and $YV$.
-   `R`: A list of rank $R_j$, where `R[j]` contains the rank of the $i$th sub connectivity matrix.

### 2. CVR_testing function

```         
CVR_testing(U, X, z, lambda1 = NULL, lambda2 = NULL)
```

- `U`: The canonical correlation directions on brain connectivity with dimension $p \times m$, i.e., `U` from `Locus_CCA`.
- `X`: Group-level brain connectivity data represented as a matrix of dimension $n \times p$, where $n$ denotes the number of subjects, and $p$ represents the number of edges in the connectivity network.
- `z`: A numeric response vector (n x 1).
- `lambda1`: A numeric value for Lasso penalty in coefficients estimation. If `NULL`, it is determined using cross-validation.
- `lambda2`: A numeric value for constrained optimization in score calculation. If `NULL`, it is  determined by our procedure.

The `CVR_testing` function characterize the significance of each  brain connectivity canonical variant (**XU**) in explaining the one overall disorder or behavior score.  The function takes the arguments including canonical correlation directions on brain connectivity (estimated U from Locus_CCA), the group-level brain connectivity data **X**, and an univariate response (usually overall evaluation, such as ADHD overall score). The function outputs the testing statsitics of CVR testing of all m canonical components. 

- `T_stats`: A length m vector containing the test statistics for each canonical variants of brain connectivity. Each entry  follows a asymptotic normal distribution. 



### 3. lambda_selection function

```
BIC_cal(X, Y, U, V)        
```

- `X`: Group-level brain connectivity data represented as a matrix of dimension $n \times p$, where $n$ denotes the number of subjects, and $p$ represents the number of edges in the connectivity network.
- `Y`: Group-level clinical/behavioral subscale scores for the same subjects in `X`, the dimension is $n\times q$.
- `U`: The canonical correlation directions on brain connectivity with  dimension $p \times m$, or the outcome `U` from `Locus_CCA`.
- `V`: The canonical correlation directions on subscale scores with  dimension $q \times m$, or the outcome `V` from `Locus_CCA`.

`BIC_cal function serves as a valuable guide for tuning the parameters $\rho$.  The function outputs a single BIC value.  A model with lower BIC value is prefered. However, it is worth noting that in certain datasets, the choice may not be straightforward solely based on BIC. Tuning parameters can also be selected based on visual inspection of the extracted connectivity traits to achieve the desired level of sparsity and appealing neuroscience interpretation.



## IV. A Toy Example

In this section, we provide a toy example to demonstrate the implementation of the package. We generated toy example data **X**, **Y**, and **z** based on real estimated lantent connectivity traits from real brain connectivity and real clinical subscale dataset on cognition. 
Specifically, we generated connectivity matrices based on the real connectivity traits , using [Power's brain atlas](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3222858/). Each connectivity  trait is symmetric with dimensions of $V \times V$, where $V = 264$ is the number of nodes.   The input $X$ matrix would be of dimension $n \times p$, where $N = 100$ subjects and $p = V(V-1)/2$ edges. The first 16 rows of the $Y$ matrix contain the concatenated dFC matrices from subject 1, and the data in the following rows are arranged similarly for the subsequent subjects.

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
