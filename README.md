<!--- Note: Enable LaTeX rendering by adding this to your browser:  
      https://github.com/orsharir/github-mathjax-extension -->  

# RFdecd: reference-free deconvolution of complex samples based on cross-cell type differential analysis

<!-- badges: start -->
<!-- badges: end -->

`RFdecd` is an R package for reference-free deconvolution based on feature selection.
It iteratively searches for cell-type-specific features by six feature selection options, 
and performs composition estimation.

## 1. Installation

### Install R Dependencies
```R
# Install development tools
if (!require("remotes")) install.packages("remotes")
if (!require("devtools")) install.packages("devtools")

# Install required dependencies
install.packages(c("quadprog","corpcor"))

# Install archived RefFreeEWAS package
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/RefFreeEWAS/RefFreeEWAS_2.2.tar.gz",
  repos = NULL,
  type = "source"
)
```

To install the **deconf** package successfully, follow these steps:

(1) Visit the GitHub repository link:
https://github.com/wwzhang-study/RFdecd/raw/main/deps/deconf_1.0.1.tar.gz
(Right-click the link → "Save Link As" to download **deconf_1.0.1.tar.gz**).

(2) Install from local file:

```R
install.packages(
  "PATH_TO_FILE/deconf_1.0.1.tar.gz", # Replace with your actual path
  repos = NULL,
  type = "source"
)
```
### Install RFdecd
```R
devtools::install_github(repo = "wwzhang-study/RFdecd",dependencies = TRUE,build_vignettes = TRUE,upgrade = "never")
library(RFdecd)
```
## 2. Workflow

<div align="center">
  <img src="https://raw.githubusercontent.com/wwzhang-study/RFdecd/main/figures/Fig1.png" 
       alt="RFdecd Workflow" 
       width="650"/>
</div>


<p align="center">Figure 1: RFdecd workflow</p>

The RFdecd workflow optimizes feature selection and composition estimation through three phases:

**Step 1: Initialization**

1a. **Feature Selection**: Select the top 1000 features with the largest coefficient of variation (CV) from the raw data matrix $Y$, generating a reduced matrix $Y_{M_0}$.

1b. **Initial Deconvolution**: Perform reference-free deconvolution on $Y_{M_0}$ to estimate initial cell-type profiles ($W_{1}$) and proportions ($H_{1}$).

1c. **Error Calculation**: Compute the root mean squared error (RMSE[1]) between the reconstructed data ($\hat{Y} = W_{1}H_{1}$) and the original data $Y$.

**Step 2: Iterative Optimization**

For each iteration $i$ (1 ≤ $i$ ≤ TotalIter):

2a. **Feature Update**: Use six feature selection strategies (VAR, CV, SvC, DvC, PwD, RFdecd) to update the feature list $M_{i}$ based on current proportions $H_{i}$.
   
2b. **Deconvolution**: Re-estimate profiles ($W_{i+1}$) and proportions ($H_{i+1}$) using the updated feature matrix $Y_{M_i}$.

2c. **Error Update**: Recalculate RMSE[i+1] for the new estimates.

**Step 3: Termination**

3a. **Optimal Solution**: Select the proportion matrix $H_{id}$ corresponding to the iteration with minimal RMSE as the final output.

## 3. Example
### Gene expression example
```R
N <- 100 # Simulate a gene expression dataset with 100 samples
K <- 4 # 4 cell types
P <- 5000 # 5000 features

## Simulate proportion matrix
Hmat <- matrix(runif(N*K), ncol = K)
Hmat <- sweep(Hmat, 1, rowSums(Hmat), FUN="/")

## Simulate reference matrix
Wmat <- matrix(abs(rnorm(P*K, 4, 4)), P, K)

## Simulate mixed expression profiles
Y_raw <- Wmat %*% t(Hmat) + abs(rnorm(N*K))
rownames(Y_raw) = paste0("gene",1:nrow(Y_raw))
colnames(Y_raw) <- paste0("Sample", 1:100)

# Run RFdecd for gene expression data
res_RFdecd <- RFdeconv(DataType = "Gene expression",
                       Y_raw,
                       K,
                       CTSoption = DEVarSelect_RFdecd,
                       nMarker = 1000,
                       InitMarker = NULL,
                       TotalIter = 30)
```

### DNA methylation example
```R
N <- 100 # Simulate a DNA methylation dataset with 100 samples
K <- 4 # 4 cell types
P <- 5000 # 5000 CpG sites

## Simulate proportion matrix
Hmat <- matrix(runif(N*K), ncol = K)
Hmat <- sweep(Hmat, 1, rowSums(Hmat), FUN="/")

## Simulate reference matrix
Wmat <- matrix(rbeta(P*K,2,5), P, K)

## Simulate mixed methylation data
Y_raw <- Wmat %*% t(Hmat) + rnorm(N*K,sd = 0.1)
rownames(Y_raw) = paste0("CpG",1:nrow(Y_raw))
colnames(Y_raw) <- paste0("Sample", 1:N)

# Run RFdecd for DNA methylation data
res_RFdecd <- RFdeconv(DataType = "DNA methylation",
                       Y_raw,
                       K,
                       CTSoption = DEVarSelect_RFdecd,
                       nMarker = 1000,
                       InitMarker = NULL,
                       TotalIter = 30)
```
Parameter specifications

- **DataType (required):** Data type specification, **"Gene expression" or "DNA methylation"**.

- **Y_raw:** A raw data matrix of complex samples. 

- **K:** The number of cell types.

- **CTSoption (required):** Feature selection options.

**Options:**

(1) DEVarSelect_CV -- Coefficient of variation; 

(2) DEVarSelect_VAR -- Variance; 

(3) DEVarSelect_1VSother -- Single vs. Composite (SvC);

(4) DEVarSelect_2VSother -- Dual vs. Composite (DvC);

(5) DEVarSelect_pairwise -- Pairwise Direct (PwD); 

(6) DEVarSelect_RFdecd -- RFdecd.

- **nMarker:** Number of markers to select (default: 1000). 

- **InitMarker:** Initial markers.

- **TotalIter:** Maximun iterations (default: 30).

We can have the estimated cell-type proportion matrix:

```R
head(res_RFdecd$estProp)
```

## Citation
Coming soon.
