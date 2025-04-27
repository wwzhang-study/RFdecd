
<!-- README.md is generated from README.Rmd. Please edit that file -->

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

# Install archived RefFreeEWAS package
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/RefFreeEWAS/RefFreeEWAS_2.2.tar.gz",
  repos = NULL,
  type = "source"
)
```

To install the deconf package successfully, you can manually download the .tar.gz file to your local machine and install it from there. Here's how you can do it:

(1) Visit the GitHub repository link:
https://github.com/wwzhang-study/RFdecd/raw/main/deps/deconf_1.0.1.tar.gz
(Right-click the link → Save Link As to download the file.)

(2) Save deconf_1.0.1.tar.gz to a known location on your computer (e.g., ~/Downloads/).

(3) Run the following R code, replacing PATH_TO_FILE with the actual path to your downloaded file:

```R
install.packages(
  "~/Downloads/deconf_1.0.1.tar.gz",
  repos = NULL,
  type = "source"
)
```
### Install and library RFdecd
```R
devtools::install_github(repo = "wwzhang-study/RFdecd",dependencies = TRUE,build_vignettes = TRUE,upgrade = "never")
library(RFdecd)
```
## 2. Workflow
The RFdecd algorithm follows a three-phase workflow to iteratively optimize feature selection and cell composition estimation. The diagram below illustrates the key steps:

![Figure 1: RFdecd Workflow](https://raw.githubusercontent.com/wwzhang-study/RFdecd/main/figures/Fig1.png)


**Step 1: Initialization**

1. **Feature Selection**: Select the top 1000 features with the largest coefficient of variation (CV) from the raw data matrix \( Y \), generating a reduced matrix \( Y_{M_0} \).

2. **Initial Deconvolution**: Perform reference-free deconvolution on \( Y_{M_0} \) to estimate initial cell-type profiles (\( W_1 \)) and proportions (\( H_1 \)).

3. **Error Calculation**: Compute the root mean squared error (RMSE[1]) between the reconstructed data (\( \widehat{Y} = W_1H_1 \)) and the original data \( Y \).

**Step 2: Iterative Optimization**

For each iteration \( i \) (1 ≤ \( i \) ≤ `TotalIter`):

4. **Feature Update**: Use six feature selection strategies (VAR, CV, SvC, DvC, PwD, RFdecd) to update the feature list \( M_i \) based on current proportions \( H_i \).
   
5. **Deconvolution**: Re-estimate profiles (\( W_{i+1} \)) and proportions (\( H_{i+1} \)) using the updated feature matrix \( Y_{M_i} \).

6. **Error Update**: Recalculate RMSE[i+1] for the new estimates.

**Step 3: Termination**

7. **Optimal Solution**: Select the proportion matrix \( H_{id} \) corresponding to the iteration with minimal RMSE as the final output.

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
Some explanations about the parameters:

- **DataType:** A covariate representing data type, either Gene expression or DNA methylation.

- **Y_raw:** A raw data matrix of complex samples. 

- **K:** The number of cell types.

- **CTSoption:** Feature selection options. 

(1) DEVarSelect_CV for selecting the top 1000 features with the largest coefficient of variation (i.e., CV) in the estimated cell-type profiles; 

(2) DEVarSelect_VAR for selecting the top 1000 features with the largest variation (i.e., VAR) in the estimated cell-type profiles; 

(3) DEVarSelect_1VSother for selecting the top 1000 features between one cell type and the other cell types (i.e., SvC);

(4) DEVarSelect_2VSother for selecting the top 1000 features between two cell types and the other cell types (i.e., DvC);

(5) DEVarSelect_pairwise for selecting the top 1000 features between one cell type and another (i.e., PwD); 

(6) DEVarSelect_RFdecd for selecting the top 1000 features between one cell type and the other cell types, as well as between two cell types and the other cell types (i.e., RFdecd).

- **nMarker:** The number of cell-type specific markers. 

- **InitMarker:** Initial markers.

- **TotalIter:** The number of iterations.

We can have the estimated cell-type proportion matrix:

```R
head(res_RFdecd$estProp)
```

## Citation
Coming soon.
