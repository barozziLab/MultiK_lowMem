## Minimal, memory conserving and parallelized implementation of MultiK

#### References
[MultiK publication](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02445-5) \
[original implementation](https://github.com/siyao-liu/MultiK/tree/main) \
[parallel implementation](https://github.com/ewilliams-uhn/MultiKParallel)

### General Information
This is a more flexible re-implementation of MultiK and does not require a seurat object as input, but directly runs on a matrix (cells x components) of any dimensionality reduction. In the seurat analysis workflow, these are PCA components, but you can also use NMF components (as for example identified by the DECIPHER-seq workflow).
In this implementation, leiden clustering is used, but it can easily be changed to other methods/algorithms by changing the main function. We are actively working on implementing this among other things as a parameter directly in the main function.
Note, that as a tradeoff for the memory-efficient behaviour of this implementation not all outputs from the original workflow can be produced. The code is under active development and we are working on adding as much of these features as possible. However, the diagnostic plots to choose the optimal k are available and those are the main aim of the method.
Below are two examples for using this MultiK implementation with a seurat object or directly on a matrix.

```
source("multiK_conserveMem.fct.R")
source("clustering.fct.R")
```

### Running MultiK on seurat objects

```
#~~~~~~~~~~~~~~~~~~~~~~
# Run on seurat objects
#######################
library(Seurat)

# load a dimensionality reduced seurat object
so <- readRDS("some_seurat_object.rds") # this is a processed version of the p3cl object, supplied by the original MultiK tutorial

# extract PCA embeddings and choose the number of PCs on which you want to cluster your cells
emb <- Embeddings(so, reduction = "pca")
emb <- emb[, 1:50]

# run multiK
result <- MultiK_par_conserveMem(input.matr = emb,
                         reps = 100, # decrease this number for a test run
                         pSample = 0.8,
                         numCores_first = 10,
                         numCores_second = 10,
                         numChunks = 4)

print(result$plots)
```

The parameters numCores_first and numCores_second let you choose the number of cores for parallelization of the 2 parts of the multiK workflow. The first part includes the clustering (which for large datasets is usually more RAM consuming, therefore a lower number might be useful for this part). The second part is the scoring of different solutions for k. This is quite memory conserving and should also run with the maximum number of available cores.
Note that, due to the implementation of the foreach (parallelization) package in R, the number of cores can usually not exceed ~120.
The numChunks parameter lets you choose into how many chunks the result matrix will be separated. The scoring will then be performed in $p = (numChunks * (numChunks + 1)) / 2$ pieces, therefore p should exceed the number of allocated cores.
However, for small datasets, using many cores for that part will not lead to a gain in speed.

### Running MultiK directly on a dimensionality reduced matrix

You can directly use the matrix as input for MultiK in that case.
```
#~~~~~~~~~~~~~~
# Run on matrix
###############

# load input matrix
your_input_matrix <- readRDS("some_input_matrix.rds") # these are the PCA embeddings of the above mentioned seurat object

# If you already have a dimensionality reduced dataset with a matrix of the format (cells x components), you can directly run multiK on that
result <- MultiK_par_conserveMem(input.matr = your_input_matrix,
                                 reps = 100, # decrease this number for a test run
                                 pSample = 0.8,
                                 numCores_first = 10,
                                 numCores_second = 10,
                                 numChunks = 4)

print(result$plots)
```
