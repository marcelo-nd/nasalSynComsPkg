# nasalSynComsPkg

Tools for analyzing nasal synthetic community (SynCom) sequencing and metabolomics data for the "Human Nasal Microbiome Synthetic Communities" paper.

This package provides a collection of functions to import, clean, align, analyze, and visualize
nasal SynCom datasets, including feature tables, metadata, clustering results, and downstream
statistical analyses.

## Installation

The recommended way to install `nasalSynComsPkg` is via **pak**, which automatically resolves
CRAN and Bioconductor dependencies:

```r
install.packages("pak")
pak::pak("github::marcelo-nd/nasalSynComsPkg")
```

This package depends on several Bioconductor packages (e.g. `limma`, `phyloseq`, `ComplexHeatmap`),
which are handled automatically by `pak`.

### System requirements
- **Windows**: Rtools
- **macOS**: Xcode Command Line Tools
- **Linux**: system compilers (usually already installed)

## Loading the package

```r
library(nasalSynComsPkg)
```

## Package scope

`nasalSynComsPkg` includes functions for:

- Importing and cleaning feature tables and metadata
- Filtering low-abundance or unwanted features
- Aligning metabolomics and sequencing data to sample metadata
- Clustering samples based on abundance profiles
- Computing distances and temporal dynamics
- Identifying cluster-specific markers using limma
- Generating publication-ready visualizations

## Examples

Many functions include examples in their documentation. To view them:

```r
?align_samples_attr
?cluster_samples
?filter_features_by_col_counts
```

Some examples are marked with `\dontrun{}` because they require real experimental data.

## Author

**Marcelo Navarro Diaz**/t
Postdoctoral researcher, University of Tübingen. Contact: marcelo.n.d@ciencias.unam.mx

## License

MIT License
