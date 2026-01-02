# nasalSynComsPkg

Tools for analyzing nasal synthetic community (SynCom) sequencing and metabolomics data.

This package provides a collection of functions to import, clean, align, analyze, and visualize
nasal SynCom datasets, including feature tables, metadata, clustering results, and downstream
statistical analyses.

## Installation

The recommended way to install `nasalSynComsPkg` is via **pak**, which automatically resolves
CRAN and Bioconductor dependencies:

```r
install.packages("pak")
pak::pak("MarceloNavarroDiaz/nasalSynComsPkg")
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

The package is designed for **research workflows** and is actively developed.

## Examples

Many functions include examples in their documentation. To view them:

```r
?align_samples_attr
?cluster_samples
?filter_features_by_col_counts
```

Some examples are marked with `\dontrun{}` because they require real experimental data.

## Development status

This package is under active development.  
The API may evolve as analyses expand and new datasets are incorporated.

## Author

**Marcelo Navarro Diaz**  
Postdoctoral researcher, microbiology  
University of Tübingen

## License

MIT License
