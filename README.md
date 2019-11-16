---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# HCAmtxzip

<!-- badges: start -->
<!-- badges: end -->

This package is not mature.

This package implements tools to access and manipulate Human Cell
Atlas '.mtx.zip' archives describing single cell expression
studies. The focus is on the pre-computed matrices available by
navigating projects on the Human Cell Atlas Data Portal [project
summary page][hca_projects].

[hca_projects]: https://data.humancellatlas.org/explore/projects

## Installation

To install this package and it's dependencies, use BiocManager


```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cran.r-project.org")

BiocManager::install("mtmorgan/HCAmtxzip")
```

Load the package and discover project titles for which pre-computed
mtx.zip archives are available


```r
suppressPackageStartupMessages({ library(HCAmtxzip) })
dd = discover()
dd
## # A tibble: 13 x 7
##    projectTitle      entryId    hits   fileFormat   size file    path      
##    <chr>             <chr>      <list> <chr>       <dbl> <chr>   <chr>     
##  1 A single-cell re… 4a95101c-… <name… mtx.zip    4.94e8 4a9510… https://s…
##  2 A single-cell tr… 8185730f-… <name… mtx.zip    1.24e8 818573… https://s…
##  3 Assessing the re… 005d611a-… <name… mtx.zip    3.23e8 005d61… https://s…
##  4 Census of Immune… cc95ff89-… <name… mtx.zip    2.39e9 cc95ff… https://s…
##  5 Dissecting the h… 4d6f6c96-… <name… mtx.zip    1.67e8 4d6f6c… https://s…
##  6 Ischaemic sensit… c4077b3c-… <name… mtx.zip    1.82e8 c4077b… https://s…
##  7 Profiling of CD3… 091cf39b-… <name… mtx.zip    2.54e8 091cf3… https://s…
##  8 Reconstructing t… f83165c5-… <name… mtx.zip    1.32e9 f83165… https://s…
##  9 Single cell prof… 116965f3-… <name… mtx.zip    6.49e7 116965… https://s…
## 10 Single cell tran… cddab57b-… <name… mtx.zip    6.24e7 cddab5… https://s…
## 11 Single-cell RNA-… 2043c65a-… <name… mtx.zip    7.67e7 2043c6… https://s…
## 12 Spatio-temporal … abe1a013-… <name… mtx.zip    5.22e8 abe1a0… https://s…
## 13 Structural Remod… f8aa201c-… <name… mtx.zip    1.41e8 f8aa20… https://s…
```

Find the smallest archive, and pull the path to the mtx.zip archive


```r
path <- dd %>%
    filter(row_number() == which.min(size)) %>%
    pull(path)
```

Download, unzip, and import the archive into _R_; the download is
cached (using [BiocFileCache][]) so only necessary once. The imported
data is represented as a [SingleCellExperiment][], a data structure at
the heart of _Bioconductor_ packages [orchestrating single cell
analysis][osca]. The `dim` output shows that there are 58347 features
assayed across 2544 cells.

[BiocFileCache]: https://bioconductor.org/packages/BiocFileCache
[SingleCellExperiment]: https://bioconductor.org/packages/SingleCellExperiment
[osca]: https://osca.bioconductor.org


```r
sce <- import.mtxzip(path, verbose = TRUE)
## download
## unzip
## rowData
## colData
## assays
## SingleCellExperiment
sce
## class: SingleCellExperiment 
## dim: 58347 2544 
## metadata(0):
## assays(1): counts
## rownames(58347): ENSG00000000003 ENSG00000000005 ...
##   ENSG00000284747 ENSG00000284748
## rowData names(8): V2 V3 ... V8 V9
## colnames(2544): 00ca0d37-b787-41a4-be59-2aff5b13b0bd
##   0103aed0-29c2-4b29-a02a-2b58036fe875 ...
##   fdb8ed17-e2f0-460a-bb25-9781d63eabf6
##   fe0d170e-af6e-4420-827b-27b125fec214
## colData names(43): genes_detected file_uuid ...
##   analysis_protocol.protocol_core.protocol_id
##   analysis_working_group_approval_status
## reducedDimNames(0):
## spikeNames(0):
## altExpNames(0):
```

The figure below provides a schematic view of the SingleCellExperiment
-- there is feature x cell matrix of counts, with row and column
annotations.

![sce][sce_image]

[sce_image]: https://raw.githubusercontent.com/Bioconductor/SummarizedExperiment/master/vignettes/SE.svg?sanitize=true

The object contains a matrix of feature x sample counts, and it's easy
to discover, e.g., that about 92\% the cells in the matrix are zeros.


```r
mean(assay(sce) == 0)
## [1] 0.9227145
```

Information about each feature can be extracted with `rowData(sce)`,
and about each cell with `colData(sce)`. Much of the `colData()` is
constant across all cells, e.g., they all come from the same
experiment. A helper function `colData_common()` summarizes these
common features, where we learn for instance that the experiment
involves the pancreas islet of Langerhans.


```r
colData_common(sce) %>% print(n=Inf)
## # A tibble: 27 x 2
##    key                                  value                              
##    <chr>                                <chr>                              
##  1 total_umis                           <NA>                               
##  2 emptydrops_is_cell                   <NA>                               
##  3 barcode                              <NA>                               
##  4 derived_organ_ontology               UBERON:0001264                     
##  5 derived_organ_label                  pancreas                           
##  6 derived_organ_parts_ontology         UBERON:0000006                     
##  7 derived_organ_parts_label            islet of Langerhans                
##  8 cell_suspension.genus_species.ontol… NCBITaxon:9606                     
##  9 cell_suspension.genus_species.ontol… Homo sapiens                       
## 10 donor_organism.diseases.ontology     PATO:0000461                       
## 11 donor_organism.diseases.ontology_la… normal                             
## 12 donor_organism.is_living             no                                 
## 13 specimen_from_organism.organ.ontolo… UBERON:0001264                     
## 14 specimen_from_organism.organ.ontolo… pancreas                           
## 15 specimen_from_organism.organ_parts.… UBERON:0000006                     
## 16 specimen_from_organism.organ_parts.… islet of Langerhans                
## 17 library_preparation_protocol.proven… 3ab6b486-f900-4f70-ab34-98859ac5f7…
## 18 library_preparation_protocol.input_… OBI:0000869                        
## 19 library_preparation_protocol.input_… polyA RNA extract                  
## 20 library_preparation_protocol.librar… EFO:0008931                        
## 21 library_preparation_protocol.librar… Smart-seq2                         
## 22 library_preparation_protocol.end_bi… full length                        
## 23 library_preparation_protocol.strand  unstranded                         
## 24 project.provenance.document_id       cddab57b-6868-4be4-806f-395ed9dd63…
## 25 project.project_core.project_short_… Single cell transcriptome analysis…
## 26 project.project_core.project_title   Single cell transcriptome analysis…
## 27 analysis_working_group_approval_sta… blessed
```

The helper function `colData_distinct()` extracts the columns that
differ between cells


```r
colData_distinct(sce)
## # A tibble: 2,544 x 16
##    genes_detected file_uuid file_version cell_suspension… specimen_from_o…
##             <int> <chr>     <chr>        <chr>            <chr>           
##  1           6924 6770c8ea… 2019-05-15T… 00ca0d37-b787-4… 9c1445a1-7287-4…
##  2           3171 68efbcbe… 2019-05-30T… 0103aed0-29c2-4… 14875995-58ca-4…
##  3           3838 209d615a… 2019-05-14T… 01a5dd09-db87-4… 56b6cd1e-7c2c-4…
##  4           4111 143a7005… 2019-05-14T… 020d39f9-9375-4… a1b35ebb-b79e-4…
##  5           5834 45e82c84… 2019-05-14T… 02583626-682b-4… 1f43dc7a-3f89-4…
##  6           2564 cb0b32b2… 2019-05-14T… 041637f8-d5c9-4… 14875995-58ca-4…
##  7           3152 2b5338c5… 2019-05-14T… 044472bd-588a-4… 67601c53-f042-4…
##  8           6375 66515cc6… 2019-05-14T… 046c1a85-77f7-4… 67601c53-f042-4…
##  9           4650 45e89c83… 2019-05-14T… 04f60cb7-5ced-4… 1f43dc7a-3f89-4…
## 10           4420 476dc6ce… 2019-05-14T… 061f92bf-fcfc-4… 56b6cd1e-7c2c-4…
## # … with 2,534 more rows, and 11 more variables:
## #   donor_organism.provenance.document_id <chr>,
## #   donor_organism.human_specific.ethnicity.ontology <chr>,
## #   donor_organism.human_specific.ethnicity.ontology_label <chr>,
## #   donor_organism.development_stage.ontology <chr>,
## #   donor_organism.development_stage.ontology_label <chr>,
## #   donor_organism.sex <chr>,
## #   analysis_protocol.provenance.document_id <chr>, dss_bundle_fqid <chr>,
## #   bundle_uuid <chr>, bundle_version <chr>,
## #   analysis_protocol.protocol_core.protocol_id <chr>
```

An additional helper function shortens column names to the shortest
'word' (defined by the occurence of `_` or `.`) suffix, as illustrated
with


```r
x <- setNames(nm = c("common_prefix.a_title", "common_prefix.a_name"))
names_abbreviate(x)
## [1] "title" "name"
names_abbreviate(x, map = TRUE)
## # A tibble: 2 x 2
##   name                  abbrev
##   <chr>                 <chr> 
## 1 common_prefix.a_title title 
## 2 common_prefix.a_name  name
```

This can make the column names a little easier to use in an
interactive session.


```r
colData_distinct(sce, abbreviate_colnames = TRUE)
## # A tibble: 2,544 x 16
##    detected file_uuid file_version suspension.prov… from_organism.p…
##       <int> <chr>     <chr>        <chr>            <chr>           
##  1     6924 6770c8ea… 2019-05-15T… 00ca0d37-b787-4… 9c1445a1-7287-4…
##  2     3171 68efbcbe… 2019-05-30T… 0103aed0-29c2-4… 14875995-58ca-4…
##  3     3838 209d615a… 2019-05-14T… 01a5dd09-db87-4… 56b6cd1e-7c2c-4…
##  4     4111 143a7005… 2019-05-14T… 020d39f9-9375-4… a1b35ebb-b79e-4…
##  5     5834 45e82c84… 2019-05-14T… 02583626-682b-4… 1f43dc7a-3f89-4…
##  6     2564 cb0b32b2… 2019-05-14T… 041637f8-d5c9-4… 14875995-58ca-4…
##  7     3152 2b5338c5… 2019-05-14T… 044472bd-588a-4… 67601c53-f042-4…
##  8     6375 66515cc6… 2019-05-14T… 046c1a85-77f7-4… 67601c53-f042-4…
##  9     4650 45e89c83… 2019-05-14T… 04f60cb7-5ced-4… 1f43dc7a-3f89-4…
## 10     4420 476dc6ce… 2019-05-14T… 061f92bf-fcfc-4… 56b6cd1e-7c2c-4…
## # … with 2,534 more rows, and 11 more variables:
## #   donor_organism.provenance.document_id <chr>, ethnicity.ontology <chr>,
## #   ethnicity.ontology_label <chr>, stage.ontology <chr>,
## #   stage.ontology_label <chr>, sex <chr>,
## #   protocol.provenance.document_id <chr>, fqid <chr>, bundle_uuid <chr>,
## #   bundle_version <chr>, protocol_id <chr>
```

Some exploration suggests that the column
`donor_organism.provenance.document_id` defines each biological
sample, so we see that the experiment consists of 8 individuals with
between 178 and 505 cells per individual.


```r
colData_distinct(sce, TRUE) %>%
    count(donor_organism.provenance.document_id)
## # A tibble: 8 x 2
##   donor_organism.provenance.document_id     n
##   <chr>                                 <int>
## 1 42e60811-4a08-45db-8db8-579f718f1675    331
## 2 6f1fd690-f44b-414e-ab7e-ea555b84942c    286
## 3 6fff3e7d-416e-4256-b33a-3448127e505e    474
## 4 a2675857-89d2-41a7-9178-f7c821cbc456    272
## 5 d361a5a0-19c0-4d5c-be21-117e93926292    505
## 6 e89af40b-6ef9-4b6a-8b6a-a51d1d72589c    221
## 7 f5b67f76-92f0-4426-aa6c-888b8865c403    277
## 8 fc0c4a2b-af93-42ec-8b68-10f68a1f622a    178
```

There are 6 male and 2 female samples with diverse ethnicity.


```r
colData_distinct(sce, TRUE) %>%
    count(donor_organism.provenance.document_id, sex, ethnicity.ontology_label)
## # A tibble: 8 x 4
##   donor_organism.provenance.docu… sex    ethnicity.ontology_label         n
##   <chr>                           <chr>  <chr>                        <int>
## 1 42e60811-4a08-45db-8db8-579f71… male   European                       331
## 2 6f1fd690-f44b-414e-ab7e-ea555b… male   Asian                          286
## 3 6fff3e7d-416e-4256-b33a-344812… female European                       474
## 4 a2675857-89d2-41a7-9178-f7c821… male   European                       272
## 5 d361a5a0-19c0-4d5c-be21-117e93… male   European                       505
## 6 e89af40b-6ef9-4b6a-8b6a-a51d1d… male   African American or Afro-Ca…   221
## 7 f5b67f76-92f0-4426-aa6c-888b88… female African American or Afro-Ca…   277
## 8 fc0c4a2b-af93-42ec-8b68-10f68a… male   ""                             178
```

The provenance `document_id` serves as a link to addition information
about the samples. The additional infomration can be extracted by
querying the HCA DSS. The information is returned as JSON, which is
relatively difficult to deal with in _R_; the even-less-mature
function `.files()` tries to parse this into a sensible format.


```r
uuids <- colData_distinct(sce, TRUE) %>%
    count(donor_organism.provenance.document_id) %>%
    pull(donor_organism.provenance.document_id)
donor <- HCAmtxzip:::.files(uuids)
donor
## # A tibble: 8 x 21
##   describedBy schema_type biomaterial_cor… biomaterial_cor… is_living sex  
##   <chr>       <chr>       <chr>            <chr>            <chr>     <chr>
## 1 https://sc… biomaterial DID_scRSq02      9606             no        male 
## 2 https://sc… biomaterial DID_scRSq05      9606             no        male 
## 3 https://sc… biomaterial DID_scRSq06      9606             no        fema…
## 4 https://sc… biomaterial DID_scRSq08      9606             no        male 
## 5 https://sc… biomaterial DID_scRSq04      9606             no        male 
## 6 https://sc… biomaterial DID_scRSq01      9606             no        male 
## 7 https://sc… biomaterial DID_scRSq07      9606             no        fema…
## 8 https://sc… biomaterial DID_scRSq03      9606             no        male 
## # … with 15 more variables: organism_age <chr>,
## #   organism_age_unit.text <chr>, organism_age_unit.ontology <chr>,
## #   organism_age_unit.ontology_label <chr>,
## #   human_specific.body_mass_index <chr>,
## #   human_specific.ethnicity.text <chr>,
## #   human_specific.ethnicity.ontology <chr>,
## #   human_specific.ethnicity.ontology_label <chr>,
## #   death.cause_of_death <chr>, development_stage.text <chr>,
## #   development_stage.ontology_label <chr>,
## #   development_stage.ontology <chr>, provenance.document_id <chr>,
## #   provenance.submission_date <chr>, provenance.update_date <chr>
```

We thus learn the sex, age and cause of death of each individual


```r
donor %>%
    select(provenance.document_id, sex, organism_age, death.cause_of_death)
## # A tibble: 8 x 4
##   provenance.document_id            sex    organism_age death.cause_of_dea…
##   <chr>                             <chr>  <chr>        <chr>              
## 1 42e60811-4a08-45db-8db8-579f718f… male   5            auto accident      
## 2 6f1fd690-f44b-414e-ab7e-ea555b84… male   22           head trauma        
## 3 6fff3e7d-416e-4256-b33a-3448127e… female 38           stroke             
## 4 a2675857-89d2-41a7-9178-f7c821cb… male   54           anoxia             
## 5 d361a5a0-19c0-4d5c-be21-117e9392… male   21           anoxia             
## 6 e89af40b-6ef9-4b6a-8b6a-a51d1d72… male   1            anoxia             
## 7 f5b67f76-92f0-4426-aa6c-888b8865… female 44           stroke             
## 8 fc0c4a2b-af93-42ec-8b68-10f68a1f… male   6            head trauma
```

Similarly detailed information about other aspects of the experiment,
such as contributors, funding, and publications, can be learned by
querying `.files()` with additional document identifiers from
elsewhere in `colData(sce)`. For instance, here's the (truncated, in
the HCA) description of the project.


```r
uuid <- colData_common(sce) %>% 
    filter(key %in% "project.provenance.document_id") %>% 
    pull(value)
project <- HCAmtxzip:::.files(uuid)
project
## # A tibble: 1 x 11
##   describedBy schema_type project_core.pr… project_core.pr…
##   <chr>       <chr>       <chr>            <chr>           
## 1 https://sc… project     Single cell tra… Single cell tra…
## # … with 7 more variables: project_core.project_description <chr>,
## #   supplementary_links <chr>, insdc_project_accessions <chr>,
## #   geo_series_accessions <chr>, provenance.document_id <chr>,
## #   provenance.submission_date <chr>, provenance.update_date <chr>
pull(project, "project_core.project_description") %>%
    strwrap(width = 80) %>%
    cat(sep="\n")
## As organisms age, cells accumulate genetic and epigenetic changes that
## eventually lead to impaired organ function or catastrophic failure such as
## cancer. Here we describe a single-cell transcriptome analysis of 2544 human
## pancreas cells from donors, spanning six decades of life. We find that islet
## cells from older donors have increased levels of disorder as measured both by
## noise in the transcriptome and by the number of cells which display
## inappropriate hormone expression, revealing a transcriptional instability
## associated with aging. By analyzing the spectrum of somatic mutations in single
## cells from previously-healthy donors, we find a specific age-dependent
## mutational signature characterized by C to A and C to G transversions,
## indicators of oxidative stress, which is absent in single cells from human
## brain tissue or in a tumor cell line. Cells carrying a high load of such
## mutations also express higher levels of stress and senescence markers,
## including FOS, JUN, and the cytoplasmic superoxide dismutase SOD1, markers
## previously linked to pancreatic diseases with substantial age-dependent risk,
## such as type 2 diabetes mellitus and adenocarcinoma. Thus, our single-cell
## approach unveils gene expression changes and somatic mutations acquired in
## aging human tissue, and identifies molecular pathways induced by these genetic
## changes that could influence human disease. Also, our results demonstrate the
## feasibility of using single-cell RNA-seq data from primary cells to derive
## meaningful insights into the genetic processes that operate on aging human
## tissue and to determine which molecular mechanisms are coordinated with these
## processes. Examination of single cells from primary human pancreas tissue
```

Here's information about the software in use to produce the above.


```r
sessionInfo()
## R Under development (unstable) (2019-11-05 r77368)
## Platform: x86_64-apple-darwin17.7.0 (64-bit)
## Running under: macOS High Sierra 10.13.6
## 
## Matrix products: default
## BLAS:   /Users/ma38727/bin/R-devel/lib/libRblas.dylib
## LAPACK: /Users/ma38727/bin/R-devel/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] HCAmtxzip_0.0.1             dplyr_0.8.3                
##  [3] SingleCellExperiment_1.9.0  SummarizedExperiment_1.17.0
##  [5] DelayedArray_0.13.0         BiocParallel_1.21.0        
##  [7] matrixStats_0.55.0          Biobase_2.47.0             
##  [9] GenomicRanges_1.39.1        GenomeInfoDb_1.23.0        
## [11] IRanges_2.21.1              S4Vectors_0.25.0           
## [13] BiocGenerics_0.33.0        
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_0.2.5       xfun_0.10              purrr_0.3.3           
##  [4] lattice_0.20-38        vctrs_0.2.0            BiocFileCache_1.11.2  
##  [7] utf8_1.1.4             blob_1.2.0             rlang_0.4.1           
## [10] pillar_1.4.2           glue_1.3.1             DBI_1.0.0             
## [13] rappdirs_0.3.1         bit64_0.9-7            dbplyr_1.4.2          
## [16] GenomeInfoDbData_1.2.2 stringr_1.4.0          zlibbioc_1.33.0       
## [19] evaluate_0.14          memoise_1.1.0          knitr_1.25            
## [22] curl_4.2               fansi_0.4.0            Rcpp_1.0.3            
## [25] backports_1.1.5        jsonlite_1.6           XVector_0.27.0        
## [28] bit_1.1-14             digest_0.6.22          stringi_1.4.3         
## [31] grid_4.0.0             cli_1.1.0              tools_4.0.0           
## [34] bitops_1.0-6           magrittr_1.5           RCurl_1.95-4.12       
## [37] tibble_2.1.3           RSQLite_2.1.2          crayon_1.3.4          
## [40] pkgconfig_2.0.3        zeallot_0.1.0          Matrix_1.2-17         
## [43] xml2_1.2.2             assertthat_0.2.1       httr_1.4.1            
## [46] R6_2.4.0               compiler_4.0.0
```
