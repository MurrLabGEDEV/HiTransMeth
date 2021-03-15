# HiTransMeth resources.

First run the generate_ref scripts to generate all reference sequences.
The data folder [here](/data/) contain the adapter sequences, barcode sequences, reference sequences and motif dataframe use.
Refer to those files for your construction.

### Dependencies

Dependencies needed for the script `generate_ref.py`:
* python > as 3.5.2
    * `pandas`
    * `biopython`

Dependencies needed for the main pipeline:
Pipeline use conda environment system to perform analysis with the requiered software version.
* [Miniconda3](https://docs.conda.io/en/latest/miniconda.html#linux-installers)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

To receive emails at start and end of the analysis you need:
* [mailx](http://www.linuxfromscratch.org/blfs/view/svn/basicnet/mailx.html)

---
### References

Refer to the article on [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/).

Data are available at GEO repository [GSE144524](https://www.ncbi.nlm.nih.gov/geo/submission/update/?acc=GSE144524).

For support with this analysis contact [Victor Ythier](victor.ythier@hcuge.ch)
