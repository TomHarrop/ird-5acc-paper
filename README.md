**This repository contains the manuscript source for the following article**:

> Thomas W. R. Harrop, Otho Mantegazza, Ai My Luong, Kevin Béthune, Mathias Lorieux, Stefan Jouannic, and Hélène Adam. A Set of *AP2*-like Genes Is Associated with Inflorescence Branching and Architecture in Domesticated Rice.” Journal of Experimental Botany, 2019. [10.1093/jxb/erz340](https://doi.org/10.1093/jxb/erz340).


### Data availability

- Sequence data from this article have been deposited with the National Center for Biotechnology Information Sequence Read Archive (SRA) under accession [PRJNA518559](http://www.ncbi.nlm.nih.gov/bioproject/518559).

### Reproducibility

- The code we used to analyse the RNAseq data and panicle phenotype data is hosted at [tomharrop/5acc](https://github.com/tomharrop/5acc), and the code for qPCR analysis is at [othomantegazza/5acc-qpcr](https://github.com/othomantegazza/5acc-qpcr).
- We used [`snakemake`](https://snakemake.readthedocs.io/en/stable/) to arrange analysis steps into workflows and monitor dependencies, and [`Singularity`](https://sylabs.io/singularity/) to capture the computing environment.
- The final results and all intermediate steps can be reproduced from the raw data using a single command using `snakemake` and `Singularity` with the two analysis repos.

### Contact

- For questions, contact the corresponding authors or open an issue in the GitHub repository.
