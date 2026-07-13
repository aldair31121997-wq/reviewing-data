# RNA structure prediction

## Transient (“trapping”) RNA structures

Putative transient RNA structures, also referred to as trapping structures, were proposed in the studies by Funk *et al.* and French *et al.* These structures are predicted to form through interactions between the 3′ and 5′ regions of the template RNA while it is being copied by the viral RNA-dependent RNA polymerase (RdRp).

To calculate these structures using the dataset from this study, first download and install the `slidingfold.py` tool developed by Mathis Funk, available at:

https://github.com/dr-funk/trapped-RdRp/tree/main/slidingfold

After installation, download the `allviruses.fasta` file provided in this repository and place it in the same directory as `slidingfold.py`.

Then run:

```bash
python3 slidingfold.py -t allviruses.fasta
```

This command generates a TSV file containing the predicted transient RNA structures that may occur during HA replication at each analysed position. The analyses were performed using the version of `slidingfold.py` available in the github repository (v1.0).

## Stationary RNA structures

Stationary RNA structures were predicted using the ViennaRNA RNAfold web server:

http://rna.tbi.univie.ac.at/

Select the **RNAfold Web Server** option and analyse individually the complete sequences contained in `sequences-stationarystructures.fasta`.

Each sequence corresponds to an 80-nt window analysed in this study. Folding predictions were performed using the default parameters available on the RNAfold web server.

## References

1. Funk, M., Spronken, M. I., Hutchinson, R. M., Arragain, B., Juyoux, P., Bestebroer, T. M., de Bruin, A. C. M., Gultyaev, A. P., Fouchier, R. A. M., Cusack, S., te Velthuis, A. J. W., and Richard, M. (2026). *Polymerase trapping as the mechanism of H5 highly pathogenic avian influenza virus genesis*. **Science, 391**(6790), eadr6632. doi: 10.1126/science.adr6632.

2. Lorenz, R., Bernhart, S. H., Höner zu Siederdissen, C., Tafer, H., Flamm, C., Stadler, P. F., and Hofacker, I. L. (2011). *ViennaRNA Package 2.0*. **Algorithms for Molecular Biology, 6**(1), 26. doi: 10.1186/1748-7188-6-26.



