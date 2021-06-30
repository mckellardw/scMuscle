# Downloads of the output files for the Visium database

Each .zip contains the outputs from the [spaceranger count](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/count) pipeline (from 10x Genomics), minus the .bam files to reduce download size.

Raw fastqs can be downloaded from GEO [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161318).

**Spatial RNA sequencing library preparation.** Tibialis anterior muscles of adult (5 mo) C57BL6/J mice were  injected with 10µl notexin (10 µg/ml) at 2, 5, and 7 days prior to collection. Upon collection, tibialis anterior muscles were isolated, embedded in OCT, and frozen fresh in liquid nitrogen. Spatially tagged cDNA libraries 365 were built using the Visium Spatial Gene Expression 3 Library Construction v1 Kit (10x Genomics, 366 Pleasanton, CA). Optimal tissue permeabilization time for 10 µm thick sections was found to be 15 minutes using the 10x Genomics Visium Tissue Optimization Kit. H&E stained tissue sections were imaged using Zeiss PALM MicroBeam laser capture microdissection system and the images were stitched and processed 369 using Fiji ImageJ software. cDNA libraries were sequenced on an Illumina NextSeq 500 using 150 cycle high output kits (Read 1=28bp, Read 2=120bp, Index 1=10bp, and Index 2=10bp). Frames around the capture area on the Visium slide were aligned manually and spots covering the tissue were selected using Loop Browser v4.0.0 software (10x Genomics). Sequencing data was then aligned to the mouse reference genome (mm10) using the spaceranger v1.0.0 pipeline to generate a feature-by-spot-barcode expression matrix (10x Genomics).

Additional info:
- Vis5A (2dpi), Vis7B (5dpi), Vis9A (7dpi)
- Visium slide number: V19T19-051
- Vis5A and Vis9A are left TA muscles, Vis7B is from a right TA
