# C-TALE-Normalization
Code and examples for C-TALE normalization

*All scripts are provided as is, without any warrenty and for use at your own risk. This is not the release of a software package. We are only providing this information and code in addition to a description of methods for making it easier to reproduce our analyses. We are not providing any support for these scripts.*

# Installation

Install using pip, e.g.
`pip install https://github.com/ArtemLuzhin/C-TALE-Normalization/archive/master.zip`

# Usage:

Minimal command ot run your C-TALE normalization:
ctale_normalize cool_file ROI_coordinates output_cool_file

cool_file - Cool file with C-TALE HiC map.

ROI_coordinates - Genomics coordinates in bp for Region of Interest in UCSC format.
Chromosome name should match the chromosome name in the cool_file.
For example, chr1:150,000-151,000;chr2:320000-420000 - if capturing multiple regions (not more than 1 per chromosome!), separate them with a semicolon.

output_cool_file - name of new cool file that will be generated.

Additional arguments you can add:

--mult_factor - float that will be used for multiplicating around ROI (see original manuscript, default 1.54).

--IC_steps - maximum number of iterative correction steps that can be applied (default 20)

--tolerance - float for required variance in coverage required for iterative correction convergence (default 10^-5)

## To cite, please use:
*Golov AK, Ulianov SV, Luzhin AV, Kalabusheva EP, Kantidze OL, Flyamer IM,
Razin SV, Gavrilov AA. C-TALE, a new cost-effective method for targeted
enrichment of Hi-C/3C-seq libraries. Methods. 2019 Jun 26. pii:
S1046-2023(18)30478-X. doi: 10.1016/j.ymeth.2019.06.022.*
