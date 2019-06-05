# C-TALE-Normalization
Code and examples for C-TALE normalization

*All scripts are provided as is, without any warrenty and for use at your own risk. This is not the release of a software package. We are only providing this information and code in addition to a description of methods for making it easier to reproduce our analyses. We are not providing any support for these scripts.*
# Usage:
python2 C-TALE-normalize.py cool_file chr ROI_start ROI_end resolution coefficient output_cool_file


cool_file - Cool file with C-TALE HiC map.

chr, ROI_start, ROI_end - Genomics coordinates in bp of Region of interest. <chr> can be number or letter if you generate cool using hiclib.

resolution -resolution of C-TALE HiC map in bp.

coefficient - float that will be used for multiplicating around ROI (see original manuscript).

output_cool_file - name of new cool file that will be generated.

# IMPORTANT
Script generate new file with normalized matrix, but delete raw matrix. So for parsing new cool file you should use balance=False.
We fix this in future.
