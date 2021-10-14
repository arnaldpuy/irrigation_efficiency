[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5551973.svg)](https://doi.org/10.5281/zenodo.5551973)

# Large variations in global irrigation withdrawals caused by uncertain irrigation efficiencies

Arnald Puy, Bruce Lankford, Jonas Meier, Saskia van der Kooij and Andrea Saltelli 

This is the R code of the paper, whose abstract is the following: 

*A proper assessment of the human impact on the global water cycle requires estimating the volume of water withdrawn for irrigation agriculture. A key parameter in this calculation is the irrigation efficiency, which represents the fraction of water lost to the crop due to bad management or conveyance losses. Here we use sensitivity auditing and uncertainty/sensitivity analysis to show that the irrigation efficiency values used in global hydrological modelling are spuriously accurate.  They downplay important ambiguities in partial efficiencies, irrigation technologies, the definition of "large-scale"" irrigated areas or managerial factors. If quantifiable uncertainties are accounted for, irrigation efficiencies turn from point-estimates to ranges spanning more than half the unit interval, making most countries largely undistinguishable. This uncertainty propagates to the computation of irrigation water withdrawals by making values vary by a factor of 3 and up to a factor of 20 at the country level. The irrigation module of global hydrological models should no longer overlook deep uncertainties if their estimates aim at informing irrigation and water-related policies in the real world.*

## Information
We provide the code in `.R`, `.rmd` and `.pdf` along with the system requirements. Contact Arnald Puy (arnald.puy@pm.me) for further information.

## Datasets included
We include all `.csv` datasets needed to replicate our work. The `.nc` files weight too much and should be downloaded from either ISIMIP (https://www.isimip.org).

`rohwer_data_all.csv`: produced from Rohwer et al. 2007.

`bos_data.csv`: produced from Bos and Nugteren 1990.

`solley_data.csv`: produced from Solley et al. 1998.

`fao_1997.csv`: produced from FAO 1997.

`ranges_efficiencies.csv`: produced from Clemmens and Molden 2007, Van Halsema and Vincent 2012, Rohwer et al. 2007, Rogers et al. 1997 and Brouwer et al. 1989.

`efficiency_10.csv`: The proportion of large-scale irrigated areas in each of the four irrigated area maps when the threshold used to identify large-scale irrigated areas is 10 km<sup>2</sup> (1,000 ha).

`efficiency_30.csv`: The proportion of large-scale irrigated areas in each of the four irrigated area maps when the threshold used to identify large-scale irrigated areas is 30 km<sup>2</sup> (3,000 ha).

`efficiency_100.csv`: The proportion of large-scale irrigated areas in each of the four irrigated area maps when the threshold used to identify large-scale irrigated areas is 100 km<sup>2</sup> (10,000 ha).

## Bibliography
M. Bos and J. Nugteren. On Irrigation Efficiencies. Tech. rep. Wageningen: International Institute for Land Recamation and Improvement / ILRI, 1990.

C. Brouwer, K. Prins, and M. Heibloem. Irrigation Water Management: Irrigation
Scheduling. Training Manual No 4. Tech. rep. Rome: FAO Land and Water Development
Division, 1989.

A. J. Clemmens and D. J. Molden. Water uses and productivity of irrigation systems".
Irrigation Science 25.3 (2007), 247:261. doi: 10.1007/s00271-007-0067-y

FAO. Irrigation potential in Africa. A basin approach. FAO Land and Water Bulletin
4. Tech. rep. Rome: Food and Agriculture Organization of the United Nations, 1997.

D. Rogers, F. Lamm, M. Alam, T. Trooien, G. C. P. Barnes, and K. Mankin. Efficiencies
and water losses of irrigation systems". Irrigation management series February
(1997), 1:6.

J. Rohwer, D. Gerten, and W. Lucht. Development of functional irrigation types for
improved global crop modelling". PIK Report 104 (2007), 1:61.

W. B. Solley, R. R. Pierce, and H. Perlman. Estimated Use of Water in the United
States in 1995. Tech. rep. US Geological Survey Circular 1200, 1998, 50.

G. E. Van Halsema and L. Vincent. Efficiency and productivity terms for water management: A matter of contextual relativism versus general absolutism". Agricultural Water Management 108 (2012), 9:15.





