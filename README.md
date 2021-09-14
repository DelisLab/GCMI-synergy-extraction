# GCMI-synergy-extraction


Prerequisites:<br />
Matlab version R2019 or older..<br/>
GCMI repository (https://github.com/robince/gcmi)<br />
Generalized Louvain repository (https://github.com/GenLouvain/GenLouvain)



The following scripts (Supporting functions folder) have been adapted directly from other github repositories:<br />
'opnmf.m' and 'NNDSVD.M' (https://github.com/asotiras/brainparts)<br />
'community_louvain.m','get_components.m' and 'threshold_absolute.m' (https://github.com/brainlife/BCT)<br />
'number_connected_components.m' and 'threshold_by_giant_component.m' (https://github.com/CarloNicolini/communityalg)<br/>


## Implementation
-Set GCMI_synergies, GCMI and genLouvain folders on Matlab path.<br />
-Type 'main' into command line and hit enter.<br />
-Import data in the following format: A Mat file containing a matrix of timepoints x trials as rows and EMG channels as columns.<br />
-Fill in the relevant information and tick the desired boxes.<br />
-Click the Start button.

## Output
Mat_thresholded: Matrix of GCMI values that have been thresholded using a modified percolation analysis.<br />
Mat_unthresholded: Matrix of GCMI values that have not been thresholded. <br />
Thresholds: A vector specifying the threshold value identified within each layer of the multiplex network. <br />
Opt_rank: The optimal model rank identified using the generalised community detection protocol. <br />
Q: The Q-statistic for maximal modularity for Opt_rank. <br />
S: A matrix indicating the hard cluster assignment produced by the community detection protocol. <br />
PNMF: A structure with the results of the PNMF dimensionality reduction using Opt_rank as the input parameter.



## Optional Input
Manual model rank selection: Instead of finding the optimal model rank, the input parameter for PNMF can be specified manually. To do so, firstly import the data then tick the checkbox and input the desired model rank in the enabled field.<br />
Manual threshold selection: Instead of finding network layer specific threshold values, the threshold can be set manually. To do so, firstly import the data and then tick the checkbox and fill in the desired threshold value that will be used across all network layers


##  References
1. IN PRESS <br />
2. Ince, R. A., Giordano, B. L., Kayser, C., Rousselet, G. A., Gross, J., & Schyns, P. G. (2017). A statistical framework for neuroimaging data analysis based on mutual information estimated via a gaussian copula. Human brain mapping, 38(3), 1541-1573. <br />
3. Yuan, Z., & Oja, E. (2005, June). Projective nonnegative matrix factorization for image compression and feature extraction. In Scandinavian Conference on Image Analysis (pp. 333-342). Springer, Berlin, Heidelberg.<br />
4. Bordier, C., Nicolini, C., & Bifone, A. (2017). Graph analysis and modularity of brain functional connectivity networks: searching for the optimal threshold. Frontiers in neuroscience, 11, 441.<br />
