# Reduced varying coefficients of regional quantile for multiple responses (2025+)

## Overview
RVCRQ applies regional quantile regression to the reduced VC model for multiple responses, leveraging K-nearest neighbors fused Lasso to capture underlying clustered patterns effectively. By representing VC with a few principal component functions and their coefficients, the reduced structure of the model contributes to the improved effectiveness of the proposed method. A key advantage of this approach is its ability to provide insights into the varying relationships between multiple outcomes and risk factors across time and quantile levels. To facilitate this, we developed an ADMM algorithm to estimate the principal component functions and their coefficients.

## Main functions
- [demo_simulation.m](https://github.com/woorimjung/rvcrq/edit/main/demo_simulation.m)
 : Toy example running the proposed method.

- [RVC.m](https://github.com/woorimjung/rvcrq/edit/main/RVC.m)
 : Main ADMM algorithm solving optimization problem.

- [upd_D.m](https://github.com/woorimjung/rvcrq/edit/main/upd_D.m)
 : Embedded ADMM algorithm estimating the coefficient matrix for principal component functions.

- [BICforRVC.m](https://github.com/woorimjung/rvcrq/edit/main/BICforRVC.m)
 : Code to calculate the BIC of the proposed method.

- [MVC.m](https://github.com/woorimjung/rvcrq/edit/main/MVC.m)
 : Code for estimating varying coefficients using the KNN method for multiple responses.

- [BICforKNN.m](https://github.com/woorimjung/rvcrq/edit/main/BICforKNN.m)
 : Code to calculate the BIC of the KNN method.

- [supplementary_code](https://github.com/woorimjung/rvcrq/edit/main/supplementary_code)
 : Folder containing all the source files and functions associated with KNN fused Lasso

## Note
Due to the non-convex nature of the optimization objective, we utilized [VC_qt_knn_admm.m](https://github.com/woorimjung/rvcrq/edit/main/supplementary_code/VC_qt_knn_admm.m), which performs VC regional quantile regression via K-nearest neighbors fused Lasso for single response data, to obtain a suitable initial matrix for our approach. For more details, visit https://github.com/younghhk/software/tree/master/MATLAB.

We employ the parametric max-flow algorithm presented in the paper "On Total Variation Minimization and Surface Evolution Using Parametric Maximum Flows" by Antonin Chambolle and Jérôme Darbon (https://link.springer.com/article/10.1007/s11263-009-0238-9).


## Authors
[**Woorim Jung**](https://www.linkedin.com/in/우림-정-202875330)

M.S. Graduate, Department of Statistics, Sungkyunkwan University 
  
[**Eun Ryung Lee**](https://sites.google.com/view/eunryunglee/home)

Department of Statistics, Sungkyunkwan University 

[**Seyoung Park**](https://sites.google.com/view/seyoungpark/home)

Department of Applied Statistics, Yonsei University

[**Hyokyoung (Grace) Hong**](https://dceg.cancer.gov/about/staff-directory/hong-grace)
 
 Division of Cancer Epidemiology & Genetics, National Cancer Institute, NIH


* Please contact at [c7012evol@g.skku.edu] for any inquiries regarding the code.
