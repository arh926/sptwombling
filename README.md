# sptwombling: An R-package for Bayesian spatiotemporal Wombling

![Maintainer](https://img.shields.io/badge/maintainer-arh926-blue)

Reference to the paper titled, "Bayesian spatiotemporal Wombling". (arXiv link: here).


<p align="center">
  <img width="600" height="550" src=https://github.com/user-attachments/assets/9a99542c-e6ec-413e-9324-d40fade26355>
  <img width="600" height="550" src=https://github.com/user-attachments/assets/e06d2221-d72d-4d23-887f-faa6cfafb713>
  <img width="600" height="550" src=https://github.com/user-attachments/assets/57403b92-aac2-4935-9dd2-f66a56fda072>
<p>


Illustration

Code for performing Bayesian spatiotemporal wombling

## Contents

1. Load the data containing (a) space-time co-ordinates (b) response (c) covariates
2. Fit a spatiotemporal Bayesian hierarchical model to the data
3. Perform spatotemporal differential process analysis
4. Locate or annotate planar curves of interest
5. Perform surface wombling over trangulated surface

Demonstration with synthetic data.

### Load the data and separate (a) co-ordinates (b) response (c) covariates
```
if(!require(devtools)) install.packages("devtools")
devtools::install_github('arh926/sptwombling')
require(sptwombling)
```

## Authors

| Name   | Email       |              |
|:------ |:----------- | :----------- |
| Aritra Halder (maintainer)| aritra.halder@drexel.edu   | Assistant Professor, Department of Biostatistics, Drexel University| 
| Didong Li (maintainer)| didongl@unc.edu   | Assistant Professor, Department of Biostatistics, University of North Carolina|
| Sudipto Banerjee | sudipto@ucla.edu   | Professor and Chair, Department of Biostatistics,  UCLA |
<!--- --->
