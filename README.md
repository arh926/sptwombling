# sptwombling: An R-package for Bayesian spatiotemporal Wombling

![Maintainer](https://img.shields.io/badge/maintainer-arh926-blue)

Reference to the paper titled, "Bayesian Spatiotemporal Wombling". (arXiv link: here).

<table>
  <tr>
    <td> <img width="600" height="260" src="https://github.com/user-attachments/assets/15eb2dac-21ea-462c-a900-9131ca906fff"/> </td>
    <td> <img width="300" height="260" src="https://github.com/user-attachments/assets/2f1321de-0d28-44b6-9d7b-de352288e67e"/> </td>
  </tr>
</table>

<table>
  <tr>
    <td> <img width="300" height="260" src="https://github.com/user-attachments/assets/9a99542c-e6ec-413e-9324-d40fade26355"/> </td>
    <td> <img width="300" height="260" src="https://github.com/user-attachments/assets/e06d2221-d72d-4d23-887f-faa6cfafb713"/> </td>
    <td> <img width="300" height="260" src="https://github.com/user-attachments/assets/57403b92-aac2-4935-9dd2-f66a56fda072"/> </td>
  </tr>
</table>

## Steps to follow for analysis

1. Load the data containing (a) space-time co-ordinates (b) response (c) covariates
2. Fit a spatiotemporal Bayesian hierarchical model to the data
3. Perform spatiotemporal differential process analysis
4. Locate or annotate planar curves of interest
5. Perform surface wombling over triangulated surface

Refer to the vignettes folder for further details and steos to reproduce the analysis. All spatiotemporal plotting relies on `ggplot2` and `rgl`.
```
if("devtools" %in% rownames(installed.packages()) == FALSE) install.packages("devtools")
devtools::install_github('arh926/sptwombling')
require(sptwombling)
```

## Authors

| Name   | Email       |              |
|:------ |:----------- | :----------- |
| Aritra Halder (maintainer)| aritra.halder@drexel.edu   | Asst. Professor, Dept. of Biostatistics, Drexel University| 
| Didong Li | didongli@unc.edu   | Asst. Professor, Dept. of Biostatistics, University of North Carolina|
| Sudipto Banerjee | sudipto@ucla.edu   | Professor & Past Chair, Dept. of Biostatistics,  UCLA |
<!--- --->
