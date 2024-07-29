# sptwombling
## An R-package for Bayesian spatiotemporal Wombling

![Maintainer](https://img.shields.io/badge/maintainer-arh926-blue)

An `R`-package for performing Bayesian spatiotemporal wombling. Further details can be found in the paper titled, "Bayesian Spatiotemporal Wombling". (arXiv link: https://arxiv.org/abs/2407.17804).

<table>
  <tr>
    <td> <img width="600" height="260" src="https://github.com/user-attachments/assets/e656e6e6-a6ff-4995-a490-fce50b547dfb"/> </td>
    <td> <img width="300" height="260" src="https://github.com/user-attachments/assets/2f1321de-0d28-44b6-9d7b-de352288e67e"/> </td>
  </tr>
</table>

<table>
  <tr>
    <td> <img width="300" height="260" src="https://github.com/user-attachments/assets/9a99542c-e6ec-413e-9324-d40fade26355"/> </td>
    <td> <img width="300" height="260" src="https://media.giphy.com/media/v1.Y2lkPTc5MGI3NjExbTllc3N4YmlkYzUzZzVhZm85OGpta2ZhbnRoYzF5ZmoyM3FxZjgwdCZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/WUhewkgbu4SyAVRolR/giphy-downsized.gif"/> </td>
    <td> <img width="300" height="260" src="https://media.giphy.com/media/v1.Y2lkPTc5MGI3NjExMXd2emNnbTJrdDl5NDhpOHdycDJzNDFraXVzMnc4ZHI4OXdxOXB2cCZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/GJUuNwQRWiJkMrKjEl/giphy-downsized.gif"/> </td>
  </tr>
</table>

## Steps to follow for analysis

1. Load the data containing (a) space-time co-ordinates (b) response (c) covariates
2. Fit a spatiotemporal Bayesian hierarchical model to the data
3. Perform spatiotemporal differential process analysis
4. Locate or annotate planar curves of interest
5. Perform surface wombling over triangulated surface

Refer to the vignettes folder for further details on the above steps to reproduce the analysis. All spatiotemporal plotting relies on `ggplot2` and `rgl`.
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
