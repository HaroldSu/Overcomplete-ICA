# Overcomplete-ICA

This codes are based on R, to reproduce some figures in the paper "Overcomplete Independent Component Analysis via SDP" by Podosinnikova et al. (2019). 

The original paper and some other references are included in `papers/`.

Some reproduced figures are included in `img/`.

## Warnings
Some R packages, including "pracma", "Matrix", "R.matlab" and "imager", would be installed to run these codes. Since some of the function need to call MATLAB server, MATLAB should be installed and set following the instruction of R package "R.matlab" reference manual (please check https://cran.r-project.org/web/packages/R.matlab/R.matlab.pdf).

## Fig.2
Running the codes of [Fig_2.R](Fig_2.R) can reproduce the first two figures of Fig.2. It may run for hours.

## Fig.3
Running the codes of [Fig_3.R](Fig_2.R) can reproduce the first two columns of figures in Fig.3. It may run for hours.

## Fig.4
Running the codes of [Fig_4.R](Fig_2.R) will transform each image in [CIFAR-10 training batch 1](data_batch_1.mat)into grayscale and form 7×7 patches to extract the X matrix. Then the corresponding mixing matrix would be estimated by several algorithms. Hwoever, due to the limitation of R, the mixing components will not be plot as Fig.4, i.e only the numeric results but not the figure itself can be output. Since this scripts need to call MATLAB server for data processing, it can take a very long time to run.

## Fig.7
Running the codes of [Fig_7.R](Fig_7.R) can reproduce the first two figures of Fig.7. It may run for hours.