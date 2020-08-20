# omgee
This R package provides an implementation of the GEE method for estimating multinomial probabilties for clustered data, for independent and exchangeable correlation structure.

## How to cite
```
Landsman V., Landsman D., Li C.S., Bang H. Overdispersion models for correlated multinomial data: Applications to blinding assessment. Statistics in Medicine. 2019; 38: 4963– 4976. https://doi.org/10.1002/sim.8344
```

## Installation
There are a few methods available for installing the package. You need to have the `devtools` package installed to use these methods:
```r
> install.packages('devtools')
```

### Using GitHub
Use the `install_github` method from devtools to install directly from this repository:
```r
> devtools::install_github('https://github.com/Express50/omgee')
```

### Using pre-built tar file
Download the latest released .tar.gz file from the [Releases](https://github.com/Express50/omgee/releases) page into an empty folder. Then, navigate to the location of the file in your command line, and run the following from an R shell:
```r
> devtools::install('.', dependencies=TRUE)
```

## Building
You can also build the code yourself by cloning the repository and using `R CMD BUILD`. Make sure to install the `rootSolve` dependecy before installing `omgee`, by running `install.packages('rootSolve')`
```sh
> git clone https://github.com/Express50/omgee omgee
> R CMD build omgee
> R CMD INSTALL omgee_1.0.tar.gz
```
