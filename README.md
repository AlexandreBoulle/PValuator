<h1 align="center"><i>PV</i>aluator</h1>


&nbsp;


## Presentation
R package for determining the optimal p-value using mathematical methods

&nbsp;


## Installation
### Linux installation
Open R console or RStudio and use this command :

```
install.packages("/ton_chemin/PValuator_1.0.0.tar.gz", repo = NULL, type = "source")
```

### Windows installation


&nbsp;


## Package using

* Load the package : 

```
library(PValuator)
```

* Prepare a dataframe containing two columns (value of each p-value tested and the associated number of occurrences / percentage)

### Method 1 : Bézier curves

* Display the graphical result for the method using the geometric properties of Bézier curves :

```
bezier_plot(df)
```

* Show the determined p-value :

```
bezier_pvalue(df)
```

### Method 2 : Euclidean distance

* Display the graphical result for the method using the Euclidean distance between a line and the curve :

```
euclidist_plot(df)
```

* Show the determined p-value :

```
euclidist_pvalue(df)
```


&nbsp;


<div align="center"><img src="All_logos.png" width="500"></div>
