<h1 align="center"><i>PV</i>aluator</h1>


&nbsp;


## Presentation
R package for determining the optimal p-value using mathematical methods.

The methods (Euclidean distance, geometry of Bézier curves) aim to determine the area where the curvature of the curve, representing the number of values obtained as a function of the p-value, is most pronounced.

&nbsp;


## Installation
### Linux installation
Open R console or RStudio and use this command :

```
install.packages("/your_path/PValuator_1.0.0.tar.gz", repo = NULL, type = "source")
```

### Windows installation


&nbsp;


## Package using

* Load the package : 

```
library(PValuator)
```

* Prepare a dataframe (df) containing two columns (value of each p-value tested and the associated number of occurrences / percentage)
* **NOTE** : if you have a table containing elements (e.g. genomic regions) associated with a p-value so you can prepare a dataframe usable by the package :

```
prepare_data(data, column, nbr.pv = 20)
```

### Method 1 : Bézier curves

* Display the graphical result for the method using the geometric properties of Bézier curves :

```
bezier_plot(df, cv = "no", sp = 0.5, fraction.var = 5e-3, main.title = "Automatic detection of p-value : Bezier curve", x.title = "-log10(p-value)", y.title = "Number / Percentage")
```

* Show the determined p-value :

```
bezier_pvalue(df, cv = "no", sp = 0.5, fraction.var = 5e-3)
```

### Method 2 : Euclidean distance

* Display the graphical result for the method using the Euclidean distance between a line and the curve :

```
euclidist_plot(df, cv = "no", sp = 0.5, fraction.var = 5e-3, main.title = "Detect the furthest point of the line", x.title = "-log10(p-value)", y.title = "Number / Percentage")
```

* Show the determined p-value :

```
euclidist_pvalue(df, cv = "no", sp = 0.5, fraction.var = 5e-3)
```


&nbsp;


<div align="center"><img src="All_logos.png" width="500"></div>
