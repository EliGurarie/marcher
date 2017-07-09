# marcher
Migration and Range Change Estimation in R (Development version)

> UPDATE 2017.04.12: A stable version of `marcher` is now on CRAN at https://cran.r-project.org/web/packages/marcher/index.html

The `marcher` package provides functions and tools for mechanistic range shift analysis decribed in Gurarie et al. (*in press*).  The methods are designed to estimate parameters of range shifting, including coordinates of the centroids of (2 or 3) ranges, the times of initiation and duration of each shift, ranging areas and time scales of location and velocity autocorrelation.  Because the estimates are likelihood based, there are several handy inferential tools including confidence intervals around all estimates and a sequence of hypothesis tests, including: (a.) What is the appropriate (highest) level of autocorrelation in the data? (b.) Is an estimated range shift significant? (c.) Is there a stop-over during the migration? (d.) Is a return migration a strict return migration? 
 
The vignette introduces the family of range shift models and illustrates methods to simulate, visualize, estimate and conduct the hypothesis tests. 


# References

Gurarie, E., Francesca, C.,  Peters, W., Fleming, C., Calabrese, J., Müller, T., & Fagan, W. (in press) A framework for modeling range shifts and migrations: Asking whether, whither, when, and will it return. *Journal of Animal Ecology*. 