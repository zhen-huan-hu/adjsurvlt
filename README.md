# adjsurvlt

A SAS macro for estimating direct adjusted survival functions for time-to-event data with or without left truncation

### Author

Zhen-Huan Hu <novelkeny@gmail.com>

### Copyright

(C) 2019-2021 Zhen-Huan Hu

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.

SAS and all other SAS Institute Inc. product or service names are registered trademarks
or trademarks of SAS Institute Inc. in the USA and other countries.

## Specifications

### Function

`%adjsurvlt`

### Parameters

| Syntax | Description |
| :--- | :--- |
| `indata` | Specifies the input data set, which requires the following variables: a continuous variable for the observed time, an event indicator, a categorical variable indicating treatment groups, and covariates for the multi-variate model. For left-truncated data, the input data set also should contain the time to left truncation. Observations with missing values need to be excluded beforehand. |
| `time` | Specifies the observed time from the start of follow-up to the event of interest or right censoring, whichever happens first. |
| `event` | Specifies the event indicator variable. It takes `1` as event of interest and `0` as right censoring. |
| `strata` | Specifies the treatment group variable. |
| `covlst` | Specifies the covariate list. Multiple covariates can be separated by spaces. A `\|` sign can be used to separate continuous/ordinal variables from discrete categorical variables. For example, if variable `A` is categorical and variable `B` is continuous/ordinal, `covlst = A \| B` needs to be specified; if both are categorical, `covlst = A B` needs to be specified; and if both are continuous/ordinal, `covlst = \| A B` needs to be specified. |
| `ltime` | Specifies the left truncation time. The macro assumes no left truncation if it is left unspecified. |
| `seed` | Specifies the seed for random number generation. By default, `seed = 0` which makes the seed obtained from the Intel RdRand instruction if feasible, or the current time value. |
| `nsim` | Specifies the number of simulation processes for constructing confidence bands. By default, `nsim = 1000`. |
| `mintime` | Specifies the lower boundary of the time interval between which the confidence bands are constructed. By default, the macro takes the minimum observed event time for each treatment group. For confidence bands of the survival difference estimates, the macro uses the greater of the minimum observed event times between the treatment group pairs. If the specified value is less than the minimum observed event time, the macro falls back to the default. |
| `maxtime` | Specifies the upper boundary of the time interval between which the confidence bands are constructed. By default, the macro takes the maximum observed event time for each treatment group. For confidence bands of the survival difference estimates, the macro uses the greater of the maximum observed event times between the treatment group pairs. If the specified value exceeds the maximum observed event time, the macro falls back to the default. |
| `alpha` | Specifies the confidence level. By default, `alpha = 0.05`. |
| `timelist` | Specifies the time points at which to display the estimates. Setting it only affects the displayed results. |
| `outdata` | Specifies the name of the main output data set containing direct adjusted survival estimates. It also defines the prefixes of the names of three additional output data sets that contain summary information of the sample, survival difference estimates and simultaneous *p* values. The default value is `outdata = adjout`, which creates the following data sets: `adjout` which contains information associated with the adjusted survival estimates for individual treatment groups including event time, treatment groups, numbers at risk, survival estimates, standard error estimates, confidence limits and confidence bands; `adjoutinfo` which contains summary of the numbers of events and censored subjects as well as boundaries of the confidence bands; `adjoutdiff` which contains information associated with the survival difference estimates between treatment group pairs including event time, treatment group pairs, survival difference estimates, standard error estimates, confidence limits and confidence bands; `adjouttest` which contains information regarding the hypothesis testing of the difference between treatment group pairs including number of simulation processes, critical values, and simultaneous p values. The confidence limits and confidence bands for the adjusted survival estimates are calculated based on log-log transformation, while those for the survival difference estimates are based on linear transformation. |
| `outsurvplot` | Enables the macro to plot direct adjusted survival curves. When `outsurvplot = 1` is set, the macro outputs a plot file named `adjout-plot` that displays a grid of adjusted survival curves, and when `outsurvplot = 2` is set, a single plot with overlapped adjusted survival curves is produced instead. |
| `outdiffplot` | Enables the macro to plot direct adjusted survival difference curves. When `outdiffplot = 1` is specified, an image file named `adjoutdiff-plot` is created that displays a grid of survival difference curves. |
| `showci` | Toggles whether to display confidence limits in plots. The macro displays confidence limits by default. Setting `showci = 0` hides the confidence limits. |
| `showcb` | Toggles whether to display confidence bands in plots. The macro displays confidence bands by default. Setting `showcb = 0` hides the confidence bands. |
| `tickvalues` | Specifies values of the tick marks on the *x* axis. |
| `width` | Specifies width of the plots. By default, `width = 800px`. |
| `height` | Specifies height of the plots. By default, `height = 600px`. |
| `style` | Specifies style of the plots. By default, `style = statistical`. |
| `imagefmt` | Specifies image format of the plots. By default, `imagefmt = png`. |
| `noprint` | Hides output results. |
 
## Reference

Hu, ZH., Wang, HL., Gale, R.P. *et al*. A SAS macro for estimating direct adjusted survival functions for time-to-event data with or without left truncation. *Bone Marrow Transplant* (2021). <https://doi.org/10.1038/s41409-021-01435-2>
