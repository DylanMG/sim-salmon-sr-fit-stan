# fitting state-space spawner-recruitment models to simulated data with Stan

To create html summary knit the .Rmd file

Upcoming features:  
- improved visual representation of how the model indexes the raw ad transformed data/parameters  
- simulate and fit time-varying productivity (Ricker $\alpha$)  
- simulate and fit multiple stocks  
- simulate and fit multiple stocks with hierarchical $\alpha$


### A brief note on code style  
There are **a lot** of variables in this code; due in part to the nature of the more complex models, the amount of data required for these models, and gnarly indexing. 
We strive to use consistent naming conventions and practices as described in DFO's 2022 Technical Expertise in Stock Assessment (TESA) [Good Coding Practices workshop](https://cgrandin.github.io/good-code-practices/) (see [here](https://cgrandin.github.io/good-code-practices/presentation/Good-code-practices.html#36) for naming conventions). 
In short, we use lower *snake_case* for objects and functions, where names (or parts of names) may be capitalized if they are a matrix or acronym ([Edwards & Auger-Methe, 2018](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13105)), and characters after the first underscore typically qualify the first part of the name (e.g. `S_MSY` for spawners that produce maximum sustainable yield, `a_min` for the minimum age class, `R0` for recruits at time 0, etc.).
We attempt to be as descriptive as possible when introducing a new object or custom function, but will not spend time describing base and tidyverse functions; use R's built in help (i.e. typing `?function_name` into console).  
