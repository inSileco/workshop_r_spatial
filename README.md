# workshop_r_spatial

Workshop to learn how to turn R into a powerful SIG. 

## Install packages 

Clone or download the repository, then 

```R
install.packages("remotes")
remotes::install_deps()
```

## Render the presentation locally

Presentation build with [Quarto](https://quarto.org/) (see [Revealjs](https://quarto.org/docs/presentations/revealjs/)). 


```R
quarto preview index.Qmd
```