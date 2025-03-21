Phylo location-scale plot trials
================
Szymon Drobniak
12/12/24

``` r
library(here)
library(tidyverse)
library(plot3D)
library(ggplot2)


data <- read.table(
    here("data", "Psittaciformes", "Psittaciformes_354spp.csv"),
    sep = ",", header = TRUE
)

data <- data %>%
    group_by(Phylo) %>%
    mutate(Genus = strsplit(Phylo, "_")[[1]][1]) %>%
    ungroup()

length(unique(data$Genus))
```

    [1] 84

``` r
data$beak_l_resid <- residuals(lm(log(Beak.Length_Culmen) ~ log(Mass), data = data))
data$beak_w_resid <- residuals(lm(log(Beak.Width) ~ log(Mass), data = data))

plot(beak_l_resid ~ log(Mass), data = data)
```

![](plot_experiments_files/figure-commonmark/unnamed-chunk-1-1.png)

``` r
plot(beak_l_resid ~ beak_w_resid, data = data)
```

![](plot_experiments_files/figure-commonmark/unnamed-chunk-1-2.png)

``` r
plot(log(Beak.Length_Culmen) ~ log(Beak.Width), data = data)
```

![](plot_experiments_files/figure-commonmark/unnamed-chunk-1-3.png)

``` r
ggplot(data = data, aes(
    x = log(Beak.Length_Culmen),
    y = log(Beak.Width), colour = Genus
)) +
    geom_point() +
    theme_bw() +
    stat_ellipse() +
    guides(colour = FALSE)
```

![](plot_experiments_files/figure-commonmark/unnamed-chunk-1-4.png)

``` r
ggplot(data = data, aes(
    x = beak_l_resid,
    y = beak_w_resid, colour = Genus
)) +
    geom_point() +
    theme_bw() +
    stat_ellipse() +
    guides(colour = FALSE)
```

![](plot_experiments_files/figure-commonmark/unnamed-chunk-1-5.png)

``` r
ngenera <- rev(sort(table(data$Genus)))
ngenera <- ngenera[ngenera > 8]

data <- data %>%
    group_by(Genus) %>%
    mutate(
        beak_l_v = var(beak_l_resid),
        beak_w_v = var(beak_w_resid),
        beak_v = beak_l_v + beak_w_v
    )

data1 <- data %>%
    filter(Genus %in% names(ngenera))

ggplot(data = data1, aes(
    x = beak_l_resid,
    y = beak_w_resid, colour = Genus
)) +
    geom_point() +
    theme_bw() +
    stat_ellipse() +
    guides(colour = FALSE)
```

![](plot_experiments_files/figure-commonmark/unnamed-chunk-1-6.png)

``` r
ggplot(data = data1, aes(
    x = beak_l_resid,
    y = beak_w_resid, colour = log(beak_v)
)) +
    geom_point() +
    theme_bw() +
    stat_ellipse(aes(group = Genus)) +
    guides(colour = FALSE) +
    scale_colour_gradient(low = "#0e91b9", high = "red")
```

![](plot_experiments_files/figure-commonmark/unnamed-chunk-1-7.png)
