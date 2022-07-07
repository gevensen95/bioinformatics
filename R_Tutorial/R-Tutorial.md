Introduction to R Tutorial
================
Garrett Evensen
5/1/2022

## My First R Code

This is a general tutorial for R if you have never used R before. Make
sure you install both R (<https://www.r-project.org>) and RStudio
(<https://www.rstudio.com/products/rstudio/download/>). When you
download R, you will be asked to pick a CRAN. Make sure to pick a CRAN
that is closest to you because this will speed up the time it takes to
download packages.

You do not need to open R; you will only need to use RStudio. When you
open RStudio, you will see three windows: Console on the left side,
Environment in the top right, and Files in the bottom right. In the
console, you can type in all the commands listed in this document. You
can also click File at the top of the screen and open an R script, which
will open a new window. You can paste code in the R script and run your
code. This is a good place to keep your code because you can save it and
come back to it later. This document was written in R Markdown. If you
want to learn more about it, check out this tutorial:
<https://rmarkdown.rstudio.com/lesson-1.html>

To run code in your R script, you can highlight the code you want, and
then click Run at the top of the R script. You can also have your cursor
anywhere in a line of code and click Run. If you don’t want to click Run
every time, then you can use ctrl + enter.

## Basic Mathematics in R

``` r
2 + 2 
```

    ## [1] 4

``` r
2 * 2
```

    ## [1] 4

``` r
2 + 2 / 2
```

    ## [1] 3

``` r
2^3
```

    ## [1] 8

``` r
pi
```

    ## [1] 3.141593

``` r
exp(2)
```

    ## [1] 7.389056

``` r
sin(89)
```

    ## [1] 0.8600694

\##Making Objects in R

We covered a bit of basic mathematics in R. Now let’s make objects.
Objects are used all the time in R, so it is important to know how to
make and use them. You can also use = instead of \<-

``` r
num.four <- 4
num.four
```

    ## [1] 4

That was a simple object, but frequently you will store a string of
numbers, or characters, into an object. This is called a vector. You
have to use c() and put the numbers between the parentheses, separated
by commas. Also, it is important to know the kind of object, or what
class an object it is. You can find this out using the class function

``` r
numbers <- c(1, 2, 3, 4)
numbers
```

    ## [1] 1 2 3 4

``` r
class(numbers)
```

    ## [1] "numeric"

Alternatively, you can use : to specify a string of integers in order.

``` r
numbers2 <- 1:4
numbers2
```

    ## [1] 1 2 3 4

You can do the same thing for words or letters. You’ll notice that the
class of letters is character, rather than integer, which means that the
object has letters, or words, in it.

``` r
letters <- c("A", "B", "C", "D")
letters
```

    ## [1] "A" "B" "C" "D"

``` r
class(letters)
```

    ## [1] "character"

Now let’s make a matrix in R. But how do I make a matrix? The function
for making a matrix is called matrix, and if you put a ? before matrix,
then you can find out how to use the function. You have to give it a
string of numbers, or characters, and specify the number of rows and
columns.

``` r
?matrix

my.mat <- matrix(c(1:9), nrow = 3, ncol = 3, byrow = TRUE)
my.mat
```

    ##      [,1] [,2] [,3]
    ## [1,]    1    2    3
    ## [2,]    4    5    6
    ## [3,]    7    8    9

``` r
class(my.mat)
```

    ## [1] "matrix" "array"

You’ll notice that the class function told us that my.mat is a matrix
and also an array. An array is a generic term encompasses N dimensional
data (e.g. a 4D array). Let’s say we wanted to change the order of the
matrix. We can either remake it using by column, or transform it.

``` r
my.mat2 <- matrix(c(1:9), nrow = 3, ncol = 3, byrow = FALSE)
my.mat2
```

    ##      [,1] [,2] [,3]
    ## [1,]    1    4    7
    ## [2,]    2    5    8
    ## [3,]    3    6    9

``` r
my.mat3 <- t(my.mat)
my.mat3
```

    ##      [,1] [,2] [,3]
    ## [1,]    1    4    7
    ## [2,]    2    5    8
    ## [3,]    3    6    9

We made matrices, but now let’s make a data frame. It’s essentially a
table with columns that can be multiple data types. Here we will make a
data frame with both numbers and letters. To

``` r
my.df <- data.frame(numbers = c(1, 2, 3, 4, 5),
                    letters = c('A', 'B', 'C', 'D', 'E'))
my.df
```

    ##   numbers letters
    ## 1       1       A
    ## 2       2       B
    ## 3       3       C
    ## 4       4       D
    ## 5       5       E

``` r
class(my.df)
```

    ## [1] "data.frame"

To see just one column of a data frame, add a $ at the end of the name
of the data frame, then add the column name. R will autocomplete the
name of the data frame and the column.

``` r
my.df$numbers
```

    ## [1] 1 2 3 4 5

``` r
class(my.df$numbers)
```

    ## [1] "numeric"

``` r
class(my.df$letters)
```

    ## [1] "character"

## Installing and Loading Packages

There is a lot you can do in what’s called base R, but you will have to
install and load R packages. We will start with what’s called the
tidyverse. There is A LOT you can do with the tidyverse, so we will only
cover some of the basics. I urge you to check out the tidyverse website
to learn more (<https://www.tidyverse.org>). The tidyverse is a
collection of R packages, so when you install and load the tidyverse,
you will install and load several packages.

To install a package, use the install.packages function. Here is how to
install the tidyverse: install.packages(“tidyverse”) R may ask if you
want to update other packages, say to update all of them.

It may take a minute or so to finish installing. Now that it’s
installed, let’s load it using the library function. You may get an
error code, which will be in red text. Read the error codes and you
might be able to figure out the problem. For example, if it says the
package maps is not installed, then install the maps package and re-load
tidyverse.

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

You will see eight packages are loaded and their version numbers. Below
that, you’ll see Conflicts. This isn’t a problem. R is telling you that
a function in one of the loaded packaged is overriding another function.
Here, the filter function in dplyr is used instead of the filter
function in the stats package. Essentially, when you use filter, you
will be using dplyr’s filter, rather than filter from stats.

Accompanying this R Markdown file, are several csv files. If you do not
have them, then feel free to email me at either <k.evensen001@umb.edu>
or <kgsevensen@gmail.com>. Make sure all the csv files are all in the
same folder, also known as a directory. In this example, I have a folder
called R on my Desktop that has all of my csv files. To access these
files, we have to set our working directory.

To figure out the precise command, click session at the top of RStudio,
then set working directory, then choose to directory. From there you can
navigate to the folder that has your csv files and click open. Note: on
Windows, the folder will appear empty. This is normal. In the console,
you will see a command using the function setwd. You can copy and paste
that into your R script so you don’t have to manually set the working
directory every time you run this script.

``` r
setwd("~/Desktop/Github/R_Tutorial")
```

Let’s load in the penguins.csv file using the read_csv command. This is
a function in the tidyverse that loads data in the form of a tibble,
which is similar to data frame and can be used in the same way. Base R
has a function called read.csv that loads data as a data frame.

The data we are using comes from the palmerpenguins R package
(<https://allisonhorst.github.io/palmerpenguins/>).

``` r
penguins <- read.csv(file = "penguins.csv")
```

If your run just penguins like we did previously with my.df, you will
see that only the first few lines of the data are shown. Also, when you
first load the data, R tells you the classes of each of the column
because we have loaded the data as a tibble.

If you want see all of the data, as if you were in Excel or Google
Sheets, then you can use the View command: View(penguins)

To better understand the structure of your data, use the str function.
This is also known as stirring your data.

``` r
str(penguins)
```

    ## 'data.frame':    344 obs. of  8 variables:
    ##  $ species          : chr  "Adelie" "Adelie" "Adelie" "Adelie" ...
    ##  $ island           : chr  "Torgersen" "Torgersen" "Torgersen" "Torgersen" ...
    ##  $ bill_length_mm   : num  39.1 39.5 40.3 NA 36.7 39.3 38.9 39.2 34.1 42 ...
    ##  $ bill_depth_mm    : num  18.7 17.4 18 NA 19.3 20.6 17.8 19.6 18.1 20.2 ...
    ##  $ flipper_length_mm: int  181 186 195 NA 193 190 181 195 193 190 ...
    ##  $ body_mass_g      : int  3750 3800 3250 NA 3450 3650 3625 4675 3475 4250 ...
    ##  $ sex              : chr  "male" "female" "female" NA ...
    ##  $ year             : int  2007 2007 2007 2007 2007 2007 2007 2007 2007 2007 ...

Let’s see what’s in the first row of the data frame, then the first
column, and then a specific cell.

``` r
penguins[1,]
```

    ##   species    island bill_length_mm bill_depth_mm flipper_length_mm body_mass_g
    ## 1  Adelie Torgersen           39.1          18.7               181        3750
    ##    sex year
    ## 1 male 2007

``` r
penguins[,1]
```

    ##   [1] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ##   [7] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ##  [13] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ##  [19] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ##  [25] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ##  [31] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ##  [37] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ##  [43] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ##  [49] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ##  [55] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ##  [61] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ##  [67] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ##  [73] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ##  [79] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ##  [85] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ##  [91] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ##  [97] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ## [103] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ## [109] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ## [115] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ## [121] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ## [127] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ## [133] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ## [139] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ## [145] "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"    "Adelie"   
    ## [151] "Adelie"    "Adelie"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [157] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [163] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [169] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [175] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [181] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [187] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [193] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [199] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [205] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [211] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [217] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [223] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [229] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [235] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [241] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [247] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [253] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [259] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [265] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [271] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"   
    ## [277] "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap"
    ## [283] "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap"
    ## [289] "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap"
    ## [295] "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap"
    ## [301] "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap"
    ## [307] "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap"
    ## [313] "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap"
    ## [319] "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap"
    ## [325] "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap"
    ## [331] "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap"
    ## [337] "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap" "Chinstrap"
    ## [343] "Chinstrap" "Chinstrap"

``` r
penguins[2,6]
```

    ## [1] 3800

Now that we have loaded the penguins data, let’s explore it a little
more. First, we will get the count of each species in the penguins data
frame.

``` r
count(penguins, species)
```

    ##     species   n
    ## 1    Adelie 152
    ## 2 Chinstrap  68
    ## 3    Gentoo 124

Tidyverse uses pipes (%\>%), which can make things easier and are more
intuitive when trying to use different functions and filtering your
data. We will get the count of each species using pipes. First you type
in the data frame, then use a pipe (%\>%), and then you can use the
function. To make things more readable, you can hit enter after typing
the pipe. You can run the code from either line, and the code will still
run. However, if you only highlight 1 line, then it won’t work.

Try just running the first line. If you do, you’ll notice in the console
that there is a + instead of a \>. This is R’s way of telling you that
you haven’t finished typing in all your code. You can fix this by either
typing in the rest of the code or hitting ctrl + c. 

``` r
penguins %>%
  count(species)
```

    ##     species   n
    ## 1    Adelie 152
    ## 2 Chinstrap  68
    ## 3    Gentoo 124

Now let’s calculate the calculate the mean bill length for each species.
You can use nested pipes.

``` r
penguins %>%
  group_by(species) %>%
  summarise(mean.bill.length = mean(bill_length_mm))
```

    ## # A tibble: 3 × 2
    ##   species   mean.bill.length
    ##   <chr>                <dbl>
    ## 1 Adelie                NA  
    ## 2 Chinstrap             48.8
    ## 3 Gentoo                NA

OH NO! There are NA’s for both Adelie and Gentoo. This is because there
are some NA values in the original penguins data frame. We can get rid
of them with the drop_na function, piping the penguins data frame, and
saving it as a new data frame. You can also specify which column to drop
NA’s from by typing the column number in the parentheses in drop_na.

``` r
penguins2 <- penguins %>%
  drop_na()
```

Let’s calculate the mean bill length for each species again.

``` r
penguins2 %>%
  group_by(species) %>%
  summarise(mean.bill.length = mean(bill_length_mm))
```

    ## # A tibble: 3 × 2
    ##   species   mean.bill.length
    ##   <chr>                <dbl>
    ## 1 Adelie                38.8
    ## 2 Chinstrap             48.8
    ## 3 Gentoo                47.6

That’s better.

## Statistial Analyses in R

Before we can do some basic statistics, we will split the penguins2 data
frame into three, one for each species.

``` r
adelie <- penguins2 %>%
  filter(species == "Adelie")
chinstrap <- penguins2 %>%
  filter(species == 'Chinstrap')
gentoo <- penguins2 %>%
  filter(species == "Gentoo")
```

Our first analysis, will be a t test to compare the mean bill length in
adelie penguins and gentoo penguins.

``` r
t.test(adelie$bill_length_mm, gentoo$bill_length_mm, alternative = 'two.sided')
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  adelie$bill_length_mm and gentoo$bill_length_mm
    ## t = -24.286, df = 233.51, p-value < 2.2e-16
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -9.453448 -8.034741
    ## sample estimates:
    ## mean of x mean of y 
    ##  38.82397  47.56807

Based on the results of the t test, we can conclude that the means are
different, with gentoo penguins having a larger bill length.

Next, we will make a linear model (think of a trend line from Excel)
using the gentoo penguins data and analyze it.

``` r
p.lm <- lm(bill_length_mm ~ bill_depth_mm, data = gentoo) 
p.lm 
```

    ## 
    ## Call:
    ## lm(formula = bill_length_mm ~ bill_depth_mm, data = gentoo)
    ## 
    ## Coefficients:
    ##   (Intercept)  bill_depth_mm  
    ##         16.67           2.06

The coefficients indicate the y-intercept and the slope of the line.

Using the plot function, you can use some basic analyses to determine
the quality of the model. We won’t go into how to use each plot here.

``` r
par(mfrow = c(2, 2))  
plot(p.lm)
```

![](R-Tutorial_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

Finally, we will determine if the slope of the line is different from 0
using the p-value and get the r-squared of the line.

``` r
summary(p.lm)
```

    ## 
    ## Call:
    ## lm(formula = bill_length_mm ~ bill_depth_mm, data = gentoo)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -7.914 -1.445  0.125  1.315  7.904 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    16.6702     3.3110   5.035 1.75e-06 ***
    ## bill_depth_mm   2.0603     0.2203   9.352 7.34e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.36 on 117 degrees of freedom
    ## Multiple R-squared:  0.4277, Adjusted R-squared:  0.4229 
    ## F-statistic: 87.45 on 1 and 117 DF,  p-value: 7.337e-16

The linear model has an r-squared of 0.4277 and the p-value is much less
than 0.05, indicating that the slope is different from 0.

## Graphing in R

Now that we have done a few things in R, let’s graph the data from
penguins2 as a scatterplot first using the base R plot function and then
the ggplot package (part of tidyverse). We will graph bill depth on the
x axis and bill length on the y axis.

``` r
plot(penguins2$bill_depth_mm, penguins2$bill_length_mm, #plot with base r
     xlab = 'Bill Depth (mm)', ylab = 'Bill Length (mm)')
```

![](R-Tutorial_files/figure-gfm/pressure-1.png)<!-- -->

That plot isn’t the most informative because we can’t distinguish which
species is which, and it’s not the prettiest. We can fix that using
ggplot. If you click Help at the top, then hover over Cheatsheets, you
can download many different cheatsheets for working in R, including
ggplot2.

In the code chunk below, aes stands for aesthetics, which ggplot uses
for x and y axes and several other parameters, such as color and shape
of points. Geom_point tells ggplot to plot the data as points.
Theme_classic cleans up the graph a bit to make it more aesthetically
pleasing. You can use labs to specify the title of the plot and the axes
titles. Finally, stat_smooth adds a trend line for each species. This is
known as a linear model.

``` r
p1 <- ggplot(data = penguins2, aes(x = bill_depth_mm, 
                             y = bill_length_mm, 
                             color = species)) +
  geom_point() + 
  theme_classic() + 
  labs(title = "Bill Depth vs. Bill Length",
       x = "Bill Depth (mm)",
       y = "Bill Length (mm)") +
  stat_smooth(method = "lm", formula = "y~x") 
p1
```

![](R-Tutorial_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

This is a relatively simple graph, and there are many other options for
ggplot, that I won’t get into here.

Next, we will make a box plot of bill length for each species.

``` r
p2 <- ggplot(penguins2, aes(species, bill_length_mm, color = species)) + 
  geom_boxplot() +
  theme_classic() +
  labs(y = "Bill Length (mm)", x = '')
p2
```

![](R-Tutorial_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

The next plot we will make is a histogram of the body masses for the
entire dataset.

``` r
p3 <- ggplot(data = penguins2, aes(x = body_mass_g)) +
  geom_histogram(bins = 10) +
  theme_classic() + 
  labs(title = "Body Mass Histogram",
       x = "Body Mass (g)", y = "Count") 
p3
```

![](R-Tutorial_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

The last plot we will make using ggplot is a bar chart for the number of
individuals for each species. We will use code that we used earlier.
Make sure to install the RColorBrewer package so we can use a specific
color palette. We will also specify text size, color, and if it’s bold.

``` r
counts <- penguins %>%
  count(species)

library(RColorBrewer)
hmcol<-brewer.pal(3,"Dark2")

p4 <- ggplot(data = counts, aes(x = species, y = n)) +
  geom_bar(stat = 'identity', aes(fill = species)) +
  theme_classic() +
  labs(title = "Number of Each Species",
       x = "",
       y = "Number of Individuals") +
  theme(title = element_text(colour = "black", face = 'bold', size = 10),
        axis.text.x = element_text(colour = "black", face = 'bold', angle=90),
        axis.text.y = element_text(colour = "black", face = 'bold'),
        axis.title.y = element_text(colour = "black", face = "bold"),
        legend.title = element_text(colour = "black", face = 'bold', size = 10),
        legend.text = element_text(colour = "black", face = 'bold')) +
  scale_fill_manual(values = hmcol, name = "Species")
p4
```

![](R-Tutorial_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

We made four graphs from the penguins data. It can be a bit clunky to
have four separate graphs, so let’s combine them so they are in one
figure. First install gridExtra.

``` r
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
grid.arrange(p1, p2, p3, p4, widths = c(2,2))
```

![](R-Tutorial_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

You’ve successfully made four different graphs and grouped them together
into one figure. Congrats! Feel free to try out different ways of making
graphs. For example, you can make a graph for just Adelie and add color
= sex within the aes() in ggplot. What happens?

## Making Tables in R

There are a lot of different ways you can make tables (i.e., Word,
Excel, R). Here we will use the ggpubr package, so make sure to install
it using install.packages() first.

``` r
library(ggpubr)
```

The first table we will make will have summary information for each
penguin species.

``` r
params <- colnames(penguins2[c(3:6)]) #pulls out the column names for each 
mean.data <- penguins2 %>% # Specify data frame
  group_by(species) %>% # Specify column you want to group_by
  summarise_at(vars(params),              # Specify columns
               list(name = mean))               # Specify function to apply to each column
```

    ## Note: Using an external vector in selections is ambiguous.
    ## ℹ Use `all_of(params)` instead of `params` to silence this message.
    ## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This message is displayed once per session.

``` r
colnames(mean.data) <- c('Species', 'Bill Length (mm)', 'Bill Depth (mm)',
                         'Flipper Length (mm)', 'Body Mass (g)')
```

Now that we have a summary dataset. We will make our first table.

``` r
mean_table <- ggtexttable(mean.data, rows = NULL)
mean_table
```

![](R-Tutorial_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

Let’s say we are only interested in sex-specific data for the Adelie
penguin.

``` r
sex.data <- penguins2 %>% # Specify data frame
  filter(species == 'Adelie') %>%
  group_by(sex) %>% # Specify column you want to group_by
  summarise_at(vars(params),              # Specify columns
               list(name = mean))               # Specify function to apply to each column
colnames(sex.data) <- c('Sex', 'Bill Length (mm)', 'Bill Depth (mm)',
                         'Flipper Length (mm)', 'Body Mass (g)')
sex.data[2:5] <- round(sex.data[2:5], digits = 0) # rounds to nearest integer
table <- ggtexttable(sex.data, rows = NULL)
table
```

![](R-Tutorial_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

## End of Tutorial

That’s it for this tutorial. There are plenty of R tutorials out there
to try out. Also, feel free to check out TidyTuesday, which posts
datasets weekly (<https://github.com/rfordatascience/tidytuesday>)
