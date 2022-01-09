
# setting working directory
setwd("C:\\Users\\nerc-user\\OneDrive\\Documents\\GitHub\\HumDiv")

# setting seed so that I get the same result whenever I start with this seed
set.seed(12345)

# creating an "out" function cause it's useful
'%!in%' <- function(x,y)!('%in%'(x,y))

# formatting to html so that complex tables can be used
options(knitr.table.format = "html") 
