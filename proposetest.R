library(magrittr)
library(ggplot2)

data(bears)
populations <- c("BLK", "PB",
                 "Bar", "Chi1", "Chi2", "Adm1", "Adm2",
                 "Denali", "Kenai", "Sweden") 

four <- bears %>% fit_all_graphs(c("BLK", "PB", "Denali", "Sweden"))
five <- bears %>% fit_all_graphs(c("BLK", "PB", "Kenai", "Adm1", "Sweden"))



summary(five)
summary(four)




plot(five)
plot(five, measure = "no_poor_fits")
