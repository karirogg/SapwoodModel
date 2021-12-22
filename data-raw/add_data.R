source("R/sapwood_utils.R")

smaland <- read_excel(here("data-raw", "smaland.xlsx"),1)
varmland <- read_excel(here("data-raw", "varmland.xlsx"),1)
uppland <- read_excel(here("data-raw", "uppland.xlsx"),1)

usethis::use_data(smaland, overwrite = TRUE)
usethis::use_data(varmland, overwrite = TRUE)
usethis::use_data(uppland, overwrite = TRUE)
