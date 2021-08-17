source("R/sapwood_utils.R")

dat_TA <- read_excel(here("data-raw", "dat_TA.xlsx"),1) %>%
          rename(H = HW, S=SW) %>%
          select(H,S,W)
#dat_KP <- read_excel(here("data-raw", "dat_KP.xlsx"),1)
#dat_61 <- read_excel(here("data-raw", "dat_61.xls"),1)
#dat_94 <- read_excel(here("data-raw", "dat_94.xlsx"),1)

usethis::use_data(dat_TA, overwrite = TRUE)
#usethis::use_data(dat_KP, overwrite = TRUE)
#usethis::use_data(dat_61, overwrite = TRUE)
#usethis::use_data(dat_94, overwrite = TRUE)
