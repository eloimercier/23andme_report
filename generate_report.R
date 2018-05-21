#install pandoc
#install miktex
#install.packages("kableExtra")
#package DT
#setwd("C:/Users/Eloi/Dropbox/Personnal/23andme")
#library(rmarkdown)
#dir_out="html_report"


setwd("C:/Users/Eloi/Dropbox/Personnal/23andme/HTML_Script")
source("config.ini")
library(rmarkdown)
render("files/index.Rmd", output_format=html_document(),
       output_file=file.path(html_out_dir, "index.html"),
       knit_root_dir = getwd())
