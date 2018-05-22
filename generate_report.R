#this is what you need (and probably more)
#install pandoc
#install miktex
#install.packages("kableExtra")
#package DT
#library(rmarkdown)



scrip_directory="C:/my_path/23andme_report" #where the script is located,

setwd(scrip_directory)
source("config.ini")
library(rmarkdown)
render("files/index.Rmd", output_format=html_document(),
       output_file=file.path(html_out_dir, "index.html"),
       knit_root_dir = getwd())
