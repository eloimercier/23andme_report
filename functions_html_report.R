
#read 23andme raw file and get snpedia page name
read23andme <- function(file, force=F){
  out_file=gsub(".txt", ".Rds", file)
  if (force | !file.exists(out_file)){
    my_snps <- read.table(file, sep="\t", header=FALSE, colClasses=c("character", "character", "numeric", "character"), col.names=c("rsid", "chrom", "position", "genotype"))
    my_snps$chrom[my_snps$chrom=="MT"] <- "M"
    my_snps_genotype <- sapply(strsplit(my_snps$genotype,""), function(x){paste0("(",x[1],";",x[2],")")})
    my_snps_genotype <-  gsub("D","-", my_snps_genotype)
    my_snps$snpedia_id=paste0(my_snps$rsid,my_snps_genotype)
    # saveRDS(my_snps, file=out_file) 
  } else {
    my_snps=readRDS(out_file)
  }
  return(my_snps)
}

#get is_a_genotype list 
getGenotypeList <- function(save_file, force=FALSE){
  require(SNPediaR)
  #1. get snpedia info for the genoypes
  # is_genoype_file = file.path(out_dir, "Is_a_genotype.Rds")
  if (force | !file.exists(save_file)){ #check whether snpedia info already exists
    snpedia_is_a_genotype <- getCategoryElements (category = "Is_a_genotype")
    snpedia_is_a_genotype <- gsub("R","r", snpedia_is_a_genotype)
    saveRDS(snpedia_is_a_genotype, file=save_file)
  } else {
    snpedia_is_a_genotype=readRDS(save_file)
  }
  return(snpedia_is_a_genotype)
}

#get snpedia pages for a specific topic
getTopicPages <- function(topic=c("Is_a_medical_condition", "Is_a_medicine", "Topic"), save_file, force=FALSE){
  require(SNPediaR)
  if (!(topic %in% c("Is_a_medical_condition", "Is_a_medicine", "Topic"))) {stop("Topic is not accepted.")}
    # topic_snpedia_pages_file = file.path(out_dir, paste0(topic, "_snpedia_pages.Rds"))
  if (force | !file.exists(save_file)){
    snpedia_topic <- getCategoryElements (category = topic)
    topic_snpedia_pages <- getPages (titles = snpedia_topic)
    saveRDS(topic_snpedia_pages, file=save_file)
  } else {
    topic_snpedia_pages=readRDS(save_file)
  }
  return(topic_snpedia_pages)
}

#get snpedia pages for genotypes
getGenotypePages <- function(snepdia_id, save_file, force=FALSE){
  require(SNPediaR)
  # genotype_snpedia_pages_file = file.path(out_dir, paste0(topic, "_snpedia_pages.Rds"))
  if (force | !file.exists(save_file)){
    genotype_snpedia_pages <- getPages (titles = snepdia_id)
    saveRDS(genotype_snpedia_pages, file=save_file)
  } else {
    genotype_snpedia_pages=readRDS(save_file)
  }
  return(genotype_snpedia_pages)
}

#get the rsis from a snpedia page
findRsTags <- function (x) {
  x <- unlist (strsplit (x, split = "\n"))
  x <- grep ("rs", x, value = TRUE, ignore.case = T)
  x <- sapply(x, function(y){regmatches(y,regexpr("rs[0-9]+",y,ignore.case = T))})
  x <- unlist(x, use.names=F)
  unname(x)
}

#get snpedia infor for the snps in 23andme data referenced in snpedia
getMyGenotypeInfoInTopic <- function(my_snps,topic=c("Is_a_medical_condition", "Is_a_medicine", "Topic"),out_dir="my_23andme_html_report", snpedia_data="ext_data/snpedia", force=FALSE){
  require(SNPediaR)
  #1. get snpedia pages from topic
  topic_snpedia_pages <- getTopicPages(topic, save_file = file.path(snpedia_data, paste0(topic, "_snpedia_pages.Rds")), force)
  snps_topic <- lapply(topic_snpedia_pages, findRsTags)
  #2. get snpedia genotype list
  genotype_snpedia <- getGenotypeList(save_file=file.path(snpedia_data, "Is_a_genotype.Rds"), force)
  #2. get genotype per category in topic
  genotype_snpedia_filtered <- my_snps[my_snps$snpedia_id %in% genotype_snpedia , c("rsid","snpedia_id")]
  # saveRDS(genotype_snpedia_filtered, file=file.path(out_dir, "my_genotype_in_snpedia.Rds"))
  my_genotype_topic <- lapply(snps_topic, function(x){genotype_snpedia_filtered[genotype_snpedia_filtered$rsid %in% x,"snpedia_id"]  })
  my_genotype_topic <- my_genotype_topic[sapply(my_genotype_topic,length)!=0]
  #3. get genoypes info for my genotypes
  my_genotypes_in_snpedia_in_topic = unique(unlist(my_genotype_topic))
  my_genotype_pages <- getGenotypePages(my_genotypes_in_snpedia_in_topic, save_file=file.path(out_dir, paste0("my_genotypes_in_",topic, "_snpedia_pages.Rds")), force)
  names(my_genotype_pages) <- gsub("R","r", names(my_genotype_pages))
  #4. reassign genotype to the topic categories
  my_genotype_topic_pages=vector("list",length(my_genotype_topic))
  names(my_genotype_topic_pages)=names(my_genotype_topic)
  for (i in seq(my_genotype_topic_pages)){
    my_genotype_in_i <- my_genotype_topic[[i]]
    my_genotype_topic_pages[[i]]=unlist(my_genotype_pages[my_genotype_in_i] )
  }
  #5. extract tags
  my_genotype_snpedia_info <- lapply(my_genotype_topic_pages, function(l){sapply(l,extractGenotypeTags)} )
  return(my_genotype_snpedia_info)
}

#get a list of tags and create a data table
createDataTableFromGenotypeTagsList <- function(genotypeTagsList){
  category_summary <- data.frame(matrix(nrow=sum(sapply(genotypeTagsList, ncol)), ncol=5))
  colnames(category_summary)=c("Category", "Consequence", "rsid", "Genotype", "Summary")   
  category_summary[,"Category"]=rep(names(genotypeTagsList), sapply(genotypeTagsList, ncol))
  category_summary[,"Consequence"]=as.character(unlist(sapply(genotypeTagsList, function(x){x["repute",]})))
  category_summary[is.na(category_summary[,"Consequence"]),"Consequence"]="Unknown"
  category_summary[,"Magnitude"]=unlist(sapply(genotypeTagsList, function(x){x["magnitude",]}))
  category_summary[,"rsid"]=unlist(sapply(genotypeTagsList, function(x){paste0("rs",x["rsid",])}))
  category_summary[,"Genotype"]=unlist(sapply(genotypeTagsList, function(x){paste(x["allele1",],x["allele2",], sep="/")}))
  category_summary[,"Summary"]=unlist(sapply(genotypeTagsList, function(x){x["summary",]}))
  category_summary=category_summary[order(category_summary[,1],category_summary[,2]),]
  rownames(category_summary)=NULL
  return(category_summary)
}

#plot repute summary for a given category
plotSNPcount <- function(x){
  df=data.frame(rsid=paste0("rs",x["rsid",]),repute=x["repute",], magnitude=as.numeric(x["magnitude",]), stringsAsFactors=F)
  df[is.na(df[,"repute"]), "repute"]="Unknown"
  # df[is.na(df[,"magnitude"]), "magnitude"]=0
  #reassign unknown to bad or good based on keyword
  
  #ggplot2
    require(ggplot2)
    p<-ggplot(data=df, aes(x=repute, y=1, fill=magnitude, label=rsid)) +
      geom_bar(stat="identity",color="black") +
    geom_text(size = 3, position = position_stack(vjust = 0.5), color="black")+
      ylab(label="SNP count")
     p <-p + scale_fill_gradient(low="white", high="purple", limits=c(0,10), na.value = "grey80")
    p
}

#clean names to remove foreign character
cleanNames <- function(s, foreign_char="szþàáâãäåçèéêëìíîïðñòóôõöùúûüý", new_char="szyaaaaaaceeeeiiiidnooooouuuuy"){
  s1 <- chartr(foreign_char, new_char, s)
  return(s1)
}

#shorten names of categories
shortenCategoryNames <- function(x, max_char=15){
  x=cleanNames(x)
  short_names=substr(x, 1,max_char)
  if(any(duplicated(short_names)>0)){
    stop("Short names identical. Try increasing the max_char.")
  }
  return(short_names)
}

