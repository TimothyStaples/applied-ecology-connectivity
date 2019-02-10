# ################################################################### ####
# Title: Applied research is on the rise but connectivity barriers    #### 
#        persist between four major subfields                         ####
# Author: Timothy L Staples                                           ####
# Details: Full data processing and analysis script                   ####  
# ################################################################### ####
# ####
# Working directory ####
rm(list=ls())
setwd("LOCATION OF THIS SCRIPT")
# Libraries ####
library(parallel)
library(stringi)
library(lme4)
library(vegan)
library(plotrix)
library(igraph)
library(mgcv)
library(tidytext)
library(SnowballC)
library(tidyr)
library(grr)

# Global functions ####

sapply(list.files(path="./Functions", pattern=".R", full.names=TRUE),
       source)

# Make first letter of any string capital, and all other letters lower case
cap.first<-function(string){

  split<-strsplit(string, " |-")

  sapply(split, function(split.entry){
    paste0(paste0(substr(toupper(split.entry), 1, 1),
                  substr(tolower(split.entry), 2, nchar(split.entry))),
           collapse=" ")
  })
}

# get a point on the circumference of a circle, given the x and y of origin,
# radius and desired degree from origin
circ.point<-function(x, y, r, deg){
  c(x + r * cos((deg * pi) / (180)),
    y + r * sin((deg * pi) / (180)))
}

# ####
# BASIC PROCESSING CODE (USED TO CREATE IMPORTED FILES) ####

#         CONVERT RAW TXT WoS DATA TO DATA-FRAMES - TIME INTENSIVE! ####
raw.wos.folders<-list.files("LOCATION OF WEB OF SCIENCE DOWNLOADS",
                            full.names=TRUE)

z<-raw.wos.folders[grepl("2017", raw.wos.folders)]

no_cores <- 4
cl <- makeCluster(no_cores)
setDefaultCluster(cl)
# export data and load libraries in cluster
clusterExport(cl=cl,
              varlist=c("raw.wos.folders"),
              envir=environment())

col.headings<-parSapply(cl, raw.wos.folders, function(z){
  print(z)
  sapply(list.files(z, pattern=".txt"), function(x){

    WOS.temp.raw<-readChar(paste0(z,"/",x),
                           file.info(paste0(z,"/",x))$size)

    wos.replace<-gsub("\n","GARBLEBREAK",WOS.temp.raw)

    wos.split<-unlist(strsplit(wos.replace,"GARBLEBREAK"))

    return(unique(substr(wos.split[3:length(wos.split)], 1, 2)))

  })
})
stopCluster(cl=cl)

codes<-unique(unlist(col.headings))
codes<-codes[!codes %in% c("  ", "")]

files<-paste0(rep(raw.wos.folders, sapply(raw.wos.folders, function(x){length(list.files(x))})),
              "/",
              unlist(sapply(raw.wos.folders, list.files)))

# For this newest update, we only need to add on 2017 papers
files<-files[grepl("2017", files)]

# remove files we've already processed
file.output<-paste0(substr(files,
                           regexpr("downloads/", files)+10,
                           regexpr("savedrecs", files)-2),
                    " - ", substr(files,
                                  regexpr("savedrecs", files)+10,
                                  regexpr("\\.txt", files)-1), ".csv")

processed.dir<-"./Processed csv"

files<-files[!file.output %in% list.files(processed.dir)]

no_cores <- 4
cl <- makeCluster(no_cores)
setDefaultCluster(cl)
# export data and load libraries in cluster
clusterExport(cl=cl,
              varlist=c("codes", "files", "processed.dir"),
              envir=environment())
parSapply(cl, files,
          function(x){

            print(x)

            WOS.temp.raw<-readChar(x, file.info(x)$size)

            wos.replace<-gsub("\n","GARBLEBREAK",WOS.temp.raw)

            wos.split<-unlist(strsplit(wos.replace,"GARBLEBREAK"))

            # SEPARATE TEXT REFERENCES #

            frame<-as.data.frame(matrix(NA, nrow=500, ncol=length(codes)))
            colnames(frame)<-codes
            # If it's a field that we want
            for(i in 3:(length(wos.split)-2)){

              # extract out row of text string
              base.string<-wos.split[i]

              # does the string have one of our categories?
              if(substr(base.string,1,2) %in% codes){

                # if yes, we want to put it in the right spot and cut out
                # the code

                # find right column
                temp.col<-which(colnames(frame)==substr(base.string,1,2))

                # get the next empty row (has NA) in that column
                if(colnames(frame)[temp.col]=="PT"){
                  temp.row=min((1:dim(frame)[1])[is.na(frame[,"PT"])])
                } else {
                  first.empty.row<-(1:dim(frame)[1])[apply(frame, 1, function(x){sum(is.na(x))})==
                                                       dim(frame)[2]][1]

                  temp.row<-ifelse(is.na(first.empty.row),
                                   dim(frame)[1],
                                   first.empty.row-1)
                }

                # subset out code from text)
                temp.string<-substr(base.string, 4, nchar(base.string))

                # find first character that's not a space and start string
                # from there.
                prep.string<-substr(temp.string,
                                    gregexpr("[^ ]*",
                                             temp.string)[[1]][attr(gregexpr("[^ ]*",
                                                                             temp.string)[[1]],
                                                                    "match.length")>0][1],
                                    nchar(temp.string))

                frame[temp.row, temp.col]<-prep.string

              }

              # if string doesn't have a category
              if(!substr(base.string,1,2) %in% codes){

                # there's a space between each paper, skip this
                if(nchar(base.string)==0){next}

                # we want to add it to the last cell to be done
                # we'll need to cycle back through the previous strings to get one that has
                # a code in it, starting with the
                j=i-1
                repeat{
                  if(substr(wos.split[j],1,2) %in% codes){
                    temp.col=which(colnames(frame)==substr(wos.split[j],1,2))
                    break
                  } else {j=j-1}
                }

                # in this case we want the last non-empty row, but we need to be careful that
                # we're not putting data onto the wrong row (say when a paper doesn't have entry
                # for a given term)

                # the way we'll do this is rely on the fact that every paper will have at least
                # one author, or be designated as a "J" journal before we get to this point.
                # then we'll look for the first row that is completely NA, then assign the
                # row as the previous one. This should ensure that if we've moved passed a
                # paper with missing info, we don't mix up our records.
                first.empty.row<-(1:dim(frame)[1])[apply(frame, 1, function(x){sum(is.na(x))})==
                                                     dim(frame)[2]][1]

                temp.row<-ifelse(is.na(first.empty.row),
                                 dim(frame)[1],
                                 first.empty.row-1)

                # get actual string
                base.string<-wos.split[i]

                # find first character that's not a space and start string
                # from there.
                prep.string<-substr(base.string,
                                    gregexpr("[^ ]*",
                                             base.string)[[1]][attr(gregexpr("[^ ]*",
                                                                             base.string)[[1]],
                                                                    "match.length")>0][1],
                                    nchar(base.string))

                frame[temp.row, temp.col]=paste0(frame[temp.row, temp.col], "; ", prep.string)

              }
            }

            write.csv(frame, paste0(processed.dir,"/",
                                    substr(x,
                                           regexpr("downloads/", x)+10,
                                           regexpr("savedrecs", x)-2),
                                    " - ",
                                    substr(x,
                                           regexpr("savedrecs", x)+10,
                                           regexpr("\\.txt", x)-1),
                                    ".csv"))

          })
stopCluster(cl=cl)

# Combine data files together
full.data.list<-lapply(list.files(processed.dir),
                       function(x){
                         print(x)
                         return(read.csv(paste0(processed.dir, "/",x),
                                         row.names=1,
                                         colClasses = "character"))
                       })

# do our dimensions match up?
dl.dims<-do.call("rbind", lapply(full.data.list, dim))
rownames(dl.dims)<-list.files(processed.dir)
table(dl.dims)

# subset list entries with extra columns down to the base 63
base<-full.data.list[[which(dl.dims[,2]==63)[1]]]

# these have one column not present in the other tables, "BS"

base.cols<-colnames(base)[colnames(base) != "BS"]

full.data.list.syn<-do.call("rbind", lapply(full.data.list, function(x){

  if(sum(!colnames(x) %in% base.cols) == 0){
    return(x)
  } else {

  return(x[,base.cols])
  }
  }))

#full.dataframe<-do.call("rbind", full.data.list)

#frame<-full.dataframe

write.csv(full.data.list.syn, "./Data/full processed data 1.csv")

# SUBSET PROCESSED DATA ####

gc()
frame<-read.csv("./Data/full processed data 1.csv", row.names=1,
                header=TRUE, stringsAsFactors = FALSE)

# list of unique entries, with titles and counts
journals<-as.data.frame(table(frame$SO))
journals<-journals[order(journals$Freq, decreasing=TRUE),]
journals[1:20,]
hist(log(journals$Freq))

# include only publications with more than 100 papers
journals.sub<-journals[journals$Freq >= 100,]

# exclude non-ecology journals that were originally included in download
non.ecology<-c("PLOS ONE",
               "PROCEEDINGS OF THE NATIONAL ACADEMY OF SCIENCES OF THE UNITED STATES OF; AMERICA",
               "NATURE",
               "SCIENCE",
               "NATURE COMMUNICATIONS",
               "PROCEEDINGS OF THE ROYAL SOCIETY B-BIOLOGICAL SCIENCES",
               "PEERJ")

frame.sub<-frame[frame$SO %in% journals.sub$Var1 &
                 !frame$SO %in% non.ecology &
                 frame$PY >= 1990,]

rm(frame)
gc()
# ####
#         CITATION MATCHING - MEMORY INTENSIVE! ####
#                 PREP PAPER-LEVEL DATA ####

# get rid of NA rows generated from read-in process
# these are rows with NAs all the way across
summary(rowSums(is.na(frame.sub))<dim(frame.sub)[2])

frame.sub<-frame.sub[rowSums(is.na(frame.sub))<dim(frame.sub)[2],];gc()

# Get rid of abstract to save on file size
frame.sub<-frame.sub[,colnames(frame.sub) != "AB"]

# Remove papers before 1990 - some years have an issue attached to them which
# we need to remove
frame.sub<-frame.sub[frame.sub$PY>=1990,]

#     Cut out book chapters and conference proceedings #

# First off we need to cut out book chapters, conference proceedings etc.
sort(table(frame.sub$DT), decreasing=TRUE)

# So we want to keep:
# Article, Review, Editorial Material, Note and Letter
frame.sub<-frame.sub[frame.sub$DT %in% c("Article", "Review",
                                         "Note", "Letter") &
                       frame.sub$PT=="J",]

# check publication over time to find journals with changed names
write.csv(table(frame.sub$SO, frame.sub$PY),
          "./Outputs/Journal rate by year - find old names.csv")

# make all journal names and abbreviations lower case
frame.sub$SO<-tolower(frame.sub$SO)
frame.sub$J9<-tolower(frame.sub$J9)

changed.journals<-read.csv("./Data/journal name changes.csv",
                           stringsAsFactors=FALSE)
changed.journals<-as.data.frame(do.call("cbind", lapply(changed.journals, tolower)),
                                stringsAsFactors=FALSE)

changed.journals$old.J9<-tolower(frame.sub$J9[match(changed.journals$old.name,
                                                    frame.sub$SO)])
changed.journals$new.J9<-tolower(frame.sub$J9[match(changed.journals$new.name,
                                                    frame.sub$SO)])

frame.sub$SO<-ifelse(frame.sub$SO %in% changed.journals$old.name,
                     changed.journals$new.name[match(frame.sub$SO,
                                                     changed.journals$old.name)],
                     frame.sub$SO)

frame.sub$J9<-ifelse(frame.sub$J9 %in% changed.journals$old.J9,
                     changed.journals$new.J9[match(frame.sub$J9,
                                                   changed.journals$old.J9)],
                     frame.sub$J9)

# the only NAs here are for Ekologica Bratislava, which doesn't have an abbreviation
# in WoS. We may not need it, but I found one from Springer that I'll add in
frame.sub$J9[frame.sub$SO=="ekologia bratislava"]="ekologia bratisl"

# CREATE REFERENCE CODES FOR CITATIONS ####

refcode.doi<-tolower(ifelse(!is.na(frame.sub$DI),
                            frame.sub$DI,
                            "no doi"))

# where there are duplicate DOIs, we'll need to rely on page-string citation
# match instead, because I can't search through ~4000 duplicates to find which
# paper has the wrong DOI.
doi.table<-table(refcode.doi[refcode.doi!="no doi"])
doi.dupes<-names(doi.table)[doi.table>1]
refcode.doi[refcode.doi %in% doi.dupes] = "no doi"

refcode.pages<-tolower(reference.code.pages(frame.sub))

# there's also some duplicate page strings, so we'll remove any of those as well
refcode.pages[refcode.pages %in% refcode.pages[duplicated(refcode.pages)]]="duplicate pagestring"
summary(refcode.pages=="duplicate pagestring")
summary(is.na(refcode.pages))

#                 CUT OUT DUPLICATE PAPERS ####

# now we want to remove duplicate paper entries, but we want to see if they have
# non-overlapping references (like each version has it's own DOI etc).

# duplicate reference codes
dupe.doi<-duplicated(refcode.doi) & refcode.doi != "no doi"
summary(dupe.doi)

dupe.pages<-duplicated(refcode.pages) & refcode.doi == "no doi"
summary(dupe.pages)

frame.sub<-frame.sub[!dupe.pages,]
refcode.doi<-refcode.doi[!dupe.pages]
refcode.pages<-refcode.pages[!dupe.pages]

# from now on we'll be working with paper numbers, so let's set up an index
# so we don't muck up (rownames should work, but let's be careful)

frame.sub$index<-1:dim(frame.sub)[1]
names(refcode.doi) <- frame.sub$index
names(refcode.pages) <- frame.sub$index

gc()
#                 GENERATE CITED REF LEVEL DATA ####

D# Prepare dataframe of cited references - split each string by separator (";")
CR.split<-strsplit(tolower(frame.sub$CR), split=";"); gc()

# Combine these together with correct index number
CR.df<-data.frame(CR=unlist(CR.split),
                  index=rep(1:length(CR.split), sapply(CR.split, length)),
                  stringsAsFactors=FALSE); rm(CR.split); gc()

# so now we need to do some fancy work so we can just search for matches, rather
# than searching for substring matches

# first off, remove any spaces at the beginning of the string
CR.df$CR<-ifelse(substr(CR.df$CR,1,1)==" ",
                 substr(CR.df$CR,2, nchar(CR.df$CR)),
                 CR.df$CR); gc()

# first off, make a separate DOI column, otherwise NA
CR.df$CR.doi<-ifelse(regexpr("doi ", CR.df$CR) == -1,
                     NA,
                     substr(CR.df$CR,
                            regexpr("doi ", CR.df$CR)+4,
                            nchar(CR.df$CR))); gc()

# next we'll create an reference using year, journal abbreviation, volume and
# page number. We'll avoid author info because that's a big mess at the moment.
CR.df$CR<-substr(CR.df$CR,
                 1,
                 ifelse(is.na(CR.df$CR.doi),
                        nchar(CR.df$CR),
                        regexpr(", doi", CR.df$CR)-1)); gc()

colnames(CR.df)<-c("CR.page.string","index","CR.doi")

# finally, we'll delete citations with no page, volume or year information

# first volume and page
CR.df$CR.page.string[!(grepl(", p", CR.df$CR.page.string) &
                         grepl(", v", CR.df$CR.page.string))]=NA; gc()

# now year
CR.df$CR.page.string[!grepl(paste0(paste0(1800:2016, collapse="|")),
                            CR.df$CR.page.string)]=NA; gc()

# now we need to get rid of author information
CR.df$CR.page.string<-substr(CR.df$CR.page.string,
                             regexpr(paste0(1800:2016, collapse="|"),
                                     CR.df$CR.page.string),
                             nchar(CR.df$CR.page.string))
head(CR.df)

write.csv(CR.df, "./Data/citation-level dataframe.csv")
write.csv(refcode.doi, "./Data/DOI reference code.csv")
write.csv(refcode.pages, "./Data/Pagestring reference code.csv")
write.csv(frame.sub, "./Data/subset paper dataframe.csv")

# now we have our dataframe of citations, and our reference code for which paper
# numbers should match up perfectly

#                 MATCH CITATIONS ####

CR.df<-read.csv("./Data/citation-level dataframe.csv", header=TRUE, row.names=1,
                stringsAsFactors=FALSE)
refcode.doi<-read.csv("./Data/DOI reference code.csv", header=TRUE, row.names=1,
                      stringsAsFactors=FALSE)[,1]
refcode.pages<-read.csv("./Data/Pagestring reference code.csv", header=TRUE, row.names=1,
                        stringsAsFactors=FALSE)[,1]

summary(CR.df$CR.doi==tolower("10.1111/j.1526-100X.1996.tb00112.x"))

summary(as.numeric(names(refcode.doi)))
summary(as.numeric(names(refcode.pages)))


# match by DOI first
doi.cit.df<-matches(x=refcode.doi,
                    y=CR.df$CR.doi,
                    all.x=FALSE,
                    all.y=TRUE)

# then by page string
pages.cit.df<-matches(x=refcode.pages,
                      y=CR.df$CR.page.string,
                      all.x=FALSE,
                      all.y=TRUE)

# now our y indexes don't match our paper numbers, instead they are the row
# numbers from our citation database

# so now let's merge across the correct paper indexes
CR.index<-data.frame(y=rownames(CR.df),
                     index=CR.df$index)

doi.cit.df<-merge(doi.cit.df, CR.index, all.x=TRUE, all.y=FALSE)
pages.cit.df<-merge(pages.cit.df, CR.index, all.x=TRUE, all.y=FALSE)

doi.cit.df.sub<-doi.cit.df[complete.cases(doi.cit.df),]
pages.cit.df.sub<-pages.cit.df[complete.cases(pages.cit.df),]

gc()
# so we start with our DOI matches, as they're our "highest certainty" matches
doi.cit.df.sub$match.type="DOI"

cit.df<-doi.cit.df.sub[,c("x","index","match.type")]

# then we can add in our pages, only adding in references that didn't match by
# doi
summary(pages.cit.df.sub$y %in% doi.cit.df.sub$y)

pages.cit.df.unique<-pages.cit.df.sub[!pages.cit.df.sub$y %in% doi.cit.df.sub$y,]
pages.cit.df.unique$match.type="pagestring"

head(pages.cit.df.unique)
# combine together
cit.df<-rbind(cit.df, pages.cit.df.unique[,c("x","index","match.type")])
colnames(cit.df)<-c("cited","citing","match.type")

cit.df<-cit.df[order(cit.df$citing),]

cit.df<-cit.df[, c(2,1,3)]

write.csv(cit.df, "./Data/cited references dataframe.csv")

# ####
# EVALUATE CITATION MAPPING - MEMORY INTENSIVE! ####

CR.df<-read.csv("./Data/citation-level dataframe.csv", row.names=1,
                stringsAsFactors=FALSE)
refcode.doi<-read.csv("./Data/DOI reference code.csv", header=TRUE, row.names=1,
                      stringsAsFactors=FALSE)[,1]
refcode.pages<-read.csv("./Data/Pagestring reference code.csv", header=TRUE, row.names=1,
                        stringsAsFactors=FALSE)[,1]
frame.ecology<-read.csv("./Data/ecology papers dataframe.csv",
                        stringsAsFactors=FALSE, row.names=1)
#cit.df<-read.csv("./Data/cited references dataframe.csv", stringsAsFactors=FALSE,
#                  row.names=1)

gc()

CR.df<-CR.df[CR.df$index %in% frame.ecology$index,]
refcode.doi<-refcode.doi[1:length(refcode.doi) %in% frame.ecology$index]
refcode.pages<-refcode.pages[1:length(refcode.pages) %in% frame.ecology$index]

# match by DOI first
doi.match<-matches(x=refcode.doi,
                   y=CR.df$CR.doi,
                   all.x=FALSE,
                   all.y=TRUE)

pages.cit<-matches(x=refcode.pages,
                   y=CR.df$CR.page.string,
                   all.x=FALSE,
                   all.y=TRUE)

CR.df$doi.match<-doi.match[order(doi.match$y),1] ; gc()
CR.df$pages.match<-pages.cit[order(pages.cit$y),1] ; gc()

summary(CR.df)

# Okay, now the summary stats

# total citations
evaluation.stats<-dim(CR.df)[1]
names(evaluation.stats)[1]<-"all.citations"

CR.save<-CR.df

# first, get number of matched citations
matched<-!is.na(CR.df$doi.match) | !is.na(CR.df$pages.match)
summary(matched)

# remove matched citations
evaluation.stats<-c(evaluation.stats,
                    sum(matched))
names(evaluation.stats)[2]<-"matched.citations"

CR.df<-CR.df[!matched,]

# citations with incomplete information
bad.info<-is.na(CR.df$CR.doi) & is.na(CR.df$CR.page.string)
summary(bad.info)

evaluation.stats<-c(evaluation.stats,
                    sum(bad.info))
names(evaluation.stats)[3]<-"missing.info"

CR.df<-CR.df[!bad.info,]; gc()

# of citations with complete information:
# citing pre-1980 (outside my network)
pre.1980<-grepl(paste0(1800:1979, collapse="|"),
                CR.df$CR.page.string)
summary(pre.1980)

evaluation.stats<-c(evaluation.stats,
                    sum(pre.1980))
names(evaluation.stats)[4]<-"pre.1980"

CR.df<-CR.df[!pre.1980,]; gc()

# sub out remaining citations and get table of journal names
cited.journal<-substr(CR.df$CR.page.string, 7, nchar(CR.df$CR.page.string))
summary(is.na(CR.df$CR.page.string))
summary(is.na(cited.journal))

non.ecol.SO<-c("science", "plos one", "peerj", "nature", "nature communications",
               "proceedings of the national academy of sciences of the united states of; america")
non.ecol.J9<-unique(tolower(frame.ecology$J9[tolower(frame.ecology$SO) %in% non.ecol.SO]))

cited.journal.df<-data.frame(journal=substr(cited.journal, 1,
                                            regexpr(", v", cited.journal)-1))
cited.journal.df$is.ecology=cited.journal.df$journal %in% unique(tolower(frame.ecology$J9))

write.csv(CR.df[cited.journal.df$is.ecology,],
          "./Data/missed citations to ecology journals.csv")

evaluation.stats<-c(evaluation.stats,
                    sum(!cited.journal.df$is.ecology),
                    sum(cited.journal.df$is.ecology))
names(evaluation.stats)[5]<-"non.ecol.journals"
names(evaluation.stats)[6]<-"ecol.journals"

# do we have all the citations now?
sum(evaluation.stats[-1]) == evaluation.stats[1]
# YES

write.csv(as.data.frame(evaluation.stats),
          "./Outputs/Citation match evaluation.csv")

# ####
# ADDITIONAL PROCESSING <- START WITH PROCESSED DATA HERE! ####
#               READ IN DATA ####

frame.ecology<-read.csv("./Data/subset paper dataframe.csv", row.names=1,
                        stringsAsFactors = FALSE)
cit.df<-read.csv("./Data/cited references dataframe.csv", row.names=1,
                 stringsAsFactors=FALSE)

# remove some non-ecology journals that slipped through our initial process
frame.ecology <- frame.ecology[!frame.ecology$SO %in% c("peerj",
                                                        "plos one",
                                                        "nature communications") &
                                frame.ecology$DT != "Editorial Material",]

cit.df <- cit.df[cit.df$citing %in% frame.ecology$index &
                 cit.df$cited %in% frame.ecology$index, ]

# Remove papers that don't cite at least one other paper in our list
non.citers<-frame.ecology$index[!frame.ecology$index %in% cit.df$citing & 
                                !frame.ecology$index %in% cit.df$cited]

frame.ecology<-frame.ecology[!frame.ecology$index %in% non.citers, ]
cit.df<-cit.df[!cit.df$citing %in% non.citers |
               !cit.df$cited %in% non.citers, ]  

summary(unique(cit.df$citing) %in% frame.ecology$index)
summary(frame.ecology$index %in% unique(cit.df$citing))
summary(frame.ecology$index %in% unique(cit.df$cited))
summary(unique(cit.df$cited) %in% frame.ecology$index)

summary(frame.ecology$index %in% c(unique(cit.df$citing), unique(cit.df$cited)))

#               IDENTIFY PAPER GROUPS ####
#                               COMBINE TITLE AND KEYWORDS ####

# I want to combine the title and keywords of a paper together to use as a
# search string for particular groups of papers.

# for starters, I need to remove any ';' from the title field that came about
# from titles split across multiple lines
titles<-tolower(gsub(";","", frame.ecology$TI))
head(titles)

paper.words<-paste0(titles, ifelse(is.na(frame.ecology$DE),
                                   "",
                                   paste0("; ", tolower(frame.ecology$DE))))
paper.words<-gsub(";;",";", paper.words)
paper.words<-gsub("[[:punct:]]","", paper.words)

#                               HIGHLIGHT APPLIED PAPERS ####

# now we want to highligh applied papers using a number of keywords
# we'll do them separately so I can see how many papers fall into each keyword
rest.words<-c("restor")
cons.words<-c("conserv")
clim.words<-c("climate change")
invas.words<-c("invasi")

applied.words<-list(rest.words,
                    cons.words,
                    clim.words,
                    invas.words)

names(applied.words)<-c("restoration","conservation","climatechange", "invasion")

applied.words.df<-lapply(list(rest.words,
                              cons.words,
                              clim.words,
                              invas.words), function(x){

                                sapply(x, function(y){grepl(y, paper.words)})

                              })

sampled.papers<-lapply(applied.words.df, function(x){

  apply(x, 2, function(y){
    frame.ecology[y, c("TI")][sample(1:sum(y), ifelse(sum(y)>20, 20, sum(y)))]

  })
})

searchword.performance<-do.call("rbind", lapply(1:length(applied.words.df),
                                                function(x){

                                                  word.df<-applied.words.df[[x]]

                                                  temp<-do.call("rbind", lapply(1:dim(word.df)[2], function(y){

                                                    data.frame(count=sum(word.df[,y]),
                                                               unique.count=sum(word.df[rowSums(word.df[,-y])==0, y]),
                                                               searchword=colnames(word.df)[y])
                                                  }))

                                                  temp$group=c("restoration","conservation","climatechange","invasion")[x]
                                                  return(temp)
                                                }))

word.table.data<-do.call("rbind", lapply(unique(searchword.performance$group),
                                         function(x){
                                           print(x)
                                           temp<-searchword.performance[searchword.performance$group==x,]

                                           data.frame(count=paste0(temp$count, " (", temp$unique.count, ")"),
                                                      group=temp$searchword)
                                         }))

write.csv(searchword.performance, "./Outputs/searchword performance table.csv")
write.csv(word.table.data, "./Outputs/searchword data table.csv")

paper.groups<-as.data.frame(sapply(applied.words.df, function(x){ifelse(rowSums(x)>0,TRUE,FALSE)}))
colnames(paper.groups)<-c("restoration","conservation","climatechange",
                          "invasion")

paper.groups$ecology<-ifelse(rowSums(paper.groups)==0,TRUE,FALSE)

sum(rowSums(paper.groups[,1:4])>0)
colSums(paper.groups)

#               INTEGRATE SUBFIELDS INTO CITATION-LEVEL DATA ####

# the first thing we want to do is quantify how many citations are within each
# group and between each group

# let's start by getting overall citation counts for each paper, which we
# can match up to the official WOS citation counts
cited.count<-table(cit.df[,"cited"])

frame.ecology$cited.count=0
frame.ecology$cited.count[match(as.numeric(names(cited.count)),
                                frame.ecology$index)]=cited.count

citing.count<-table(cit.df[,"citing"])

frame.ecology$citing.count=0
frame.ecology$citing.count[match(as.numeric(names(citing.count)),
                                 frame.ecology$index)]=citing.count

hist(frame.ecology$citing.count, breaks=20)
hist(log(frame.ecology$cited.count), breaks=20)

# how does citing count compare to length of the "cited references" field
with(frame.ecology, plot(y=log(cited.count), x=log(TC)))
with(frame.ecology, plot(y=log(citing.count), x=log(NR)))

# Let's look at the highest cited papers and see if they make sense
frame.ecology$TI[order(frame.ecology$cited.count, decreasing=TRUE)][1:20]
frame.ecology$TC[order(frame.ecology$cited.count, decreasing=TRUE)][1:20]

# next get a dataframe of whether each citation was made to or from a
# particular group
paper.groups.citing<-paper.groups
paper.groups.cited<-paper.groups

cit.df<-merge(cit.df,
              data.frame(citing=frame.ecology$index,
                         citing.year=frame.ecology$PY),
              all.x=TRUE, all.y=FALSE)
head(cit.df)

paper.groups.citing$citing=frame.ecology$index
colnames(paper.groups.citing)[1:5]<-paste0(colnames(paper.groups)[1:5],".citing")

paper.groups.cited$cited=frame.ecology$index
colnames(paper.groups.cited)[1:5]<-paste0(colnames(paper.groups)[1:5],".cited")

cit.df.groups<-merge(cit.df, paper.groups.citing, all.x=TRUE)
cit.df.groups<-cit.df.groups[,colnames(cit.df.groups)!="index"]
cit.df.groups<-merge(cit.df.groups, paper.groups.cited, all.x=TRUE)

# ####
# MODELLING ####
#               OVERALL RATES ####
#                               PUB RATE ####

pub.rate<-do.call("rbind", lapply(1990:2017, function(x){
  temp<-rowSums(paper.groups[frame.ecology$PY==x,1:4]>0)
  return(data.frame(pub.rate= sum(temp) / length(temp),
                    PY=x))
}))

pub.rate.raw<-sapply(1990:2017, function(x){
  temp<-rowSums(paper.groups[frame.ecology$PY==x,1:4]>0)
  return(sum(temp))
})

total.pub.rate.bGAM<-gam(pub.rate ~ s(PY, bs="cr", k=4),
                         family=betar(),
                         data=pub.rate,
                         method="REML")

summary(total.pub.rate.bGAM)
gam.check(total.pub.rate.bGAM)

#                               CIT RATE ####

cit.rate<-do.call("rbind", lapply(1990:2017, function(x){
  temp<-rowSums(cit.df.groups[cit.df.groups$citing.year == x,
                              grepl("\\.cited", colnames(cit.df.groups)) &
                                !grepl("ecology", colnames(cit.df.groups))])
  return(data.frame(cit.rate= sum(temp) / length(temp),
                    PY=x))

}))

overall.cit.rate.bGAM<-gam(cit.rate ~ s(PY, bs="cr", k=4),
                           family=betar(),
                           data=cit.rate,
                           method="REML")
summary(overall.cit.rate.bGAM)
gam.check(overall.cit.rate.bGAM)

#               SUB-FIELD RATES ####
#                               PUB RATE ####

# summarise publication rate by year
group.pub.rate<-do.call("rbind", lapply(1990:2017, function(x){
  temp<-colSums(paper.groups[frame.ecology$PY==x, 1:4])
  return(data.frame(pub.rate=temp / sum(frame.ecology$PY==x),
                    group=names(temp),
                    PY=rep(x, 4)))
}))

group.pub.rate.bGAM<-gam(pub.rate ~ group + s(PY, bs="cr", k=4, by=group),
                         family=betar(),
                         data=group.pub.rate,
                         method="REML")
summary(group.pub.rate.bGAM)
gam.check(group.pub.rate.bGAM)

#                               CIT RATE ####

cited<-"restoration"
year<-2017
cit.summ<-do.call("rbind", lapply(colnames(paper.groups)[1:4], function(cited){
  print(cited)

  cited.sub<-cit.df.groups[,paste0(cited,".cited")]

  temp<-do.call("rbind",
                lapply(1990:2017, function(year){
                  data.frame(cit.prop=sum(cited.sub & cit.df.groups[,"citing.year"]==year) / 
                               sum(cit.df.groups[,"citing.year"]==year),
                             PY=year)}))

  temp$group=cited

  return(temp)
}))

cit.summ$group<-as.factor(cit.summ$group)

group.cit.rate.bGAM<-gam(cit.prop ~ group + s(PY, bs="cr", k=4, by=group),
                         family=betar(),
                         data=cit.summ,
                         method="REML")
summary(group.cit.rate.bGAM)
gam.check(group.cit.rate.bGAM)

#               CITATION PROPORTION RATES ####

# I'm going to change these models into hurdle models like I have used before,
# where we start by modelling the probability of a paper citing another paper,
# and then a secondary model with only those who do cite a paper, looking at
# the count or proportion of papers being cited.

bween.group.citrate.binary<-do.call("rbind",
                                    lapply(colnames(paper.groups)[1:4],
                                           function(x){
                     
                     # for each group
                     do.call("rbind", lapply(colnames(paper.groups)[1:4], function(xx){
                       
                       # subset out just papers from group 1
                       cit.group.sub<-cit.df.groups[cit.df.groups[, paste0(x, ".citing")], ]
                       
                       # now split papers by year, then get proportion that cite var 2
                       year.rate<-do.call("rbind", lapply(1990:2017, function(y){
                         
                         year<-y
                         # subset just the year we want
                         y<-cit.group.sub[cit.group.sub$citing.year==y,]
                         
                         cite.prop<-tapply(y[,paste0(xx, ".cited")],
                                           y$citing,
                                           function(z){
                                             ifelse(sum(z)>0, 1, 0)
                                           })
                         
                         return(data.frame(PY=year,
                                           cite.1=sum(cite.prop),
                                           cite.0=length(cite.prop)-sum(cite.prop)))
                       }))
                       
                       temp.df<-data.frame(citing.group=rep(x, dim(year.rate)[1]),
                                           cited.group=rep(xx, dim(year.rate)[1]))
                       
                       return(cbind(year.rate, temp.df))
                       
                     }))
                   }))

summary(unique(cit.df$cited) %in% frame.ecology$index)
summary(frame.ecology$index %in% unique(cit.df$cited))


bween.group.citrate.count<-do.call("rbind", lapply(colnames(paper.groups)[1:4],
                                  function(x){
                                    
                                    # for each group
                                    do.call("rbind", 
                                            lapply(colnames(paper.groups)[1:4], 
                                                   function(xx){
                                      
                                      # subset out just papers from group 1
                                      cit.group.sub<-cit.df.groups[cit.df.groups[, paste0(x, ".citing")], ]
                                      
                                      # now split papers by year, then get proportion that cite var 2
                                      year.rate<-do.call("rbind", lapply(1990:2017, function(y){
                                        
                                        year<-y
                                        # subset just the year we want
                                        y<-cit.group.sub[cit.group.sub$citing.year==y,]
                                        
                                        # get just papers from group x that cite 
                                        # a paper from group xx
                                        y<-y[y$citing %in% y$citing[y[,paste0(xx, ".cited")]],]

                                        
                                        cite.count<-do.call("rbind", lapply(unique(y$citing), function(z){
                                          
                                          data.frame(cite.count=sum(y[y$citing==z, paste0(xx, ".cited")]),
                                                     index=z)
                                           
                                        }))
                                        
                                        cite.count<-merge(cite.count,
                                                          frame.ecology[,c("index","NR")],
                                                          all.x=TRUE, all.y=FALSE)
                                        
                                        cite.prop<-cite.count$cite.count/
                                                   cite.count$NR
                                        
                                        return(data.frame(PY=rep(year, length(cite.prop)),
                                                          cite.prop=cite.prop,
                                                          index=cite.count$index))
                                      }))
                                      
                                      temp.df<-data.frame(citing.group=rep(x, dim(year.rate)[1]),
                                                          cited.group=rep(xx, dim(year.rate)[1]))
                                      
                                      return(cbind(year.rate, temp.df))
                                      
                                    }))
                                  }))


# binary gam for proportions
bween.group.citrate.binary$group<-as.factor(paste0(bween.group.citrate.binary$citing.group,
                                           ":",
                                           bween.group.citrate.binary$cited.group))

tapply(rowSums(bween.group.citrate.binary[,c("cite.0","cite.1")]),
       bween.group.citrate.binary$group, sum)

bween.citrate.bin.gam<-gam(cbind(cite.1, cite.0) ~ s(PY, bs="cr", k=4, by=group) + group,
                           data=bween.group.citrate.binary,
                           family=binomial(),
                           method="REML")
summary(bween.citrate.bin.gam)

write.csv(cbind(summary(bween.citrate.bin.gam)$p.table,
                summary(bween.citrate.bin.gam)$s.table),
                "./Outputs/binary citation probability gam summ table.csv")

# beta gam for counts
bween.group.citrate.count$group<-as.factor(paste0(bween.group.citrate.count$citing.group,
                                                  ":",
                                                  bween.group.citrate.count$cited.group))

bween.citrate.beta.gam<-gam(log(cite.prop) ~ s(PY, bs="cr", k=4, by=group) + group,
                            data=bween.group.citrate.count,
                            family=gaussian,
                            method="REML")
gam.check(bween.citrate.beta.gam)

# model with aggregated yearly means - does not differ substantially from paper-level model
bween.group.citrate.count.agg<-do.call("rbind", lapply(split(bween.group.citrate.count,
                                                             f=bween.group.citrate.count$group),
                                                       function(x){
                                                         
                  data.frame(cite.prop = sapply(1990:2017, function(y){
                    ifelse(sum(x$PY == y) == 0, 0, mean(x$cite.prop[x$PY==y]))
                  }),
                     PY=1990:2017,
                  group=x$group[1])
                                                         }))

bween.citrate.beta.gam.agg<-gam(cite.prop ~ s(PY, bs="cr", k=4, by=group) + group,
                                data=bween.group.citrate.count.agg,
                                family=betar(),
                                method="REML")
gam.check(bween.citrate.beta.gam.agg)
write.csv(cbind(summary(bween.citrate.beta.gam)$p.table,
                summary(bween.citrate.beta.gam)$s.table),
          "./Outputs/beta citation rate gam summ table.csv")

#               JOURNAL ORDINATION ####
frame.ecology$SO.fact<-as.factor(frame.ecology$SO)
group.journal<-do.call("cbind", lapply(paper.groups[,1:4], function(x){
  
  with(frame.ecology[x,], 
       tapply(SO.fact, as.factor(SO.fact), length))
}))

# replace NAs with 0s
group.journal[is.na(group.journal)]=0
group.journal <- group.journal[rowSums(group.journal)>0,]

# import impact factors
ifactor<-read.csv("./Data/ordination journals impact factor.csv", stringsAsFactors=FALSE)
ifactor<-ifactor[!is.na(ifactor$ifactor), ]

journal.prop.mat<-group.journal / sapply(rownames(group.journal), function(x){
  sum(frame.ecology$SO == x)})
journal.sub.mat<-journal.prop.mat[rowSums(journal.prop.mat[,1:4]) > 0.4,]

journal.sub.mat <- as.data.frame(journal.sub.mat)
journal.sub.mat$abbrev <- c("AC", "AVS", "AQI", "B&C", "BC", "BI", "CB", "DD", "EMR",
                            "GCB", "JNC", "NAJ", "OR", "RE")

write.csv(journal.sub.mat, "./Outputs/Table 1.csv")

journal.prop.mat<-cbind(journal.prop.mat,
                        sapply(rownames(journal.prop.mat),
                               function(x){sum(frame.ecology$SO == x)}))

rownames(journal.prop.mat) <- cap.first(rownames(journal.prop.mat))
write.csv(journal.prop.mat, "./Outputs/Table S8.csv")

rownames(journal.prop.mat) <- tolower(rownames(journal.prop.mat))

ifactor <- ifactor[ifactor$SO %in% tolower(rownames(journal.prop.mat)),]
journal.prop.mat.ord<-journal.prop.mat[tolower(rownames(journal.prop.mat)) %in% ifactor$SO,]
journal.mat.ord <- group.journal[tolower(rownames(group.journal)) %in% ifactor$SO,]

ifactor.ord <- ifactor$ifactor[match(rownames(journal.prop.mat.ord),
                                     ifactor$SO)]

# run ordination to compare journals based on group membership
j.group.ord<-metaMDS(journal.mat.ord, distance="jaccard", trymax=100)

ord.ifactor<-ordisurf(j.group.ord ~ ifactor.ord, plot=FALSE, isotropic=FALSE,
                      knots=4, bs="cr", fx=TRUE, select=FALSE)

pred.x<-seq(min(j.group.ord$points[,1]-0.025), 
            max(j.group.ord$points[,1]+0.025), length.out=100)
pred.y<-seq(min(j.group.ord$points[,2]-0.025), 
            max(j.group.ord$points[,2]+0.025), length.out=100)

pred.df<-expand.grid(pred.x, pred.y)

ifactor.preds<-calibrate(ord.ifactor, pred.df)



#               WORD ANALYSIS ####

# What I want to do here, is rather than look at the most-used words in each group
# is to use the words used in each group's most cited papers, by each group.
# Hopefully we'll be able to keep the same plot format, just with better words.

group.combos<-expand.grid(sort(colnames(paper.groups)[1:4]),
                          sort(colnames(paper.groups)[1:4]))
#group.combos<-group.combos[group.combos[,1] != group.combos[,2],]
groups<-t(group.combos)[,1]

word.citation.list<-lapply(as.data.frame(t(group.combos)), function(groups){
  
  print(paste0(groups[2], " citing ", groups[1]))
  
  groups<-as.character(groups)
  # first we need a word x citation count, with, ideally:
  # word, word commonness, citation count
  
  head(cit.df.groups)
  
  # step 1, get a paper citation count by each group
  outgroup.cite<-with(cit.df.groups[cit.df.groups[, paste0(groups[1], ".cited")] &
                                    cit.df.groups[, paste0(groups[2], ".citing")], ],
                    as.data.frame(table(cited)))
  colnames(outgroup.cite)<-c("index", "outgroup.cite")
  
  ingroup.cite<-with(cit.df.groups[cit.df.groups[, paste0(groups[1], ".cited")] &
                                      cit.df.groups[, paste0(groups[1], ".citing")], ],
                      as.data.frame(table(cited)))
  colnames(ingroup.cite)<-c("index", "ingroup.cite")
  
  temp.frame<-merge(frame.ecology[,c("index", "TI", "DE", "PY")], outgroup.cite, 
                    by.x="index", by.y="index",
                    all.x=TRUE, all.y=FALSE)
  
  temp.frame<-merge(temp.frame, ingroup.cite, 
                    by.x="index", by.y="index",
                    all.x=TRUE, all.y=FALSE)
  
  frame<-temp.frame[paper.groups[,groups[1]] &
                    temp.frame$PY >=2007,]
  frame$outgroup.cite[is.na(frame$outgroup.cite)]=0
  frame$ingroup.cite[is.na(frame$ingroup.cite)]=0

  summary(frame$ingroup.cite>0)
  sum(frame$ingroup.cite)
  
  # get a citation count for each word for each paper
  
  frame$words<-tolower(paste0(frame$TI, " ", ifelse(is.na(frame$DE), "", frame$DE)))
  frame$words<-gsub("[[:punct:]]", " ", frame$words)
  
  listed.words<-strsplit(frame$words, " ")
  stemmed.words<-lapply(listed.words, wordStem)
  
  # quick reference table for back-transforming words
  word.tables<-cbind(listed=unlist(listed.words), 
                     stemmed=unlist(stemmed.words))
  
  unlisted.words<-data.frame(words=unlist(stemmed.words),
                             outgroup.cite=rep(frame$outgroup.cite,
                                               sapply(listed.words, length)),
                             ingroup.cite=rep(frame$ingroup.cite,
                                              sapply(listed.words, length)))
  sum(unlisted.words$ingroup.cite)

  word.agg<-data.frame(word=levels(unlisted.words$word),
                       outgroup.cite=tapply(unlisted.words$outgroup.cite, 
                                            unlisted.words$word, sum),
                       ingroup.cite = tapply(unlisted.words$ingroup.cite,
                                             unlisted.words$word, sum))
  
  word.count<-data.frame(word=as.data.frame(table(unlisted.words$words))$Var1,
                         count=as.data.frame(table(unlisted.words$words))$Freq)
  
  word.agg<-merge(word.agg, word.count, all.x=TRUE, all.y=FALSE,
                  by.x="word", by.y="word")
  
  word.agg$outgroup.weight<-word.agg$outgroup.cite / word.agg$count
  word.agg$ingroup.weight<-word.agg$ingroup.cite / word.agg$count
  
  word.agg<-word.agg[word.agg$count >= floor(sum(paper.groups[,groups[1]])*0.025),]
  
  word.agg<-word.agg[!word.agg$word %in% stop_words$word &
                     word.agg$word != "", ]
  
  word.agg<-word.agg[order(word.agg$ingroup.weight, decreasing=TRUE),]

  # back-transform stemmed words with their most common extension
  word.agg$stemmed.words<-word.agg$word
  word.agg$word <- sapply(word.agg$word, function(x){
    
    temp.table<-table(word.tables[word.tables[,"stemmed"]==x, "listed"])
    
    return(names(temp.table[which.max(temp.table)]))
    
  })
  
  temp.applied.words <- c("restor", "climate", "change", "conserv",
                          "invas")
  
  word.agg<-word.agg[!grepl(paste0(temp.applied.words, collapse="|"),
                            word.agg$word),]
  
  return(word.agg)

})  
names(word.citation.list) = 1:length(word.citation.list)

lapply(word.citation.list, dim)

word.citation.cor<-sapply(word.citation.list, function(x){
  cor(x[,c("outgroup.weight", "ingroup.weight")])[1,2]
})

temp<-cbind(sapply(word.citation.list, function(x){
  cor(x[,c("outgroup.weight", "ingroup.weight")])[1,2]
}),
sapply(word.citation.list, function(x){
  cor(log(x[,c("outgroup.weight", "ingroup.weight")]))[1,2]
}))

word.citation.cormat<-matrix(word.citation.cor, nrow=4, ncol=4)

word.citation.cormat=matrix(sprintf("%.2f", round(word.citation.cormat, 2)),
                            nrow=4, ncol=4)
  
rownames(word.citation.cormat)=sort(colnames(paper.groups)[1:4])
colnames(word.citation.cormat)=sort(colnames(paper.groups)[1:4])

write.csv(word.citation.cormat, 
          "./Outputs/word citation rate correlation matrix.csv")

test <- c(1,5,9,13)

do.call("rbind", lapply(word.citation.list[test], function(x){
  x[1,]
}))


# ####
# FIGURES
# FIGURES ####

colour.rgb<-as.data.frame(col2rgb(c("blue","red","darkgreen","orange","black")))
colnames(colour.rgb)<-c("restoration","climatechange","conservation","invasion","ecology")

#              FIGURE 1 ####
?pdf
pdf(paste0("Staples_et_al_2018-Figure1 ",
           Sys.Date(), ".pdf"), height=2.95, width=3.25, useDingbats=FALSE,
    fonts="Helvetica")

split.screen(rbind(c(0.125,0.535,0.525,0.95),
                   c(0.535,0.945,0.525,0.95),
                   c(0.125,0.535,0.1,0.525),
                   c(0.535,0.945,0.1,0.525)))

# get matrix of colours as rgb
colour.rgb<-as.data.frame(col2rgb(c("red","darkgreen", "black", "orange","blue")))
colnames(colour.rgb)<-sort(colnames(paper.groups))

#                               OVERALL PUB RATE ####
screen(1)

# predict CIs around splines
CIs<-spline.CIs(model=total.pub.rate.bGAM,
                num.data=pub.rate$PY,
                spline.point.n=200)

raw.counts<-sapply(1990:2017, function(x){
  sum(rowSums(paper.groups[frame.ecology$PY == x ,1:4]>0))
})

prop.rates <- (raw.counts[-1] / raw.counts[-length(raw.counts)]) - 1
mean(prop.rates)

par(mar=c(0,0,0,0), ps=6, tck=-0.025, mgp=c(3,0.3,0), las=1)
plot(y=NULL, x=NULL, type="n",
     ylim=c(-0.0075,0.275), xlim=c(1989.25,2022), yaxs="i", xaxs="i", xlab="", ylab="", axes=FALSE)
axis(side=1, at=1990:2017, tck=-0.01, labels=NA)
axis(side=1, at=c(1990,2000,2010,2017), mgp=c(3,0.15,0), labels=NA)
 axis(side=2, at=seq(0,0.3,0.1), labels=c("0.00", "0.10", "0.20", "0.30"))
 axis(side=2, at=seq(0.05,0.25,0.1), tck=-0.01, labels=NA)
mtext(side=2, line=1.25, at=par("usr")[3],
      text="Proportion of ecology", las=3)

# raw points
points(y=pub.rate$pub.rate,
       x=pub.rate$PY, pch=16, cex=0.5, lwd=0.5,
       col=rgb(0.5,0.5,0.5,0.4))

# plot GAM slopes
polygon(y=c(plogis(CIs$fit + 1.96*CIs$se.fit), 
            rev(plogis(CIs$fit - 1.96*CIs$se.fit))),
        x=c(CIs$PY, rev(CIs$PY)),
        border=NA, col=rgb(0,0,0,0.2))

# do falses across all range
points(y=plogis(CIs$fit),
       x=CIs$PY,
       type="l", lwd=1.5,
       col="black")

box()
text(x=relative.axis.point(0.075, "x"),
     y=relative.axis.point(0.925, "y"),
     labels="(a)", font=2, cex=1)

mtext(side=3, line=0, text="New publications", font=2)

close.screen(1)


#                               OVERALL CIT RATE ####

screen(2)

# predict CIs around splines
CIs<-spline.CIs(model=overall.cit.rate.bGAM,
                num.data=cit.rate$PY,
                spline.point.n=200)

par(mar=c(0,0,0,0), ps=6, tck=-0.025, mgp=c(3,0.3,0), las=1)
plot(y=NULL, x=NULL, type="n",
     ylim=c(-0.0075,0.275), xlim=c(1989.25,2022), yaxs="i", xaxs="i", xlab="", ylab="", axes=FALSE)
axis(side=1, at=1990:2017, tck=-0.01, labels=NA)
axis(side=1, at=c(1990,2000,2010,2017), mgp=c(3,0.15,0), labels=NA)

axis(side=2, at=seq(0,0.3,0.1), labels=NA)
axis(side=2, at=seq(0.05,0.25,0.1), tck=-0.01, labels=NA)

# raw points
points(y=cit.rate$cit.rate,
       x=cit.rate$PY, pch=16, cex=0.5, lwd=0.5,
       col=rgb(0.5,0.5,0.5,0.4))


# significant of GAM slopes

# plot GAM slopes
polygon(y=c(plogis(CIs$fit + 1.96*CIs$se.fit), 
            rev(plogis(CIs$fit - 1.96*CIs$se.fit))),
        x=c(CIs$PY, rev(CIs$PY)),
        border=NA, col=rgb(0,0,0,0.2))

# do falses across all range
points(y=plogis(CIs$fit),
       x=CIs$PY,
       type="l", lwd=1.5,
       col="black")

box()
text(x=relative.axis.point(0.075, "x"),
     y=relative.axis.point(0.925, "y"),
     labels="(b)", font=2, cex=1)


mtext(side=3, line=0, text="Citations received", font=2, las=0)

mtext(side=4, line=-0.2, text="All subfields combined", font=2, las=0)

close.screen(2)
#                               GROUP PUB RATE ####
screen(3)

# predict CIs around splines
group.CIs<-grouped.spline.CIs(model=group.pub.rate.bGAM,
                        num.data=group.pub.rate$PY,
                        group.data=group.pub.rate$group,
                        spline.point.n=200)

par(mar=c(0,0,0,0), ps=6, tck=-0.025, mgp=c(3,0.3,0), las=1)
plot(y=NULL, x=NULL, type="n",
     ylim=c(-0.002,0.10), xlim=c(1989.25,2022), yaxs="i", xaxs="i", xlab="", ylab="", axes=FALSE)
axis(side=1, at=1990:2017, tck=-0.01, labels=NA)
axis(side=1, at=c(1990,2000,2010,2017), mgp=c(3,0.15,0), labels=NA)
axis(side=1, at=c(seq(1990,2010, 10), 2017), mgp=c(3,-0.1,0))
axis(side=1, at=c(2017), mgp=c(3,-0.1,0))

axis(side=2, at=seq(0,0.08,0.02))
axis(side=2, at=seq(0.01,0.09,0.02), tck=-0.01, labels=NA)

mapply(x=1:4,
       temp.col=colour.rgb[,levels(group.pub.rate$group)],
       end.names=c("Cl","Co","I","R"),
       CIs=group.CIs,
       function(x, temp.col, end.names, CIs){

         var<-levels(group.pub.rate$group)[x]

         # raw points
         points(y=group.pub.rate$pub.rate[group.pub.rate$group==var],
                x=group.pub.rate$PY[group.pub.rate$group==var], 
                pch=16, cex=0.5, lwd=0.5,
                col=rgb(temp.col[1]/255,
                        temp.col[2]/255,
                        temp.col[3]/255,
                        0.4))

         # plot GAM slopes
         polygon(y=c(plogis(CIs$fit + 1.96*CIs$se.fit), 
                     rev(plogis(CIs$fit - 1.96*CIs$se.fit))),
                 x=c(CIs$PY, rev(CIs$PY)),
                 border=NA,   col=rgb(temp.col[1]/255,
                                      temp.col[2]/255,
                                      temp.col[3]/255,
                                      0.2))
         
         # do falses across all range
         points(y=plogis(CIs$fit),
                x=CIs$PY,
                type="l", lwd=1.5,
                col=rgb(temp.col[1]/255,
                        temp.col[2]/255,
                        temp.col[3]/255,
                        1))
         
         #  right-axis labels
         text(x=2018, y=plogis(CIs$fit[200]),
              adj=0, labels=end.names, font=2,
              col=rgb(temp.col[1]/255,
                      temp.col[2]/255,
                      temp.col[3]/255,
                      1))

       })
box()
text(x=relative.axis.point(0.075, "x"),
     y=relative.axis.point(0.925, "y"),
     labels="(c)", font=2, cex=1)

mtext(side=1, at=par("usr")[2], line=0.5, text="Year")

close.screen(3)

#                               CIT RATE ####
screen(4)

# predict CIs around splines
group.CIs<-grouped.spline.CIs(model=group.cit.rate.bGAM,
                        num.data=cit.summ$PY,
                        group.data=cit.summ$group,
                        spline.point.n=200)

par(mar=c(0,0,0,0), ps=6, tck=-0.025, mgp=c(3,0.3,0), las=1)
plot(y=NULL, x=NULL, type="n",
     ylim=c(-0.002,0.10), xlim=c(1989.25,2022), yaxs="i", xaxs="i", xlab="", ylab="", axes=FALSE)
axis(side=1, at=1990:2017, tck=-0.01, labels=NA)
axis(side=1, at=c(1990,2000,2010,2017), mgp=c(3,0.15,0), labels=NA)
axis(side=1, at=c(seq(1990,2010, 10), 2017), mgp=c(3,-0.1,0))
axis(side=1, at=c(2017), mgp=c(3,-0.1,0))


axis(side=2, at=seq(0,0.08,0.02), labels=NA)
axis(side=2, at=seq(0.01,0.09,0.02), tck=-0.01, labels=NA)

mapply(x=1:4,
       temp.col=colour.rgb[,levels(group.pub.rate$group)],
       end.names=c("Cl","Co","I","R"),
       CIs=group.CIs,
       function(x, temp.col, end.names, CIs){

         var<-levels(cit.summ$group)[x]

         if(var=="ecology"){return(NULL)}

         # raw points
         points(y=cit.summ$cit.prop[cit.summ$group==var],
                x=cit.summ$PY[cit.summ$group==var],
                pch=16, cex=0.5, lwd=0.5,
                col=rgb(temp.col[1]/255,
                        temp.col[2]/255,
                        temp.col[3]/255,
                        0.4))

         # plot GAM slopes
         polygon(y=c(plogis(CIs$fit + 1.96*CIs$se.fit), 
                     rev(plogis(CIs$fit - 1.96*CIs$se.fit))),
                 x=c(CIs$PY, rev(CIs$PY)),
                 border=NA,   col=rgb(temp.col[1]/255,
                                      temp.col[2]/255,
                                      temp.col[3]/255,
                                      0.2))
         
         # do falses across all range
         points(y=plogis(CIs$fit),
                x=CIs$PY,
                type="l", lwd=1.5,
                col=rgb(temp.col[1]/255,
                        temp.col[2]/255,
                        temp.col[3]/255,
                        1))
         
                  #  right-axis labels
         text(x=2018, y=plogis(CIs$fit[200]),
              adj=0, labels=end.names, font=2,
              col=rgb(temp.col[1]/255,
                      temp.col[2]/255,
                      temp.col[3]/255,
                      1))
         
       })
box()
text(x=relative.axis.point(0.075, "x"),
     y=relative.axis.point(0.925, "y"),
     labels="(d)", font=2, cex=1)

mtext(side=4, line=-0.2, text="Subfields separate", font=2, las=0)

close.screen(4)
close.screen(all.screens=TRUE)
dev.off()























#             FIGURE 2 ####

pdf(paste0("Staples_et_al_2018-Figure2-1 ",
           Sys.Date(), ".pdf"), height=3.05, width=6)

split.screen(rbind(c(0.1,0.32,0.55,0.95),
                   c(0.32,0.54,0.55,0.95),
                   c(0.54,0.76,0.55,0.95),
                   c(0.76,0.98,0.55,0.95),

                   c(0.1,0.32,0.15,0.55),
                   c(0.32,0.54,0.15,0.55),
                   c(0.54,0.76,0.15,0.55),
                   c(0.76,0.98,0.15,0.55)))

colour.rgb<-as.data.frame(col2rgb(c("red","darkgreen", "black", "orange","blue")))
colnames(colour.rgb)<-sort(colnames(paper.groups))
colour.rgb<-colour.rgb[,match(levels(bween.group.citrate.binary$citing.group), 
                              colnames(colour.rgb))]
end.names<-c("R","Co","Cl","I")
names(end.names)=levels(bween.group.citrate.binary$citing.group)

CIs<-grouped.spline.CIs(model=bween.citrate.bin.gam,
                        num.data=bween.group.citrate.binary$PY,
                        group.data=bween.group.citrate.binary$group,
                        spline.point.n=200)

names(CIs)<-substr(levels(bween.group.citrate.binary$group),
                   1,
                   regexpr(":", levels(bween.group.citrate.binary$group))-1)

mapply(group=unique(names(CIs)),
       letter=letters[1:4],
       screen=1:4,
       full.names=c("Climate change", "Conservation", "Invasion", "Restoration"),
       function(group, letter, screen, full.names){

  screen(screen)
  par(mar=c(0,0,0,0), ps=6, tck=-0.025, mgp=c(3,0.3,0), las=1)
  plot(y=NULL, x=NULL, type="n",
       ylim=c(-0.025,1.1), xlim=c(1989.25,2022), yaxs="i", xaxs="i",
       xlab="", ylab="", axes=FALSE)
  sub.pred<-CIs[names(CIs)==group]

  lapply(sub.pred, function(group.preds){
  
  group.preds$main.group<-substr(group.preds$group, 
                                 1, 
                                 regexpr(":", group.preds$group)-1)
  
  group.preds$spline.groups<-substr(group.preds$group, 
                                    regexpr(":", group.preds$group)+1,
                                    nchar(as.character(group.preds$group)))
  
  temp.col<-colour.rgb[,group.preds$spline.groups[1]]
  
  
  with(bween.group.citrate.binary[bween.group.citrate.binary$group==group.preds$group[1],],
    points(y=cite.1/(cite.1+cite.0),
           x=PY, pch=16, cex=0.5, lwd=0.5,
           col=rgb(temp.col[1]/255,
                   temp.col[2]/255,
                   temp.col[3]/255,
                   0.4)))
  
  group.preds$upper<-plogis(group.preds$fit + 1.96*group.preds$se.fit)
  group.preds$lower<-plogis(group.preds$fit - 1.96*group.preds$se.fit)
  group.preds$fit <- plogis(group.preds$fit)
  
  # do falses across all range
  points(y=group.preds$fit,
         x=group.preds$PY,
         type="l", lwd=1.5,
         col=rgb(temp.col[1]/255,
                 temp.col[2]/255,
                 temp.col[3]/255,
                 1))
  
  polygon(y=c(group.preds$upper, rev(group.preds$lower)),
          x=c(group.preds$PY, rev(group.preds$PY)),
         col=rgb(temp.col[1]/255,
                 temp.col[2]/255,
                 temp.col[3]/255,
                 0.2), border=NA)
    
  
  #  right-axis labels
  text(x=2018, y=group.preds$fit[200],
       adj=0, 
       labels=end.names[substr(group.preds$group[1],
                               regexpr(":", group.preds$group[1])+1,
                               nchar(as.character(group.preds$group[1])))],
       font=2,
       col=rgb(temp.col[1]/255,
               temp.col[2]/255,
               temp.col[3]/255,
               1))
  })

  if(screen == 1){
    axis(side=2, at=seq(0,1,0.2))
    axis(side=2, at=seq(0.1,0.9,0.1), tck=-0.01, labels=NA)
    mtext(side=2, line=1.85,
          text=bquote("Probability of citing"), las=3)
    mtext(side=2, line=1.25,
          text=bquote("" >= "one paper from subfield"), las=3)
  } else {
    axis(side=2, at=seq(0,1,0.2), labels=NA)
    axis(side=2, at=seq(0.1,0.9,0.1), tck=-0.01, labels=NA)
    
  }
  
  axis(side=1, at=c(1990,2000,2010,2017), tck=-0.025, labels=NA,
       mgp=c(3,-0.2,0))
  axis(side=1, at=2017, tck=-0.025, labels=NA,
       mgp=c(3,-0.2,0))
  axis(side=1, at=1990:2016, tck=-0.01, labels=NA)
  
  box()
  text(x=relative.axis.point(0.075, "x"),
       y=relative.axis.point(0.925, "y"),
       labels=paste0("(", letter,")"), font=2, cex=1)
  
  text(x=relative.axis.point(0.125, "x"),
       y=relative.axis.point(0.925, "y"),
       labels=full.names, las=3, font=2, adj=0,
       col=rgb(colour.rgb[1,group]/255,
               colour.rgb[2,group]/255,
               colour.rgb[3,group]/255))
  
  
  box()
  close.screen(screen)
  
       })  

CIs<-grouped.spline.CIs(model=bween.citrate.beta.gam,
                        num.data=bween.group.citrate.count.agg$PY,
                        group.data=bween.group.citrate.count.agg$group,
                        spline.point.n=200)
names(CIs)<-substr(levels(bween.group.citrate.count.agg$group),
                   1,
                   regexpr(":", levels(bween.group.citrate.count.agg$group))-1)

names(CIs)<-substr(levels(bween.group.citrate.count$group),
                   1,
                   regexpr(":", levels(bween.group.citrate.count$group))-1)

mapply(group=unique(names(CIs)),
       letter=letters[5:8],
       screen=5:8,
       full.names=c("Climate change", "Conservation", "Invasion", "Restoration"),
       function(group, letter, screen, full.names){
         
         screen(screen)
         par(mar=c(0,0,0,0), ps=6, tck=-0.025, mgp=c(3,0.3,0), las=1)
         
         plot(y=NULL, x=NULL, type="n",
              ylim=c(0,0.13), xlim=c(1989.25,2022), yaxs="i", xaxs="i",
              xlab="", ylab="", axes=FALSE)
         
         sub.pred<-CIs[names(CIs)==group]
         
         lapply(sub.pred, function(group.preds){
           
           group.preds$main.group<-substr(group.preds$group, 
                                          1, 
                                          regexpr(":", group.preds$group)-1)
           
           group.preds$spline.groups<-substr(group.preds$group, 
                                             regexpr(":", group.preds$group)+1,
                                             nchar(as.character(group.preds$group)))
           
           temp.col<-colour.rgb[,group.preds$spline.groups[1]]
           
           # with(bween.group.citrate.count.agg[bween.group.citrate.count.agg$group==group.preds$group[1],],
           #      points(y=cite.prop,
           #             x=sort(unique(PY)), pch=16, cex=0.5, lwd=0.5,
           #             col=rgb(temp.col[1]/255,
           #                     temp.col[2]/255,
           #                     temp.col[3]/255,
           #                     0.4)))
          # plot GAM slopes
           
           group.preds$upper<-plogis(group.preds$fit + 1.96*group.preds$se.fit)
           group.preds$lower<-plogis(group.preds$fit - 1.96*group.preds$se.fit)
           group.preds$fit <- plogis(group.preds$fit)
           
           # do falses across all range
           points(y=group.preds$fit,
                  x=group.preds$PY,
                  type="l", lwd=1.5,
                  col=rgb(temp.col[1]/255,
                          temp.col[2]/255,
                          temp.col[3]/255,
                          1))
           
           polygon(y=c(group.preds$upper, rev(group.preds$lower)),
                   x=c(group.preds$PY, rev(group.preds$PY)),
                   col=rgb(temp.col[1]/255,
                           temp.col[2]/255,
                           temp.col[3]/255,
                           0.2), border=NA)
           
           
          #  right-axis labels
           text(x=2018, y=group.preds$fit[200],
                adj=0, 
                labels=end.names[substr(group.preds$group[1],
                                        regexpr(":", group.preds$group[1])+1,
                                        nchar(as.character(group.preds$group[1])))],
                font=2,
                col=rgb(temp.col[1]/255,
                        temp.col[2]/255,
                        temp.col[3]/255,
                        1))
         })
         
         if(screen == 5){
           axis(side=2, at=seq(0,0.3,0.05))
           axis(side=2, at=seq(0,0.3,0.025), tck=-0.01, labels=NA)
           mtext(side=2, line=1.85,
                 text="Proportion of citations", las=3)
           mtext(side=2, line=1.25,
                 text="made to subfield", las=3)
         } else {
           axis(side=2, at=seq(0,0.3,0.05), labels=NA)
           axis(side=2, at=seq(0,0.3,0.025), tck=-0.01, labels=NA)
           
         }
         
         axis(side=1, at=c(1990,2000,2010,2017), tck=-0.025,
              mgp=c(3,-0.2,0))
         axis(side=1, at=1990:2016, tck=-0.01, labels=NA)
         axis(side=1, at=2017, tck=-0.025,
              mgp=c(3,-0.2,0))
         
         if(screen == 7){
           mtext(side=1, text="Year", line=0.5, at=par("usr")[1])
         }
         
         
         box()
         text(x=relative.axis.point(0.075, "x"),
              y=relative.axis.point(0.925, "y"),
              labels=paste0("(", letter,")"), font=2, cex=1)
         
         text(x=relative.axis.point(0.125, "x"),
              y=relative.axis.point(0.925, "y"),
              labels=full.names, las=3, font=2, adj=0,
              col=rgb(colour.rgb[1,group]/255,
                      colour.rgb[2,group]/255,
                      colour.rgb[3,group]/255))
         
         
         box()
         close.screen(screen)
         
       })  


close.screen(all.screens=TRUE)
dev.off()

#                             CROSS-FIELD PAPERS ####

mapply(model=cross.field.bGAM,
       data=cross.field.rates,
       letter=letters[1:4],
       screen=1:4,
       full.names=c("Climate change", "Conservation", "Invasion", "Restoration"),
       function(model, data, letter, screen, full.names){

         screen(screen)

         CIs<-grouped.spline.CIs(model=model,
                                 num.data=data$PY,
                                 group.data=data$group,
                                 spline.point.n=200)

         par(mar=c(0,0,0,0), ps=6, tck=-0.025, mgp=c(3,0.3,0), las=1)
         plot(y=NULL, x=NULL, type="n",
              ylim=c(-0.005,0.19), xlim=c(1989.25,2021), yaxs="i", xaxs="i",
              xlab="", ylab="", axes=FALSE)

         axis(side=1, at=1990:2016, tck=-0.01, labels=NA)
         axis(side=1, at=c(1990,2000,2010,2016), mgp=c(3,0.15,0), labels=NA)

         if(screen == 1){
           axis(side=2, at=seq(0,0.2,0.05))
           axis(side=2, at=seq(0,0.2,0.025), tck=-0.01, labels=NA)
           mtext(side=2, line=1.5,
                 text="Proportional of citations made to subfield", las=3)
         } else {
           axis(side=2, at=seq(0,0.2,0.05), labels=NA)
           axis(side=2, at=seq(0,0.2,0.025), tck=-0.01, labels=NA)
           axis(side=2, at=seq(0.1,0.9,0.1), tck=-0.01, labels=NA)
         }

         mapply(x=1:3,
                temp.col=colour.rgb[,levels(data$group)],
                end.names=end.names[levels(data$group)],
                tempcis=CIs,
                function(x, temp.col, end.names, tempcis){

                  var<-levels(data$group)[x]

                  # raw points
                  points(y=data$rate[data$group==var],
                         x=1990:2016, pch=16, cex=0.5, lwd=0.5,
                         col=rgb(temp.col[1]/255,
                                 temp.col[2]/255,
                                 temp.col[3]/255,
                                 0.4))

                  # significant of GAM slopes

                  # break up spline into significant and non-significant sections
                  cisave<-tempcis
                  groups.df<-cbind(0,0)
                  index=1
                  repeat{

                    temp.start<-min(which(tempcis$diff0==TRUE))
                    if(temp.start == Inf) {break}
                    temp.end<-min(which(tempcis$diff0[temp.start:dim(cisave)[1]]==FALSE))

                    temp.start<- temp.start + groups.df[index,2]
                    temp.end<- ifelse(temp.end == Inf | temp.end>dim(cisave)[1],
                                      dim(cisave)[1],
                                      temp.end + temp.start + groups.df[index,2])

                    groups.df<-rbind(groups.df,
                                     c(temp.start, temp.end))
                    if(dim(cisave)[1]==temp.end){break}

                    tempcis<-tempcis[(temp.end+1):dim(cisave)[1],]
                    index = index+1

                  }

                  # plot GAM slopes

                  # do falses across all range
                  points(y=cisave$fit,
                         x=cisave$PY,
                         type="l", lwd=1.5, lty="21",
                         col=rgb(temp.col[1]/255,
                                 temp.col[2]/255,
                                 temp.col[3]/255,
                                 1))

                  # and do TRUEs for each specific region
                  sapply(2:dim(groups.df)[1], function(n){

                    y<-groups.df[n,]

                    points(y=cisave$fit[y[1]:y[2]],
                           x=cisave$PY[y[1]:y[2]],
                           type="l", lwd=1.5,
                           col=rgb(temp.col[1]/255,
                                   temp.col[2]/255,
                                   temp.col[3]/255,
                                   1))
                  })

                  #  right-axis labels
                  text(x=2017, y=cisave$fit[200],
                       adj=0, labels=end.names, font=2,
                       col=rgb(temp.col[1]/255,
                               temp.col[2]/255,
                               temp.col[3]/255,
                               1))

                })

         box()
         text(x=relative.axis.point(0.075, "x"),
              y=relative.axis.point(0.925, "y"),
              labels=paste0("(", letter,")"), font=2, cex=1)

         text(x=relative.axis.point(0.125, "x"),
              y=relative.axis.point(0.925, "y"),
              labels=full.names, las=3, font=2, adj=0,
              col=rgb(colour.rgb[1,screen]/255,
                      colour.rgb[2,screen]/255,
                      colour.rgb[3,screen]/255))

         close.screen(screen)
       })

#                             BETWEEN SUBFIELD CITING RATES ####

mapply(model=bween.group.citrate.bGAM,
       data=bween.group.citrate,
       letter=letters[5:8],
       screen=5:8,
       full.names=c("Climate change", "Conservation", "Invasion", "Restoration"),
       function(model, data, letter, screen, full.names){

         screen(screen)

         CIs<-grouped.spline.CIs(model=model,
                                 num.data=data$PY,
                                 group.data=data$group,
                                 spline.point.n=200)

         par(mar=c(0,0,0,0), ps=6, tck=-0.025, mgp=c(3,0.3,0), las=1)
         plot(y=NULL, x=NULL, type="n",
              ylim=c(-0.025,1.1), xlim=c(1989.25,2021), yaxs="i", xaxs="i",
              xlab="", ylab="", axes=FALSE)
         axis(side=1, at=1990:2016, tck=-0.01, labels=NA)

         if(screen == 5){
           axis(side=2, at=seq(0,1,0.2))
           axis(side=2, at=seq(0.1,0.9,0.1), tck=-0.01, labels=NA)
           mtext(side=2, line=1.5,
                 text="Proportional citing rate", las=3)
         } else {
           axis(side=2, at=seq(0,1,0.2), labels=NA)
           axis(side=2, at=seq(0.1,0.9,0.1), tck=-0.01, labels=NA)

         }

         axis(side=1, at=c(1990,2000,2010,2016), tck=-0.025,
              mgp=c(3,-0.2,0))
         axis(side=1, at=2016, tck=-0.025,
              mgp=c(3,-0.2,0))

         if(screen == 7){
           mtext(side=1, text="Year", line=0.5, at=par("usr")[1])
         }

         mapply(x=1:4,
                temp.col=colour.rgb[,levels(data$group)],
                end.names=c("R","Co","Cl","I"),
                tempcis=CIs,
                function(x, temp.col, end.names, tempcis){

                  var<-levels(data$group)[x]

                  # raw points
                  points(y=data$cit.prop[data$group==var],
                         x=1990:2016, pch=16, cex=0.5, lwd=0.5,
                         col=rgb(temp.col[1]/255,
                                 temp.col[2]/255,
                                 temp.col[3]/255,
                                 0.4))

                  # significant of GAM slopes

                  # break up spline into significant and non-significant sections
                  cisave<-tempcis
                  groups.df<-cbind(0,0)
                  index=1
                  repeat{

                    temp.start<-min(which(tempcis$diff0==TRUE))
                    if(temp.start == Inf) {break}
                    temp.end<-min(which(tempcis$diff0[temp.start:dim(cisave)[1]]==FALSE))

                    temp.start<- temp.start + groups.df[index,2]
                    temp.end<- ifelse(temp.end == Inf | temp.end>dim(cisave)[1],
                                      dim(cisave)[1],
                                      temp.end + temp.start + groups.df[index,2])

                    groups.df<-rbind(groups.df,
                                     c(temp.start, temp.end))
                    if(dim(cisave)[1]==temp.end){break}

                    tempcis<-tempcis[(temp.end+1):dim(cisave)[1],]
                    index = index+1

                  }

                  # plot GAM slopes

                  # do falses across all range
                  points(y=cisave$fit,
                         x=cisave$PY,
                         type="l", lwd=1.5, lty="21",
                         col=rgb(temp.col[1]/255,
                                 temp.col[2]/255,
                                 temp.col[3]/255,
                                 1))

                  # and do TRUEs for each specific region
                  sapply(2:dim(groups.df)[1], function(n){

                    y<-groups.df[n,]

                    points(y=cisave$fit[y[1]:y[2]],
                           x=cisave$PY[y[1]:y[2]],
                           type="l", lwd=1.5,
                           col=rgb(temp.col[1]/255,
                                   temp.col[2]/255,
                                   temp.col[3]/255,
                                   1))
                  })

                  #  right-axis labels
                  text(x=2017, y=cisave$fit[200],
                       adj=0, labels=end.names, font=2,
                       col=rgb(temp.col[1]/255,
                               temp.col[2]/255,
                               temp.col[3]/255,
                               1))
                })

         box()
         text(x=relative.axis.point(0.075, "x"),
              y=relative.axis.point(0.925, "y"),
              labels=paste0("(", letter,")"), font=2, cex=1)

         text(x=relative.axis.point(0.125, "x"),
              y=relative.axis.point(0.925, "y"),
              labels=full.names, las=3, font=2, adj=0,
              col=rgb(colour.rgb[1,screen-4]/255,
                      colour.rgb[2,screen-4]/255,
                      colour.rgb[3,screen-4]/255))

         close.screen(screen)
       })

close.screen(all.screens=TRUE)
dev.off()



#             FIGURE 3 ####

colour.rgb<-as.data.frame(col2rgb(c("red","darkgreen", "black", "orange","blue")))
colnames(colour.rgb)<-sort(colnames(paper.groups))
x<-9
pdf(paste0("./Plots/group word comparison weighted citations 2010 onwards",
           Sys.Date(), ".pdf"), height=5.7, width=6, useDingbats=FALSE)

par(mfrow=c(4,4), mar=c(0,0,0,0), oma=c(4,4,1,1), 
    las=1, mgp=c(3,0.5,0), ps=8, tck=-0.025)

tri.grad<-function(xlimits, ylimits, rgba, alpha, position, grid.size) {
  
  temp<-as.matrix(seq(1,0, length=grid.size)) %*%
    t(as.matrix(seq(1,0, length=grid.size)))
  
  if(position=="upper"){
    temp.rev<-do.call("cbind", lapply(grid.size:1, function(x){temp[x,]}))
    temp.rev[lower.tri(temp.rev)]<-0
  }
  
  if(position=="lower"){
    temp.rev<-do.call("cbind", lapply(1:grid.size, function(x){rev(temp[x,])}))
    temp.rev[upper.tri(temp.rev)]<-0
  }
  
  image(list(x=seq(xlimits[1], xlimits[2], length=grid.size),
             y=seq(ylimits[1], ylimits[2], length=grid.size),
             z=temp.rev),
        col=colorRampPalette(c(rgb(1,1,1,0),
                               rgb(1,1,1,0),
                               rgb(1,1,1,0),
                               rgb(rgba[1],rgba[2],rgba[3],0.25*alpha),
                               rgb(rgba[1],rgba[2],rgba[3],0.5*alpha),
                               rgb(rgba[1],rgba[2],rgba[3],0.75*alpha),
                               rgb(rgba[1],rgba[2],rgba[3],alpha)),
                             alpha=TRUE, bias=1.75)(grid.size),
        add=TRUE, useRaster=TRUE)
  
  return(temp.rev)
}

tri.grad.diag<-function(xlimits, ylimits, rgba, alpha, grid.size) {
  
  temp.rev<-as.matrix(seq(0,1, length=grid.size)) %*%
    t(as.matrix(seq(0,1, length=grid.size)))
  
  image(list(x=seq(xlimits[1], xlimits[2], length=grid.size),
             y=seq(ylimits[1], ylimits[2], length=grid.size),
             z=temp.rev),
        col=colorRampPalette(c(rgb(1,1,1,0),
                               rgb(1,1,1,0),
                               rgb(1,1,1,0),
                               rgb(rgba[1],rgba[2],rgba[3],0.25*alpha),
                               rgb(rgba[1],rgba[2],rgba[3],0.5*alpha),
                               rgb(rgba[1],rgba[2],rgba[3],0.75*alpha),
                               rgb(rgba[1],rgba[2],rgba[3],alpha)),
                             alpha=TRUE, bias=1.75)(grid.size),
        add=TRUE, useRaster=TRUE)
  
  return(temp.rev)
}

upper.half.circle <- function(x,y,r,nsteps=100,...){  
  
  plot.window<-dev.size(units="cm")
  x.axis<-par("usr")[1:2]
  y.axis<-par("usr")[3:4]
  
  # transform each point into cm scores, using the origin of the plot as our
  # 0cm point.
  x<-(x-x.axis[1]) / (x.axis[2]-x.axis[1]) * plot.window[1]
  y<-(y-y.axis[1]) / (y.axis[2]-y.axis[1]) * plot.window[2]
  
  rs <- seq(0,pi,len=nsteps) 
  xc <- x+r*cos(rs) 
  yc <- y+r*sin(rs) 
  
  # now we need to convert our points along the arc back into plot units
  xc<-(xc / plot.window[1]) * (x.axis[2]-x.axis[1]) + x.axis[1]
  yc<-(yc / plot.window[2]) * (y.axis[2]-y.axis[1]) + y.axis[1]
  
  polygon(xc,yc,...) 
} 

lower.half.circle <- function(x,y,r,nsteps=100,...){ 
  
  plot.window<-dev.size(units="cm")
  x.axis<-par("usr")[1:2]
  y.axis<-par("usr")[3:4]
  
  # transform each point into cm scores, using the origin of the plot as our
  # 0cm point.
  x<-(x-x.axis[1]) / (x.axis[2]-x.axis[1]) * plot.window[1]
  y<-(y-y.axis[1]) / (y.axis[2]-y.axis[1]) * plot.window[2]
  
  
  rs <- seq(0,pi,len=nsteps) 
  xc <- x-r*cos(rs) 
  yc <- y-r*sin(rs) 
  
  # now we need to convert our points along the arc back into plot units
  xc<-(xc / plot.window[1]) * (x.axis[2]-x.axis[1]) + x.axis[1]
  yc<-(yc / plot.window[2]) * (y.axis[2]-y.axis[1]) + y.axis[1]
  
  polygon(xc,yc,...) 
} 

deg2rad <- function(deg) {(deg * pi) / (180)}

name.vect<-c(climatechange="Climate change",
             conservation="Conservation",
             invasion="Invasion",
             restoration="Restoration")

xlimits<-c(-2.4,3.5)
ylimits<-c(-3.4,2.5)

exp(xlimits)
exp(ylimits)

lapply(1:16, function(x){
  
  if(!x %in% names(word.citation.list)){
    plot.new()
    return(NULL)
  }
  
  plot(x=NULL, y=NULL, xlim=xlimits, ylim=ylimits, axes=FALSE, xlab="", ylab="")
  box()
  
  index<-which(names(word.citation.list) == x)
  word.list<-word.citation.list[[index]]
  word.groups<-as.character(unlist(group.combos[index,]))
  
  if(x %in% c(1,6,11,16)){
    tri.grad.diag(xlimits=par("usr")[1:2], ylimits=par("usr")[3:4],
                  rgba=colour.rgb[,word.groups[1]]/255,
                  alpha=0.6, grid.size=500)
  } else {
    
    # gradients
    tri.grad(xlimits=par("usr")[1:2], ylimits=par("usr")[3:4],
             rgba=colour.rgb[,word.groups[1]]/255,
             alpha=0.6,
             position="lower", grid.size=500)
    
    tri.grad(xlimits=par("usr")[1:2], ylimits=par("usr")[3:4],
             rgba=colour.rgb[,word.groups[2]]/255,
             alpha=0.6,
             position="upper", grid.size=500)
  }
  
  # Axes & axis labels #
  
  large.ticks<-data.frame(ticks=log(c(0.01,0.1,1,10)),
                          labels=c(0.01,0.1,1,10))
  small.ticks<-data.frame(ticks=log(c(seq(0.01, 0.1, 0.01), seq(0.1,1,0.1), 
                                      seq(1,10,1), seq(10,100,10))),
                          labels=NA)
  
  axis(side=2, at=small.ticks$ticks,
       labels=NA, tck=-0.01)
  axis(side=1, at=small.ticks$ticks,
       labels=NA, tck=-0.01)
  
  if(x %in% c(1,5,9,13)){
    axis(side=2, at=large.ticks$ticks,
         labels=large.ticks$labels)
    
    mtext(side=2, line=1.5, text=name.vect[as.character(group.combos[as.character(x),2])],
          col=rgb(colour.rgb[1,as.character(group.combos[as.character(x),2])]/255,
                  colour.rgb[2,as.character(group.combos[as.character(x),2])]/255,
                  colour.rgb[3,as.character(group.combos[as.character(x),2])]/255),
          cex=1, font=1, las=0)
  }  else {axis(side=2, at=large.ticks$ticks, labels=NA)}
  
  if(x %in% c(13:16)){
    axis(side=1, at=large.ticks$ticks,
         labels=large.ticks$labels, 
         mgp=c(3,0,0))
    
    mtext(side=1, line=1, text=name.vect[as.character(group.combos[as.character(x),1])],
          col=rgb(colour.rgb[1,as.character(group.combos[as.character(x),1])]/255,
                  colour.rgb[2,as.character(group.combos[as.character(x),1])]/255,
                  colour.rgb[3,as.character(group.combos[as.character(x),1])]/255),
          cex=1, font=1)
  }
  
  if(x %in% c(5:16)){
    axis(side=3, at=large.ticks$ticks, labels=NA, tck=0.025)
    
    axis(side=3, at=small.ticks$ticks,
         labels=NA, tck=0.01)
  }
  
  if(x==5){mtext(side=2, 
                 at=par("usr")[3],
                 line=3, las=0,
                 text="Between-subfield weighted citation rate (citation count / word commonness)")
  }
  
  if(x==14){mtext(side=1,
                  at=par("usr")[2],
                  line=2.25,
                  text="Within-subfield weighted citation rate (citation count / word commonness)")
  }
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.915, "y"),
       labels=substitute(bold("("*a*")")*": "*plain("n = "*b),
                         list(a = letters[which(names(word.citation.list) == x)],
                              b=dim(word.list)[1])), 
       adj=0, cex=1.25)
  
  # Points #
  
  points(log(word.list$outgroup.weight) ~ log(word.list$ingroup.weight),
         pch=16, col="grey60")
  
  top5<-word.list[order(word.list$outgroup.weight, decreasing=TRUE)[1:10],]
  right5<-word.list[order(word.list$ingroup.weight, decreasing=TRUE)[1:10],]
  
  # Words
  
  right5<-right5[order(log(right5[,"outgroup.weight"]) + log(right5[,"ingroup.weight"])), ]
  
  # separate joint words to plot separately
  joint.words<-top5[top5$word %in% right5$word,]
  
  right5<- right5[!right5$word %in% top5$word, ]
  
  upper.angles<-seq(200, 50, length=dim(top5)[1])
  
  radii<-c(relative.axis.point(0.80, "x"),
           relative.axis.point(0.80, "y"))
  
  upper.pos<-data.frame(x=mean(log(top5$ingroup.weight)) + radii[1]*cos(deg2rad(upper.angles)),
                        y=mean(log(top5$outgroup.weight)) +
                          radii[2]*sin(deg2rad(upper.angles)))
  
  top5<-top5[order(log(top5[,"outgroup.weight"]) + log(top5[,"ingroup.weight"])), ]
  
  # the lower go from 270 to 360 (or 0)
  lower.angles<-seq(230, 370, length=dim(right5)[1])
  
  lower.angles<-ifelse(lower.angles >360, lower.angles-360, lower.angles)
  
  lower.pos<-data.frame(x=mean(log(right5$ingroup.weight)) +
                          radii[1]*cos(deg2rad(lower.angles)),
                        y=mean(log(right5$outgroup.weight)) +
                          radii[2]*sin(deg2rad(lower.angles)))
  
  upper.offset<-sapply(upper.angles, function(x){
    if(x > 90){return(2)}
    if(x < 70){return(4)}
    return(3)
  })
  
  lower.offset<-sapply(lower.angles, function(x){
    if(x > 280 | x < 50){return(4)}
    if(x < 270){return(2)}
    return(1)
  })
  
  par(xpd=NA)
  text(x=upper.pos[,"x"], y=upper.pos[,"y"],
       labels=top5$word,
       pos=upper.offset,
       offset=0.2,
       cex=0.9)
  
  segments(x0=upper.pos[,"x"], y0=upper.pos[,"y"],
           x1=log(top5[,"ingroup.weight"]), y1=log(top5[,"outgroup.weight"]),
           col="black", lwd=0.5)
  
  if(!x %in% c(1,6,11,16)){
    text(x=lower.pos[,"x"], y=lower.pos[,"y"],
         labels=right5$word, pos=lower.offset,
         offset=0.2, cex=0.9)
    
    segments(x0=lower.pos[,"x"], y0=lower.pos[,"y"],
             x1=log(right5[,"ingroup.weight"]), y1=log(right5[,"outgroup.weight"]),
             col="black", lwd=0.5)
  }
  
  par(xpd=FALSE)
  # coloured points
  semi.radius<-0.3
  
  sapply(1:dim(joint.words)[1], function(y){
    
    upper.half.circle(x=log(joint.words[y,"ingroup.weight"]),
                      y=log(joint.words[y,"outgroup.weight"]),
                      r=semi.radius,
                      col=rgb(colour.rgb[1,as.character(group.combos[as.character(x),2])]/255,
                              colour.rgb[2,as.character(group.combos[as.character(x),2])]/255,
                              colour.rgb[3,as.character(group.combos[as.character(x),2])]/255),
                      border=NA)
    
    lower.half.circle(x=log(joint.words[y,"ingroup.weight"]),
                      y=log(joint.words[y,"outgroup.weight"]),
                      r=semi.radius,
                      col=rgb(colour.rgb[1,as.character(group.combos[as.character(x),1])]/255,
                              colour.rgb[2,as.character(group.combos[as.character(x),1])]/255,
                              colour.rgb[3,as.character(group.combos[as.character(x),1])]/255),
                      border=NA)
    
    points(y=log(joint.words[y,"outgroup.weight"]), x=log(joint.words[y,"ingroup.weight"]),
           pch=21, bg=rgb(0,0,0,0), lwd=0.5)                  
    
  })
  
  with(top5[!top5$word %in% joint.words$word,],
       points(y=log(outgroup.weight), x=log(ingroup.weight), pch=21, lwd=0.5,
              bg=rgb(colour.rgb[1,as.character(group.combos[as.character(x),2])]/255,
                     colour.rgb[2,as.character(group.combos[as.character(x),2])]/255,
                     colour.rgb[3,as.character(group.combos[as.character(x),2])]/255)))
  with(right5,
       points(y=log(outgroup.weight), x=log(ingroup.weight), pch=21, lwd=0.5,
              bg=rgb(colour.rgb[1,as.character(group.combos[as.character(x),1])]/255,
                     colour.rgb[2,as.character(group.combos[as.character(x),1])]/255,
                     colour.rgb[3,as.character(group.combos[as.character(x),1])]/255)))
  box()
})

dev.off()
#             FIGURE 4 ####


options(scipen=1)
pdf(paste0("./Plots/journal ordination ", Sys.Date(), ".pdf"),
    height=4, width=4)
par(ps=8, mar=c(2.5,2.5,2.5,2.5), tck=-0.015, las=1, mgp=c(3,0.5,0))

plot(j.group.ord$points, type="n", axes=FALSE, xlab="", ylab="",
     asp=1)

contour(x=pred.x, y=pred.y,
        z=matrix(ifactor.preds, nrow=100, ncol=100),
        col="grey60", add=TRUE, method="edge", nlevels=6)

size.cats<-cut(rowSums(journal.mat), breaks=c(0,10,100,1000,10000))
table(size.cats)

points(j.group.ord$points, pch=16, 
       col="grey80", 
       cex=c(0.2,0.4,0.6,1)[size.cats])

temp.colour.rgb<-colour.rgb[,match(rownames(j.group.ord$species),
                                   colnames(colour.rgb))]
temp.colour.rgb<-temp.colour.rgb[, order(colnames(temp.colour.rgb))]

arrows(x0=0, y0=0,
       x1=j.group.ord$species[order(rownames(j.group.ord$species)),1],
       y=j.group.ord$species[order(rownames(j.group.ord$species)),2],
       pch=16, col=rgb(temp.colour.rgb[1,]/255,
                       temp.colour.rgb[2,]/255,
                       temp.colour.rgb[3,]/255), lwd=3, length=0.05)

text(j.group.ord$species[order(rownames(j.group.ord$species)),], 
     labels=c("Cl", "Co", "I", "R"),
     col=rgb(temp.colour.rgb[1,]/255,
             temp.colour.rgb[2,]/255,
             temp.colour.rgb[3,]/255), font=2,
     pos=c(3, 3, 1, 3, 3), offset=0.35)

mtext(side=2, line=1.75, text="nMDS2", las=3)
mtext(side=1, line=1.25, text="nMDS1")
axis(side=1, mgp=c(3,0.2,0))
axis(side=2)

# labels for outlying journals
x.quant<-quantile(j.group.ord$points[,1], probs=c(0.025, 0.975))
y.quant<-quantile(j.group.ord$points[,2], probs=c(0.025, 0.975))

outliers<-j.group.ord$points[tolower(rownames(j.group.ord$points)) %in% rownames(journal.sub.mat),]

outlier.cats<-cut(rowSums(journal.mat[match(tolower(rownames(outliers)), 
                                            rownames(journal.mat)),]),
                  breaks=c(0,10,100,1000,10000))

segments(x0=outliers[,1], x1=outliers[,1]+0.1,
         y0=outliers[,2], y1=outliers[,2]+0.1)

points(outliers, pch=16,
       cex=c(0.2,0.4,0.6,1)[outlier.cats])

text(labels=journal.sub.mat$abbrev[match(tolower(rownames(outliers)), rownames(journal.sub.mat))],
     x=outliers[,1] + 0.1, y=outliers[,2] + 0.1, pos=2)

rect(xleft=0.595, xright=par("usr")[2],
     ytop= -0.675, ybottom=par("usr")[3],
     col="white") 

legend(x="bottomright",
       pch=16, pt.cex=c(0.2,0.4,0.6,1),
       col="black",
       legend=c("1 - 10",
                "10 - 100",
                "100 - 1000",
                "1000 +"),
       y.intersp=0.65,
       x.intersp=0.75, bg="white", bty="n")

text(x=relative.axis.point(0.845, "x"),
     y=relative.axis.point(0.23, "y"),
     labels="# papers", font=2, adj=0.5)

box()
dev.off()
# Appendix ####
#       Compare invasion to novel ecosystems

novel.words <- grepl("novel eco|novel comm", paper.words)
summary(novel.words)

novel.pub.rate<-do.call("rbind", lapply(1990:2017, function(x){
  
  temp<-sum(novel.words[frame.ecology$PY==x])
  return(data.frame(pub.rate=temp / sum(frame.ecology$PY==x),
                    group="novel",
                    PY=x))
}))

novel.invas <- droplevels(rbind(novel.pub.rate,
                     group.pub.rate[group.pub.rate$group=="invasion",]))

novel.invas.bGAM<-gam(pub.rate ~ group + s(PY, bs="cr", k=4, by=group),
                         family=betar(),
                         data=novel.invas)
plot(novel.invas.bGAM)

colour.rgb<-as.data.frame(col2rgb(c("blue","red","darkgreen","orange","black")))

colour.rgb <- as.data.frame(col2rgb(c("purple", "orange")))
colnames(colour.rgb)<-c("novel","invasion")


pdf(paste0("Staples_et_al_2018-AppendixFig1 ",
           Sys.Date(), ".pdf"), height=2.95, width=3.25, useDingbats=FALSE)

split.screen(rbind(c(0.125,0.535,0.525,0.95),
                   c(0.535,0.945,0.525,0.95),
                   c(0.125,0.535,0.1,0.525),
                   c(0.535,0.945,0.1,0.525)))


group.CIs<-grouped.spline.CIs(model=novel.invas.bGAM,
                              num.data=novel.invas$PY,
                              group.data=droplevels(novel.invas$group),
                              spline.point.n=200)

screen(1)
par(mar=c(0,0,0,0), ps=6, tck=-0.025, mgp=c(3,0.3,0), las=1)
plot(y=NULL, x=NULL, type="n",
     ylim=c(-0.002,0.10), xlim=c(1989.25,2022), yaxs="i", xaxs="i", xlab="", ylab="", axes=FALSE)
axis(side=1, at=1990:2017, tck=-0.01, labels=NA)
axis(side=1, at=c(1990,2000,2010,2017), mgp=c(3,0.15,0), labels=NA)
axis(side=1, at=c(seq(1990,2010, 10), 2017), mgp=c(3,-0.1,0), labels=NA)
axis(side=1, at=c(2017), mgp=c(3,-0.1,0), labels=NA)

axis(side=2, at=seq(0,0.08,0.02))
axis(side=2, at=seq(0.01,0.09,0.02), tck=-0.01, labels=NA)
mtext(side=2, line=1.5, at=par("usr")[3],
      text="Proportion of ecology", las=3)


x=2
temp.col=colour.rgb[, 2]
end.names=c("I")
CIs=group.CIs[[2]]
         
         var<-levels(novel.invas$group)[x]
         
         # raw points
         points(y=novel.invas$pub.rate[novel.invas$group==var],
                x=novel.invas$PY[novel.invas$group==var], 
                pch=16, cex=0.5, lwd=0.5,
                col=rgb(temp.col[1]/255,
                        temp.col[2]/255,
                        temp.col[3]/255,
                        0.4))
         
         # plot GAM slopes
         polygon(y=c(plogis(CIs$fit + 1.96*CIs$se.fit), 
                     rev(plogis(CIs$fit - 1.96*CIs$se.fit))),
                 x=c(CIs$PY, rev(CIs$PY)),
                 border=NA,   col=rgb(temp.col[1]/255,
                                      temp.col[2]/255,
                                      temp.col[3]/255,
                                      0.2))
         
         # do falses across all range
         points(y=plogis(CIs$fit),
                x=CIs$PY,
                type="l", lwd=1.5,
                col=rgb(temp.col[1]/255,
                        temp.col[2]/255,
                        temp.col[3]/255,
                        1))
         
         #  right-axis labels
         text(x=2018, y=plogis(CIs$fit[200]),
              adj=0, labels=end.names, font=2,
              col=rgb(temp.col[1]/255,
                      temp.col[2]/255,
                      temp.col[3]/255,
                      1))
         
box()
text(x=relative.axis.point(0.075, "x"),
     y=relative.axis.point(0.925, "y"),
     labels="(a)", font=2, cex=1)

close.screen(1)
screen(3)

par(mar=c(0,0,0,0), ps=6, tck=-0.025, mgp=c(3,0.3,0), las=1)
plot(y=NULL, x=NULL, type="n",
     ylim=c(-0.0002,0.01), xlim=c(1989.25,2022), yaxs="i", xaxs="i", xlab="", ylab="", axes=FALSE)
axis(side=1, at=1990:2017, tck=-0.01, labels=NA)
axis(side=1, at=c(1990,2000,2010,2017), mgp=c(3,0.15,0), labels=NA)
axis(side=1, at=c(seq(1990,2010, 10), 2017), mgp=c(3,-0.1,0))
axis(side=1, at=c(2017), mgp=c(3,-0.1,0))

axis(side=2, at=seq(0,0.008,0.002))
axis(side=2, at=seq(0.001,0.009,0.001), tck=-0.01, labels=NA)


x=1
temp.col=colour.rgb[, 1]
end.names=c("N")
CIs=group.CIs[[1]]

var<-levels(novel.invas$group)[x]

# raw points
points(y=novel.invas$pub.rate[novel.invas$group==var],
       x=novel.invas$PY[novel.invas$group==var], 
       pch=16, cex=0.5, lwd=0.5,
       col=rgb(temp.col[1]/255,
               temp.col[2]/255,
               temp.col[3]/255,
               0.4))

# plot GAM slopes
polygon(y=c(plogis(CIs$fit + 1.96*CIs$se.fit), 
            rev(plogis(CIs$fit - 1.96*CIs$se.fit))),
        x=c(CIs$PY, rev(CIs$PY)),
        border=NA,   col=rgb(temp.col[1]/255,
                             temp.col[2]/255,
                             temp.col[3]/255,
                             0.2))

# do falses across all range
points(y=plogis(CIs$fit),
       x=CIs$PY,
       type="l", lwd=1.5,
       col=rgb(temp.col[1]/255,
               temp.col[2]/255,
               temp.col[3]/255,
               1))

#  right-axis labels
text(x=2018, y=plogis(CIs$fit[200]),
     adj=0, labels=end.names, font=2,
     col=rgb(temp.col[1]/255,
             temp.col[2]/255,
             temp.col[3]/255,
             1))
box()
text(x=relative.axis.point(0.075, "x"),
     y=relative.axis.point(0.925, "y"),
     labels="(b)", font=2, cex=1)

mtext(side=1, line=0.5, text="Year")

close.screen(3)
close.screen(all.screens=TRUE)

dev.off()