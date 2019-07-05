reference.code.pages<-function(frame){

# first we need year
year<-as.character(frame$PY)

# now we need journal, which should work as is
journal<-frame$J9

# then we need volume and page, if there are some
volume<-as.character(frame$VL)
table(volume)

page<-frame$BP

page<-gsub("[a-z]|[A-Z]", "", page)
table(page)

refcode<-tolower(paste0(year, ", ",
                        journal, ", ",
                        "v", volume, ", ",
                        "p", page))

return(refcode)

}