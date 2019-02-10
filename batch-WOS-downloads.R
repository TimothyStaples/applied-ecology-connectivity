library(RSelenium)
library(XML)

sapply(list.files(path="/home/timothy/R/x86_64-pc-linux-gnu-library/3.4/RSelenium/examples/serverUtils", 
                  pattern=".R", full.names=TRUE), source)

checkForServer()
startServer()
remDr <- remoteDriver(remoteServerAddr = "localhost" 
                      , port = 4444,
                      browserName="firefox")

remDr$open(silent=TRUE)
remDr$getStatus()

# insert WOS search page here.
# base search

remDr$navigate("http://apps.webofknowledge.com/summary.do?product=WOS&search_mode=AdvancedSearch&page=1&qid=9&SID=D5ugO8cUW3BDM8Z1XLp&parentProduct=WOS&excludeEventConfig=ExcludeIfReload")

# get total number of papers
n.text<-remDr$findElement(using="class", "title4")
n=substr(n.text$getElementText(),
         regexpr(": ", n.text$getElementText())+2,
         nchar(n.text$getElementText()))
n=as.numeric(gsub(",","",n))

# break into groups of 500
start=seq(1,n,500)
end=c((start-1)[-1], n)

options(scipen=999) # make it so 10000 doesn't turn to 1e+5

for(x in 1:length(start)){

  print(x)
  
# navigate to other formats save menu
saveToMenu <-remDr$findElement(using = 'id', "select2-saveToMenu-container")
saveToMenu$clickElement()

Sys.sleep(0.3)
# once open, select "Records"
records<-remDr$findElement(using= "id", "numberOfRecordsRange")
records$clickElement()

# next, add start number to correct box
start.box<-remDr$findElement(using= "id", "markFrom")
start.box$sendKeysToElement(list(as.character(start[x])))
Sys.sleep(0.3)
# add final number to correct box
end.box<-remDr$findElement(using= "id", "markTo")
end.box$sendKeysToElement(list(as.character(end[x])))

Sys.sleep(0.3)
# change field selector to full record and cited references
output.selection.box<-remDr$findElement(using="name", "fields_selection")
output.selection.box$sendKeysToElement(list(key="enter"))
for(i in 1:2){output.selection.box$sendKeysToElement(list(key="down_arrow"))}
output.selection.box$sendKeysToElement(list(key="enter"))

Sys.sleep(0.3)
# change save options to "plain text"
save.options.box<-remDr$findElement(using="name", "save_options")
save.options.box$sendKeysToElement(list(key="enter"))
for(i in 1:3){save.options.box$sendKeysToElement(list(key="down_arrow"))}
save.options.box$sendKeysToElement(list(key="enter"))

# then find the save button and click it
send.buttons<-remDr$findElements(using="class name", "primary-button")
Sys.sleep(0.1)
send.buttons[[2]]$clickElement()

# wait until file has had time to download
Sys.sleep(0.3)
# click close option
close.text<-remDr$findElement(using="link text", "Close")
Sys.sleep(0.1)
close.text$clickElement()

Sys.sleep(1)
}
