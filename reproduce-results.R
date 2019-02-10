# ################################################################### ####
# Title: Applied research is on the rise but connectivity barriers    #### 
#        persist between four major subfields                         ####
# Author: Timothy L Staples                                           ####
# Details: Minimal code to reproduce results and figures              ####  
# ################################################################### ####

# ####
# Working directory ####
rm(list=ls())
setwd("LOCATION OF THIS SCRIPT")
# Libraries ####
library(lme4)
library(vegan)
library(plotrix)
library(mgcv)

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

#                               PUB RATE ####
pub.rate <- read.csv("./Data/overall publication rate data.csv", row.names = 1)
total.pub.rate.bGAM<-gam(pub.rate ~ s(PY, bs="cr", k=4),
                         family=betar(),
                         data=pub.rate,
                         method="REML")
summary(total.pub.rate.bGAM)
gam.check(total.pub.rate.bGAM)

#                               CIT RATE ####

cit.rate <- read.csv("./Data/overall citation rate data.csv", row.names = 1)
overall.cit.rate.bGAM<-gam(cit.rate ~ s(PY, bs="cr", k=4),
                           family=betar(),
                           data=cit.rate,
                           method="REML")
summary(overall.cit.rate.bGAM)
gam.check(overall.cit.rate.bGAM)

#               SUB-FIELD RATES ####
#                               PUB RATE ####

group.pub.rate <- read.csv("./Data/subfield publication rate data.csv", row.names = 1)
group.pub.rate.bGAM<-gam(pub.rate ~ group + s(PY, bs="cr", k=4, by=group),
                         family=betar(),
                         data=group.pub.rate,
                         method="REML")
summary(group.pub.rate.bGAM)
gam.check(group.pub.rate.bGAM)

#                               CIT RATE ####

cit.summ <- read.csv("./Data/subfield citation rate data.csv", row.names = 1)
group.cit.rate.bGAM<-gam(cit.prop ~ group + s(PY, bs="cr", k=4, by=group),
                         family=betar(),
                         data=cit.summ,
                         method="REML")
summary(group.cit.rate.bGAM)
gam.check(group.cit.rate.bGAM)

#               CITATION PROPORTION RATES ####

bween.group.citrate.binary <- read.csv("./Data/binary between group citation data.csv",
                                       row.names=1)

bween.citrate.bin.gam<-gam(cbind(cite.1, cite.0) ~ s(PY, bs="cr", k=4, by=group) + group,
                           data=bween.group.citrate.binary,
                           family=binomial(),
                           method="REML")
summary(bween.citrate.bin.gam)

# beta gam for counts

bween.group.citrate.count <- read.csv("./Data/proportion between group citation data.csv",
                                       row.names=1)

bween.citrate.beta.gam<-gam(cite.prop ~ s(PY, bs="cr", k=4, by=group) + group,
                            data=bween.group.citrate.count,
                            family=betar(),
                            method="REML")

#               JOURNAL ORDINATION ####

group.journal <- read.csv("./Data/journal subfield abundance.csv",
                          row.names=1)
journal.totals <- read.csv("./Data/journal totals.csv", row.names=1)

# replace NAs with 0s
group.journal[is.na(group.journal)]=0
group.journal <- group.journal[rowSums(group.journal)>0,]

# import impact factors
ifactor<-read.csv("./Data/ordination journals impact factor.csv", stringsAsFactors=FALSE)
ifactor<-ifactor[!is.na(ifactor$ifactor), ]

journal.prop.mat<- group.journal / journal.totals[,1]

journal.sub.mat<-journal.prop.mat[rowSums(journal.prop.mat[,1:4]) > 0.4,]

journal.sub.mat <- as.data.frame(journal.sub.mat)
journal.sub.mat$abbrev <- c("AC", "AVS", "AQI", "B&C", "BC", "BI", "CB", "DD", "EMR",
                            "GCB", "JNC", "NAJ", "OR", "RE")

write.csv(journal.sub.mat, "./Outputs/Table 1.csv")

journal.prop.mat<-cbind(journal.prop.mat, journal.totals)

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

word.citation.list <- readRDS("./Data/concept citation rate data.rds")
group.combos<-expand.grid(c("climatechange", "conservation", "invasion", "restoration"),
                          c("climatechange", "conservation", "invasion", "restoration"))

# FIGURES ####

colour.rgb<-as.data.frame(col2rgb(c("red","darkgreen", "black", "orange","blue")))
colnames(colour.rgb)<-c("climatechange", "conservation", "ecology", "invasion", "restoration")

#              FIGURE 1 ####

pdf(paste0("./Plots/Figure1 ",
           Sys.Date(), ".pdf"), height=2.95, width=3.25, useDingbats=FALSE,
    fonts="Helvetica")

split.screen(rbind(c(0.125,0.535,0.525,0.95),
                   c(0.535,0.945,0.525,0.95),
                   c(0.125,0.535,0.1,0.525),
                   c(0.535,0.945,0.1,0.525)))

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

pdf(paste0("./Plots/Figure2 ",
           Sys.Date(), ".pdf"), height=3.05, width=6)

split.screen(rbind(c(0.1,0.32,0.55,0.95),
                   c(0.32,0.54,0.55,0.95),
                   c(0.54,0.76,0.55,0.95),
                   c(0.76,0.98,0.55,0.95),

                   c(0.1,0.32,0.15,0.55),
                   c(0.32,0.54,0.15,0.55),
                   c(0.54,0.76,0.15,0.55),
                   c(0.76,0.98,0.15,0.55)))

colour.rgb<-colour.rgb[,match(levels(bween.group.citrate.binary$citing.group), 
                              colnames(colour.rgb))]
end.names<-c("Cl","Co","I","R")
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
                        num.data=bween.group.citrate.count$PY,
                        group.data=bween.group.citrate.count$group,
                        spline.point.n=200)

names(CIs)<-substr(levels(bween.group.citrate.count$group),
                   1,
                   regexpr(":", levels(bween.group.citrate.count$group))-1)

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
              ylim=c(0,0.19), xlim=c(1989.25,2022), yaxs="i", xaxs="i",
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

#             FIGURE 3 ####

pdf(paste0("./Plots/Figure 3 ",
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
pdf(paste0("./Plots/Figure 4 ", Sys.Date(), ".pdf"),
    height=4, width=4)
par(ps=8, mar=c(2.5,2.5,2.5,2.5), tck=-0.015, las=1, mgp=c(3,0.5,0))

plot(j.group.ord$points, type="n", axes=FALSE, xlab="", ylab="",
     asp=1)

contour(x=pred.x, y=pred.y,
        z=matrix(ifactor.preds, nrow=100, ncol=100),
        col="grey60", add=TRUE, method="edge", nlevels=6)

size.cats<-cut(rowSums(group.journal), breaks=c(0,10,100,1000,10000))
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

outlier.cats<-cut(rowSums(group.journal[match(tolower(rownames(outliers)), 
                                            rownames(group.journal)),]),
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