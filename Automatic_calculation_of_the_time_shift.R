# Cleans everything
rm(list = ls())

# Working directory
setwd("")

# Opening the table containing the results of the script Pairing microphones
library(data.table)
library(Hmisc)
TrajTot=read.csv2("TrajTot.csv")

ListPaires=read.csv2("ListPairs.csv")

maxdecalage=300
decalageFinal=vector()

for (h in 1:nrow(ListPaires)){
  TrajTot1=subset(TrajTot, participation %in% ListPaires$participation_1[h])
  if(nrow(TrajTot1)==0)
  {
    decalageTemp=NA
  }else{
    TrajTot1=cbind(TrajTot1,Cote=1, Point=ListPaires$point_1[h])
    TrajTot2=subset(TrajTot, participation %in% ListPaires$participation_2[h])
    if(nrow(TrajTot2)==0)
    {
      decalageTemp=NA
    }else{
      TrajTot2=cbind(TrajTot2,Cote=2, Point=ListPaires$point_2[h])
      Trav=rbind(TrajTot1,TrajTot2)
      Trav$DecDeb=as.numeric(as.character(Trav$DecDeb))
      Trav$DecFin=as.numeric(as.character(Trav$DecFin))
      Trav$Prob=as.numeric(as.character(Trav$probabilite))
      Trav$espece2=Trav$espece
      Trav$espece=as.numeric(as.factor(Trav$espece))
      Date=vector(length=nrow(Trav))
      Heure=vector(length=nrow(Trav))
      Minute=vector(length=nrow(Trav))
      Seconde=vector(length=nrow(Trav))
      MiliSec=vector(length=nrow(Trav))
      Fich=as.character(Trav$File0)
      for (i in 1:length(Fich)){
        Date[i]=substr(Fich[i],nchar(Fich[i])-18,nchar(Fich[i])-11)
        Heure[i]=substr(Fich[i],nchar(Fich[i])-9,nchar(Fich[i])-8)
        Minute[i]=substr(Fich[i],nchar(Fich[i])-7,nchar(Fich[i])-6)
        Seconde[i]=substr(Fich[i],nchar(Fich[i])-5,nchar(Fich[i])-4)
        MiliSec[i]=substr(Fich[i],nchar(Fich[i])-2,nchar(Fich[i]))
      }
      Datenum=as.numeric(Date)
      Heure=as.numeric(Heure)
      Minute=as.numeric(Minute)
      Seconde=as.numeric(Seconde)
      MiliSec=as.numeric(MiliSec)
      Temps=Heure*3600+Minute*60+Seconde+MiliSec/1000
      Trav=cbind(Trav,Heure,Minute,Seconde,MiliSec,Temps,Date,Datenum)
      Ent1=subset(Trav, DecDeb>0 & DecFin<0 & Cote==1)
      Sort1=subset(Trav, DecDeb<0 & DecFin>0 & Cote==1)
      Ent2=subset(Trav, DecDeb>0 & DecFin<0 & Cote==2)
      Sort2=subset(Trav, DecDeb<0 & DecFin>0 & Cote==2)
      
      if((((nrow(Ent1)==0)|nrow(Sort2)==0))&((nrow(Ent2)==0)|(nrow(Sort1)==0)))
      {
        decalageTemp=0
      }else{
        
        Interval=5
        NbMatch=vector()
        for (j in (-(round(maxdecalage/Interval)):round(maxdecalage/Interval)))
        {
          decalage=j*Interval
          EntTemp=Ent1$Temps+decalage
          SortTemp=Sort1$Temps+decalage
          NbMatchTemp=0
          if((((nrow(Ent1)>0)&(nrow(Sort2)>0))))
          {
            TC12=find.matches(as.matrix(cbind(Ent1$Datenum,Ent1$espece,EntTemp))
                              ,as.matrix(cbind(Sort2$Datenum,Sort2$espece,Sort2$Temps))
                              ,tol=c(0,0,(ListPaires$inter_distance[h]/5+(ListPaires$inter_distance[h]/5*0.5)+5))
                              ,maxmatch=1)
            
            NbMatchTemp=NbMatchTemp+length(subset(TC12$matches,TC12$matches>0))
          }
          if((((nrow(Ent2)>0)&(nrow(Sort1)>0))))
          {
            TC21=find.matches(cbind(Ent2$Datenum,Ent2$espece,Ent2$Temps)
                              ,cbind(Sort1$Datenum,Sort1$espece,SortTemp)
                              ,tol=c(0,0,(ListPaires$inter_distance[h]/5+(ListPaires$inter_distance[h]/5*0.5)+5))
                              ,maxmatch=1)
            NbMatchTemp=NbMatchTemp+length(subset(TC21$matches,TC21$matches>0))
          }
          NbMatch=c(NbMatch,NbMatchTemp) 
          print(paste(h, "/", nrow(ListPaires), " ", j*Interval))
        }
        DecPotentiels=(-(round(maxdecalage/Interval)):round(maxdecalage/Interval))*Interval
        barplot(NbMatch,names.arg=DecPotentiels,cex.names=0.5,main=h)
        DecSelectionnes=subset(DecPotentiels,NbMatch==max(NbMatch))
        decalageTemp=quantile(DecSelectionnes,0.5,type=1)
      }
    }
  }
  decalageFinal=c(decalageFinal,decalageTemp)
}

ListDecalage=cbind(ListPaires,decalageFinal)
write.csv2(ListDecalage,"ListPairs_Dec.csv")
