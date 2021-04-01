# Construct social networks

library(dplyr)
library(igraph)
library(lubridate)

############
# Creating the networks from raw dara-----

sessions <- read.csv("Data/tblSessions.csv", header=TRUE, stringsAsFactors=FALSE, fill=T)
sessions=sessions[sessions$location!="d",] # removing den sessions
sessions=sessions[sessions$location!="n",]
sessions=sessions[sessions$hyenas!="",] # removing empty lines
sessions=sessions[!is.na(sessions$hyenas),]
sessions=sessions[sessions$clan=="talek",]
sessions = sessions %>% arrange(date)

tblHyenas <- read.csv("Data/tblHyenas.csv", header=TRUE, stringsAsFactors=FALSE, fill=T)

fission <- read.table("Data/fission.csv", sep=",", header=TRUE, stringsAsFactors=FALSE, fill=T)
fission = fission[fission$Membership=="e",]

ny=27 # no. of years
ni=2400 # no. of individuals

first.year=2 # 2 for 1989
allhyenas=1
females=0
males=0

counts=array(0,c(ni,ni,ny))
years=array(0,c(ni,ny))
names=c("")
presence=array(0,c(ny,12000,ni))
ai=array(0,c(ni,ni,ny)) # association index
last.year=0

remove.fission=1 # remove hyenas that ended in Talek East after 2002
c=0 #count sessions

for (r in 1:length(sessions$hyenas)){
  row = unlist(strsplit(as.character(sessions$hyenas[r]),","))
  row=row[-1]
  date=as.Date(sessions$date[r],format="%Y-%m-%d")
  year=as.numeric(format(date,format="%Y"))-1987
  if (year>29) break
  if (year>last.year){
    print(paste(year,","))
    last.year=year
    obs=0
  }
  obs=obs+1
  
  if (year<first.year) next 
  if (year>first.year+(ny-1)) break
  
  c=c+1
  for (i in 1:length(row)){
    if (length(row)==0) next
    seen=match(row[i],names,nomatch=0)
    if (seen>0){
      next
    } 
    if (remove.fission==1){
      talek.east=match(row[i],fission$ID,nomatch=0)
      if (talek.east>0) next # removing Talek East hyenas
    }
    if (allhyenas==1) names=append(names,row[i]) # all only
    else if(females==1){ # females only
      isfemale = match(row[i],tblHyenas$id,nomatch=0)
      if (isfemale==0) {
        print(row[i])
        next
      }
      isfemale1 = tblHyenas$sex[isfemale]
      if (isfemale1 != 'f' && isfemale1 != 'F') next
      names=append(names,row[i])
    }
    else if (males==1){ # males only
      ismale = match(row[i],tblHyenas$id,nomatch=0)
      if (ismale==0) {
        print(row[i])
        next
      }
      ismale1 = tblHyenas$sex[ismale]
      if (ismale1 != 'm' && ismale1 != 'M') next
      names=append(names,row[i])
    }
  }
  l=length(row)-1
  if (l<0) next
  for (i in 1:l){
    a=match(row[i],names)
    years[a,year-(first.year-1)]=years[a,year-(first.year-1)]+1 #counts how many times a was observed in each year
    presence[year-(first.year-1),obs,a]=1
    if (l<1) break
    for (j in (i+1):(l+1)){
      b=match(row[j],names)
      years[b,year-(first.year-1)]=years[b,year-(first.year-1)]+1
      presence[year-(first.year-1),obs,b]=1
      counts[a,b,year-(first.year-1)]=counts[a,b,year-(first.year-1)]+1
      counts[b,a,year-(first.year-1)]=counts[b,a,year-(first.year-1)]+1
    }
  }
}

# Remove hyenas that were observed less than 20 times in a given year
ani=length(names)-1
ani
presence1=presence
presence1[presence1>1]=1
rm(presence)

presence1=presence1[,,2:(ani+1)]
counts1=counts[2:(ani+1),2:(ani+1),]

rm(counts)

j=1
seen.thre=20 #number of times a year an individual must be seen to be included
while (j <= ani){
  print(ani-j)
  ok=0
  for (i in 1:ny) {
    if (sum(presence1[i,1:12000,j])<seen.thre) {
      presence1[i,,j]=0
      counts1[j,,i]=0
      counts1[,j,i]=0
    }
    else ok=1
  }
  if (ok==0){ # was not observed at least 20 times in any year
    presence1=presence1[,,-j]
    counts1=counts1[-j,-j,]
    names=names[-(j+1)]
    ani=ani-1
    next
  }
  j=j+1
}
ani
names=names[-1]


# Calculate association index

for (i in 1:ny){
  for (a in 1:ani){
    for (b in 1:ani){
      ab.xor = sum(presence1[i,,a]>0 | presence1[i,,b]>0)
      if (ab.xor > 0) ai[a,b,i] = counts1[a,b,i]/ab.xor # simple ratio index
      else ai[a,b,i] = 0
      
    }
  }
}


# list of all yearly matrices 

all.wave=list()
for (i in 1:ny){
  all.wave[[i]]=ai[1:ani,1:ani,i]
  colnames(all.wave[[i]])=names
}

for (j in 1:ny){
  remove=c()
  for (i in 1:ani){
    if (sum(presence1[j,1:12000,i])<seen.thre){
      remove=c(remove,i) #noting who to remove
    }
  }
  all.wave[[j]]=all.wave[[j]][-remove,-remove]
}

save(presence1,file="presence1.Rdata")
rm(presence1)

for (i in 1:ny) cat(dim(all.wave[[i]])," ")
for (i in 1:ny) rownames(all.wave[[i]])=colnames(all.wave[[i]])

ai = ai[1:ani,1:ani,1:ny]
