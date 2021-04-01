# Offspring-mother correlations

library(dplyr)
library(ggplot2)
library(cowplot)
library(ggridges)
library(lubridate)
library(lmerTest)
library(glmmTMB)
library(stringr)

# Finding correlation with mothers

ranks <- read.csv("Data/tblFemaleRanks.csv", header=TRUE, stringsAsFactors=FALSE, fill=T)

all.cor.years=vector(mode="numeric")
si.df = data.frame(Mother=character(),Offspring=character(),MO.AI=double(),MO.AI.cor=double(),
                   Year=integer(), Age=double(), MRank=double(), StanMRank=double(),stringsAsFactors = F)
for (i in 1:ny){
  net=all.wave[[i]]
  net.size=dim(net)[1]
  all.cor=vector(mode="numeric",length=net.size*(net.size-1)/2)
  c=0
  for (a in 1:net.size){
    for (b in (a+1):net.size){
      if (b>net.size) break
      if (sum(net[-c(a,b),a])==0 | sum(net[-c(a,b),b])==0) next()
      c=c+1
      all.cor[c]=cor(net[-c(a,b),a],net[-c(a,b),b])
      if(is.na(all.cor[c])){
        stop()
      }
    }
  }
  all.cor.years=c(all.cor.years,all.cor)
  
  for (j in 1:net.size){
    focal = rownames(net)[j]
    focal.loc = match(focal,tblHyenas$id,nomatch=0)
    if (focal.loc==0) {
      print(focal)
    }
    mother = tblHyenas$mom[focal.loc]
    if (!is.na(mother)){
      mother.loc = match(mother,rownames(net),nomatch = 0)
      if (mother.loc != 0){
        dob = ymd(tblHyenas$birthdate[focal.loc])
        now = ymd(paste0(i+1988,"-07-01"))
        tmp.age = as.duration(dob %--% now) / ddays(1)
        if (length(tmp.age)==0) tmp.age=NA 
        tmp.mrank=ranks$rank[ranks$year==(i+1988) & ranks$id==mother]
        tmp.stan.mrank=ranks$stan_rank[ranks$year==(i+1988) & ranks$id==mother]
        if (length(tmp.mrank)==0) tmp.mrank=NA 
        if (length(tmp.stan.mrank)==0) tmp.stan.mrank=NA 
        tmp.df=data.frame(Mother=mother, Offspring=focal, MO.AI=net[j,mother.loc],
                          MO.AI.cor=cor(net[-c(j,mother.loc),j],net[-c(j,mother.loc),mother.loc]),
                          Year=i, Age=tmp.age, MRank=tmp.mrank, StanMRank=tmp.stan.mrank,stringsAsFactors = F)
        si.df=rbind(si.df,tmp.df)
      }
    }
  }
}

# Toy example
MOC = tibble(m=runif(10),o=m+runif(10,-0.05,0.05),a=runif(10),b=runif(10))
ggplot(MOC,aes(m,o)) + geom_point(size=4) + geom_point(aes(x=a,y=b),color="red",size=4) +
  geom_smooth(aes(x=a,y=b),method = "lm",color="red",se=F) + geom_smooth(method = "lm",color="black",se=F) + 
  labs(x="Association strengths of individual 1",y="Association strengths of individual 2")
ggsave("AI cor toy example 2019.pdf")

# Social similarity vs. distance in maternal tree
sim.fun = function(x, y){
  cor(x, y)
}
all.cor.net.years=list()
mul.si.df = data.frame(A=character(),B=character(),AB.AI=double(),AB.AI.cor=double(),AB.AI.cos=double(),
                       Year=integer(), Tree.dist=integer(),stringsAsFactors = F)
for (i in 1:ny){
  net=all.wave[[i]]
  net.size=dim(net)[1]
  all.cor.net=net
  all.cor.net[,]=NA
  for (a in 1:net.size){
    for (b in (a+1):net.size){
      if (b>net.size) break
      if (sum(net[-c(a,b),a])==0 | sum(net[-c(a,b),b])==0) next()
      all.cor.net[a,b]=sim.fun(net[-c(a,b),a],net[-c(a,b),b])
      if(is.na(all.cor.net[a,b])){
        stop()
      }
      # checking if they exist in the tree and calculate distance
      if (class(try(V(talek)[colnames(net)[a]],silent=T)) == "try-error") {tmp.dist=Inf
      } else if (class(try(V(talek)[colnames(net)[b]],silent=T)) == "try-error") {tmp.dist=Inf
      } else 
      {tmp.dist = distances(talek,v=V(talek)[colnames(net)[a]],to=V(talek)[colnames(net)[b]],mode="all")[1,1]}
      tmp.mul.df=data.frame(A=colnames(net)[a],B=colnames(net)[b],AB.AI=net[a,b],AB.AI.cor=all.cor.net[a,b],
                            AB.AI.cos=lsa::cosine(net[-c(a,b),a],net[-c(a,b),b]),Year=i, Tree.dist=tmp.dist,
                            stringsAsFactors = F)
      mul.si.df=rbind(mul.si.df,tmp.mul.df)
    }
  }
  all.cor.net.years[[i]]=all.cor.net
}

# Social inheritance and maternal rank

rank.si.df=si.df[!is.na(si.df$MRank),]
rank.si.df$bins <- as.integer(cut(rank.si.df$MRank,breaks = 13,dig.lab = 0,include.lowest = T))
rank.si.df$nrank = as.numeric(levels(rank.si.df$MRank))[rank.si.df$MRank]

more.sex = tblHyenas %>% dplyr::select(id,sex)
new.si = merge(rank.si.df,more.sex,by.x = "Offspring", by.y = "id",all.x = T, all.y = F)
new.si = new.si[new.si$sex!="u",]
new.si = new.si[new.si$sex!="",]
new.si = new.si[!is.na(new.si$sex),]

by_pair=group_by(new.si,Mother,Offspring) %>% mutate(al.year=Year-min(Year)+1) %>% 
  mutate(pair=paste0(Mother,".",Offspring))

y1.corr=by_pair[by_pair$al.year==1,]

sort_concat <- function(a,b){
  temp = sort(c(a,b))
  temp = paste0(temp[1],temp[2])
  temp
}

mul.si.df <- mul.si.df %>% rowwise() %>% mutate(Pair = sort_concat(A,B))
mul.si.df$New.dist = mul.si.df$Tree.dist
mul.si.df$New.dist[mul.si.df$New.dist!=1]="Other"
mul.si.df$New.dist[mul.si.df$New.dist==1]="Mother-Offspring"
mul.si.df$New.dist[mul.si.df$New.dist=="Mother-Offspring"] = "MO"
mul.si.df$New.dist = as.factor(mul.si.df$New.dist)

si.mod=lmer(AB.AI.cor ~ New.dist + (1|Pair),data = mul.si.df)
summary(si.mod)

# Association between relatedness and social inheritance

relate.si.df = data.frame(Offspring=y1.corr$Offspring,Mother=y1.corr$Mother,
                          relate.ai.cor=rep(0,length(y1.corr$Offspring)),
                          AB.AI.cor=rep(0,length(y1.corr$Offspring)))
for (i in 1:length(y1.corr$Offspring))
{
  x = y1.corr$Offspring[i]
  net = all.wave[[(y1.corr$Year[i])]]
  x.ai = net[colnames(net)==x,]
  x.df = data.frame(name=names(x.ai),ai=x.ai)
  x.df = x.df %>% filter(name != x)
  x.dist = mul.si.df[mul.si.df$A==x | mul.si.df$B==x,] %>% select(A,B,AB.AI,AB.AI.cor,Tree.dist,Year,Pair) %>% 
    group_by(Pair) %>% filter(Year==min(Year))
  tmp.A = x.df %>% left_join(x.dist, by=c("name"="A")) %>% filter(!is.na(B))
  tmp.B = x.df %>% left_join(x.dist, by=c("name"="B")) %>% filter(!is.na(A))
  x.df = bind_rows(tmp.A %>% select(-B), tmp.B %>% select(-A))
  x.df = x.df %>% mutate(Tree.dist=ifelse(Tree.dist=="Inf",15,Tree.dist))
  relate.si.df$relate.ai.cor[i] = cor(1/x.df$Tree.dist,x.df$AB.AI)
  relate.si.df$AB.AI.cor[i] = y1.corr$MO.AI.cor[i]
}


relate.fig = relate.si.df %>% ggplot(aes(relate.ai.cor,AB.AI.cor)) + geom_point() + 
  labs(x="Correlation between genetic relatdness and AI",y="Mother-offspring association similarity") + 
  theme_cowplot() 
ggsave("Relatedness_effect_on_SI_2020.pdf", relate.fig)

relate.si.df %>% select(Mother) %>% distinct() # There are 82 different mothers
summary(lmer(AB.AI.cor ~ relate.ai.cor + (1|Mother), data = relate.si.df))

# Offspring bond dependence on maternal bonds

alt.ai.df = data.frame(Mother=character(),Offspring=character(),X=character(),X.Strength=double(),
                       Off.X.AI=double(),X.Mother.AI=double(),MO.AI=double(),Year=integer(), MRank=double())

for (i in 1:length(si.df$Offspring)) {
  year = si.df$Year[i]
  net = all.wave[[year]]
  Offspring=si.df$Offspring[i]
  
  for (j in 1:dim(net)[1]) {
    if (colnames(net)[j] != Offspring & colnames(net)[j] != si.df$Mother[i]) {
      X = colnames(net)[j]
      Off.loc = which(colnames(net)==Offspring)
      X.loc = which(colnames(net)==X)
      Mother.loc = which(colnames(net) == si.df$Mother[i])
      X.Strength = sum(net[X.loc,-c(Off.loc,Mother.loc)])
      tmp.df = data.frame(Mother=si.df$Mother[i],Offspring=Offspring,X=X,X.Strength=X.Strength,
                          Off.X.AI=net[Off.loc,X.loc],X.Mother.AI=net[X.loc,Mother.loc],Year=year, 
                          MRank=si.df$MRank[i])
      
      alt.ai.df = alt.ai.df %>% bind_rows(tmp.df)
    }
  }
  cat(i," ")
}


MO.AI.df = si.df %>% select(Mother, Offspring, Year, MO.AI)
alt.ai.df = alt.ai.df %>% select(-MO.AI) %>% inner_join(MO.AI.df, by = c("Mother","Offspring","Year")) %>% 
  distinct()

m1 <- glmmTMB(Off.X.AI ~ X.Mother.AI + X.Strength +(1|Mother/Offspring), data = alt.ai.df, ziformula=~1, 
              family=beta_family(link = "logit"))

summary(m1)

# Correlation between mother's AIs and offspring's AI in next year

si.next.df = data.frame(Mother=character(),Offspring=character(),MO.AI=double(),MO.AI.cor=double(),
                        Year=integer(), stringsAsFactors = F)
for (i in 1:(ny-1)){
  net=all.wave[[i]]
  net.size=dim(net)[1]
  next.net=all.wave[[i+1]]
  next.net.size=dim(next.net)[1]
  
  for (j in 1:next.net.size){
    focal = rownames(next.net)[j]
    focal.prev = match(focal,rownames(net),nomatch = 0)
    if (focal.prev != 0) next()  # we're interested only in new offspring
    focal.loc = match(focal,tblHyenas$id,nomatch=0)
    if (focal.loc==0) {
      print(focal)
    }
    mother = tblHyenas$mom[focal.loc]
    if (!is.na(mother)){
      mother.loc = match(mother,rownames(net),nomatch = 0)
      if (mother.loc != 0){
        ordered.net = net[order(rownames(net)),order(rownames(net))]
        ordered.next.net = next.net[order(rownames(next.net)),order(rownames(next.net))]
        # those that were present in both years:
        inter.names = intersect(rownames(ordered.net),rownames(ordered.next.net)) 
        focal.loc = match(focal,rownames(ordered.next.net),nomatch = 0)
        mother.loc = match(mother,rownames(ordered.net),nomatch = 0)
        mother.next.loc = match(mother,rownames(ordered.next.net),nomatch = 0)
        focal.ai = ordered.next.net[-c(focal.loc,mother.next.loc),focal.loc]
        focal.ai = focal.ai[names(focal.ai) %in% inter.names]
        focal.ai = focal.ai[!is.na(focal.ai)]
        mother.ai = ordered.net[-c(mother.loc),mother.loc]
        mother.ai = mother.ai[names(mother.ai) %in% inter.names]
        mother.ai = mother.ai[!is.na(mother.ai)]
        tmp.df=data.frame(Mother=mother, Offspring=focal, 
                          MO.AI=ifelse(mother.next.loc==0,NA,ordered.next.net[focal.loc,mother.next.loc]),
                          MO.AI.cor=cor(focal.ai,mother.ai),
                          Year=i, stringsAsFactors = F)
        si.next.df=rbind(si.next.df,tmp.df)
      }
    }
  }
}
si.next.df=si.next.df[!is.na(si.next.df$MO.AI.cor),]

ggplot(si.next.df, aes(x=MO.AI.cor)) + geom_density(alpha=.9,fill="dark orange") + 
  xlab("Correlation in associaion index with others") + ylab("Density")  + 
  geom_vline(xintercept = 0, linetype="dashed")
ggsave("next_year_mother_offspring_AI_correlation 2019.pdf")

# Association with agemates

agemates.df = data.frame(A=character(),B=character(),AI=double(), AI.cor=double(), Year=integer(), 
                         AgeA=double(), AgeB=double())

for (i in 1:ny){
  net=all.wave[[i]]
  net.size=dim(net)[1]
  
  for (j in 1:(net.size-1)){
    A = rownames(net)[j]
    A.loc = match(A,tblHyenas$id,nomatch=0)
    if (A.loc==0) {print(A)}
    else {
      dob = ymd(tblHyenas$birthdate[A.loc])
      now = ymd(paste0(i+1988,"-07-01"))
      tmp.age.A = as.duration(dob %--% now) / ddays(1)
      if (length(tmp.age.A)==0) tmp.age.A=NA 
    }
    
    for (k in (j+1):net.size) {
      B = rownames(net)[k]
      B.loc = match(B,tblHyenas$id,nomatch=0)
      if (B.loc==0) {print(B)}
      else {
        dob = ymd(tblHyenas$birthdate[B.loc])
        now = ymd(paste0(i+1988,"-07-01"))
        tmp.age.B = as.duration(dob %--% now) / ddays(1)
        if (length(tmp.age.B)==0) tmp.age.B=NA 
      }
      if (sum(net[-c(j,k),j])==0 | sum(net[-c(j,k),k])==0) next()
      cor.temp = sim.fun(net[-c(j,k),j],net[-c(j,k),k])
      tmp.df = data.frame(A=A,B=B,AI=net[j,k], AI.cor=cor.temp, Year=i, AgeA=tmp.age.A, AgeB=tmp.age.B)
      agemates.df = rbind(agemates.df, tmp.df)
    }
  }
  cat(i,append = T)
}

agemates.df$AgeDiff = abs(agemates.df$AgeA - agemates.df$AgeB)
agemates.df = agemates.df %>% mutate(YearDiff = floor(AgeDiff/365))

agemates.df = agemates.df %>% mutate(SexA = tblHyenas$sex[match(agemates.df$A,tblHyenas$id)], 
                                     SexB = tblHyenas$sex[match(agemates.df$B,tblHyenas$id)],
                                     PairSex = paste0(SexA,SexB)) %>% 
  mutate(PairSex = ifelse(PairSex=="mf","fm",PairSex))

agemates.df = agemates.df %>% mutate(PairSex = toupper(PairSex))

agemates.df = agemates.df %>% mutate(Age.mates = ifelse(YearDiff==0,"Yes","No"))

age.ai.fig = agemates.df %>% filter(PairSex %in% c("FF","MM","FM"),!is.na(Age.mates)) %>% 
  mutate(Age.mates = factor(Age.mates)) %>% mutate(Age.mates = relevel(Age.mates,"Yes")) %>%  
  ggplot(aes(PairSex,AI.cor,fill=Age.mates)) + geom_boxplot() + theme_cowplot() + xlab("Sex combination") + 
  ylab("Correlation in AIs") + 
  scale_fill_manual("",labels = c("Agemates","Non agemates"), values = c("blue", "red")) + 
  theme(legend.position="top")
ggsave("AI.Age2020.pdf",age.ai.fig)

ex.aicor.fig = agemates.df %>% filter(PairSex %in% c("FF","MM","FM"),!is.na(Age.mates)) %>% 
  mutate(Age.mates = factor(Age.mates)) %>% mutate(Age.mates = relevel(Age.mates,"Yes")) %>%  
  ggplot(aes(AI.cor,fill=Age.mates)) + geom_density(alpha=0.7) + theme_cowplot() + ylab("Density") + 
  xlab("Correlation in AIs") + 
  scale_fill_manual("",labels = c("Agemates","Non agemates"), values = c("blue", "red")) + 
  theme(legend.position="top") + facet_wrap(~PairSex)
ggsave("Sex.CorAI.Age2020.pdf",sex.aicor.fig)

# Centrality of mothers and offspring

cent.df = data.frame(Mother=character(),Offspring=character(),MO.AI=double(),MO.AI.cor=double(),
                     Year=integer(), Age=double(), MRank=double(), Stan_MRank=double(), M.stren=double(),
                     M.eigen=double(),M.cc=double(),O.stren=double(),O.eigen=double(),O.cc=double(),
                     stringsAsFactors = F)

for (i in 1:ny){
  net=all.wave[[i]]
  net.size=dim(net)[1]
  inet=graph_from_adjacency_matrix(net,mode="undirected",weighted = T)
  
  for (j in 1:net.size){
    focal = rownames(net)[j]
    focal.loc = match(focal,tblHyenas$id,nomatch=0)
    if (focal.loc==0) {
      print(focal)
    }
    mother = tblHyenas$mom[focal.loc]
    if (!is.na(mother)){
      mother.loc = match(mother,rownames(net),nomatch = 0)
      if (mother.loc != 0){
        dob = ymd(tblHyenas$birthdate[focal.loc])
        now = ymd(paste0(i+1988,"-07-01"))
        tmp.age = lubridate::as.duration(dob %--% now) / ddays(1)
      
        if (length(tmp.age)==0) tmp.age=NA 
        tmp.mrank=ranks$rank[ranks$year==(i+1988) & ranks$id==mother]
        tmp.stan_mrank=ranks$stan_rank[ranks$year==(i+1988) & ranks$id==mother]
      
        if (length(tmp.mrank)==0) tmp.mrank=NA 
        if (length(tmp.stan_mrank)==0) tmp.stan_mrank=NA 
        tmpc.df=data.frame(Mother=mother, 
                           Offspring=focal, 
                           MO.AI=net[j,mother.loc],
                           MO.AI.cor=cor(net[-c(j,mother.loc),j],net[-c(j,mother.loc),mother.loc]),
                           Year=i, Age=tmp.age, 
                           MRank=tmp.mrank, 
                           Stan_MRank=tmp.stan_mrank, 
                           M.stren=strength(inet)[mother.loc]/net.size,
                           M.eigen=eigen_centrality(inet)$vector[mother.loc],
                           M.cc=transitivity(inet,"weighted",isolates = "zero")[mother.loc],
                           O.stren=strength(inet)[j]/net.size,O.eigen=eigen_centrality(inet)$vector[j],
                           O.cc=transitivity(inet,"weighted",isolates = "zero")[j],stringsAsFactors = F)
        cent.df=rbind(cent.df,tmpc.df)
      }
    }
  }
}

ggplot(cent.df,aes(x=M.stren,y=O.stren))+geom_point(alpha=0.4) + geom_smooth(method = "lm") + 
  labs(x="Maternal strength centrality", y="Offspring strength centrality") 
ggsave("mo_centrality_2019.pdf")

# ontogeny of centrality and social inheritance

# Calculate the last age of each individual
last.seen = data.frame(id=names,last.seen=NA)
for (i in 1:length(names)){
  id.sessions = sessions[str_detect(sessions$hyenas,names[i]),]
  id.sessions = id.sessions[!is.na(id.sessions$hyenas),]
  last.seen$last.seen[i] = tail(id.sessions,1)$date
}

last.seen$last.age = NA
last.seen$sex = NA
for (i in 1:length(names)){
  focal.loc = match(last.seen$id[i],tblHyenas$id,nomatch=0)
  dob = ymd(tblHyenas$birthdate[focal.loc])
  last = ymd(last.seen$last.seen[i])
  last.seen$last.age[i] = as.duration(dob %--% last) / ddays(1)
  last.seen$sex[i] = tblHyenas$sex[focal.loc]
}

# removing those who are still alive:
for (i in 1:length(names)){
  if (year(last.seen$last.seen[i]) >= 2015) last.seen$last.age[i] = NA  
}

cent.last.age=merge(cent.df,last.seen,by.x = "Offspring",by.y = "id",all.y=F)
cent.last.age=filter(cent.last.age,sex=="f" | sex=="m")

cent.last.age=mutate(cent.last.age,pair=paste0(Offspring,".",Mother))

cent_by_pair=group_by(cent.last.age,pair)
cent_by_pair=mutate(cent_by_pair,al.year=Year-min(Year)+1)

nlabels = c(table(as.factor(mul.si.df$Tree.dist)))
plabels = paste0(names(nlabels),"  n=",nlabels)
plabels=plabels[-c(8,9,10)]

ggplot(mul.si.df[mul.si.df$Tree.dist!="8" & mul.si.df$Tree.dist!="9" & mul.si.df$Tree.dist!="10",], 
       aes(x=AB.AI.cor,y=as.factor(Tree.dist), fill=as.factor(Tree.dist))) + 
  geom_density_ridges(rel_min_height = 0.01) + scale_x_continuous(limits=c(NA,1),expand = c(0.01, 0)) + 
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_fill_brewer(direction=-1,labels=plabels) +
  labs(x="Correlation in\n association index",y="Distance in\n maternal tree",fill="Sample size") +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14),legend.text=element_text(size=14),
        legend.position="none") + theme_cowplot()

ggplot(by_pair %>% filter(sex!="l"),aes(x=factor(al.year),y=MO.AI.cor,fill=sex))+geom_boxplot()+
  xlab("Years after offspring left den")+ylab("MO correlation in AI") + theme_cowplot() +
  theme(legend.position = c(0.8, 0.9)) + scale_fill_brewer(palette="Set1",direction = 1) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14),legend.text=element_text(size=14))

ggplot(by_pair[by_pair$sex!="l",],aes(x=factor(al.year),y=MO.AI,fill=sex))+geom_boxplot()+
  xlab("Years after offspring left den")+ylab("MO AI") + scale_fill_brewer(palette="Set1",direction = 1) + 
  theme_cowplot() + theme(legend.position = c(0.8, 0.9)) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14),legend.text=element_text(size=14))

moc.ont = lmer(MO.AI.cor ~ al.year*sex + (1|pair), data = by_pair[by_pair$al.year<7,])
summary(moc.ont)

mo.ai.ont = lmer(MO.AI ~ al.year*sex + (1|pair), data = by_pair[by_pair$al.year<7,])
summary(mo.ai.ont)

median(by_pair[by_pair$sex!="l" & by_pair$al.year==1,]$MO.AI)
median(by_pair[by_pair$sex!="l" & by_pair$al.year==6,]$MO.AI)
