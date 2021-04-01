library(igraph)
library(dplyr)
library(igraph)
library(cowplot)

# Number of observations per individual

sum.presence = tibble(Year = integer(),Sum=integer())
for (y in 1:ny){
  for (i in 1:482){
    sum.presence = add_row(sum.presence,Year=y,Sum=sum(presence1[y,,i]))
  }
}
sum.medians = sum.presence %>% filter(Sum>=20) %>% group_by(Year) %>% summarise(med = median(Sum))
plot.obs = sum.presence %>% filter(Sum>=20) %>% ggplot(aes(x=as.factor(Year),y=Sum)) + geom_boxplot() + 
  geom_hline(yintercept = 20, linetype="dashed", color="dark grey") + theme_cowplot() + 
  labs(x="Year",y="Number of times observed") + 
  geom_text(data = sum.medians, aes(x=Year,y=550, label=round(med)),size=3)

ggsave("n.obs.2021.pdf",plot.obs)

# Individuals observed less than 20 times

sum.presence0 = tibble(Year = integer(),Sum=integer())
for (y in 1:ny){
  for (i in 1:ani){
    sum.presence0 = add_row(sum.presence0,Year=y,Sum=sum(presence[y,,i]))
  }
}

sum.presence0 %>% filter(Sum<20,Sum>1) %>% 
  ggplot(aes(x=as.factor(Year),y=Sum)) + geom_boxplot() + ylab("Number of times observed if < 20") + 
  xlab("Year") + theme_cowplot()
ggsave("omitted_2021.pdf")

# Maternal pedigree

talek = make_empty_graph()
for (i in 1:length(tblHyenas$id))
{
  mom = tblHyenas$mom[i]
  if (is.na(mom) | (grepl("\\?",mom) == 1)) next()
  mom.node = mom %in% V(talek)$name
  if (!mom.node) {
    talek = talek + vertex(mom, sex="f")
  }
  offspring = tblHyenas$id[i]
  if (is.na(offspring) | (grepl("\\?",offspring) == 1)) next()
  offspring.node =  offspring %in% V(talek)$name
  if (!offspring.node) talek = talek + vertex(offspring, sex=tblHyenas$sex[i])
  talek = talek + edge(V(talek)[V(talek)$name == mom],V(talek)[V(talek)$name == offspring])
}
talek = simplify(talek)
colrs <- c("red","blue","white")
V(talek)$sex[V(talek)$sex == "f"] = as.numeric(1)
V(talek)$sex[V(talek)$sex == "m"] = as.numeric(2)
V(talek)$sex[is.na(V(talek)$sex)] = as.numeric(3)
V(talek)$sex[V(talek)$sex == "u"] = as.numeric(3)
V(talek)$sex = as.numeric(V(talek)$sex)
V(talek)$color = colrs[as.numeric(V(talek)$sex)]
talek.for.plot=delete_vertices(talek, V(talek)[V(talek)$color == "white"])
pdf(file = "Maternal pedigree 2019.pdf", width = 7,height = 7)
plot(talek.for.plot,vertex.size=3, vertex.label=NA,edge.arrow.size=.1,asp=0,
     layout=layout_as_tree(talek.for.plot, circular=TRUE))
dev.off()


# Association with age mates

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

agemates.df = agemates.df %>% mutate(YearDiff = floor(AgeDiff/365))

cor.age.fig = agemates.df %>% filter(!is.na(YearDiff)) %>% ggplot(aes(as.factor(YearDiff),AI.cor)) + 
  geom_boxplot() + theme_cowplot() + xlab("Age difference (years)") + ylab("Correlation in AIs")

ggsave("CorAI.Age2020.pdf",cor.age.fig)

agemates.df = agemates.df %>% mutate(SexA = tblHyenas$sex[match(agemates.df$A,tblHyenas$id)], 
                                     SexB = tblHyenas$sex[match(agemates.df$B,tblHyenas$id)],
                                     PairSex = paste0(SexA,SexB)) %>% 
  mutate(PairSex = ifelse(PairSex=="mf","fm",PairSex))

agemates.df = agemates.df %>% mutate(PairSex = toupper(PairSex))

agemates.df %>% filter(PairSex %in% c("FF","MM","FM"),!is.na(Age.mates)) %>%
  mutate(Age.mates = factor(Age.mates)) %>% mutate(Age.mates = relevel(Age.mates,"Yes")) %>%
  ggplot(aes(PairSex,AI.cor,fill=Age.mates)) + geom_boxplot() + theme_cowplot() + xlab("Sex combination") + 
  ylab("Correlation in AIs") + 
  scale_fill_manual("",labels = c("Agemates","Non agemates"), values = c("blue", "red")) + 
  theme(legend.position="top")

age.ai.fig = agemates.df %>% filter(PairSex %in% c("FF","MM","FM"),!is.na(Age.mates)) %>% 
  mutate(Age.mates = factor(Age.mates)) %>% mutate(Age.mates = relevel(Age.mates,"Yes")) %>% 
  ggplot(aes(PairSex,AI,fill=Age.mates)) + geom_boxplot() + theme_cowplot() + xlab("Sex combination") + 
  ylab("AI") + scale_fill_manual("",labels = c("Agemates","Non agemates"), values = c("blue", "red")) + 
  theme(legend.position="top") + scale_y_log10()

ggsave("AI.Age2020.pdf",age.ai.fig)

sex.aicor.fig = agemates.df %>% filter(PairSex %in% c("FF","MM","FM"),!is.na(Age.mates)) %>% 
  mutate(Age.mates = factor(Age.mates)) %>% mutate(Age.mates = relevel(Age.mates,"Yes")) %>% 
  ggplot(aes(AI.cor,fill=Age.mates)) + geom_density(alpha=0.7) + theme_cowplot() + ylab("Density") + 
  xlab("Correlation in AIs") + 
  scale_fill_manual("",labels = c("Agemates","Non agemates"), values = c("blue", "red")) + 
  theme(legend.position="top") + facet_wrap(~PairSex)

ggsave("Sex.CorAI.Age2020.pdf",sex.aicor.fig)




