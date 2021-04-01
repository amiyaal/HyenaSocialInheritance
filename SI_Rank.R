# Effects of maternal social rank

Fig.2A = ggplot(by_pair[by_pair$al.year==1,],aes(x=MRank,y=MO.AI.cor, group=MRank)) + geom_boxplot() +
  geom_smooth(aes(group=1),color = "blue", alpha = 0.5, method = "lm", se = T) + 
  xlab("Maternal rank") + ylab("Corrleation of mother and offspring AIs") + 
  theme(strip.background= element_rect(colour = "black", fill = "white")) + ggtitle("1st year")

Fig.2B = ggplot(by_pair[by_pair$al.year==2,],aes(x=MRank,y=MO.AI.cor, group=MRank)) + geom_boxplot() +
  geom_smooth(aes(group=1),color = "blue", alpha = 0.5, method = "lm", se = T) + 
  xlab("Maternal rank") + ylab("Corrleation of mother and offspring AIs") + 
  theme(strip.background= element_rect(colour = "black", fill = "white")) + ggtitle("2nd year")

Fig.2C = ggplot(by_pair[by_pair$al.year==3,],aes(x=MRank,y=MO.AI.cor, group=MRank)) + geom_boxplot() +
  geom_smooth(aes(group=1),color = "blue", alpha = 0.5, method = "lm", se = T) + 
  xlab("Maternal rank") + ylab("Corrleation of mother and offspring AIs") + 
  theme(strip.background= element_rect(colour = "black", fill = "white")) + ggtitle("3rd year")

Fig.2D = ggplot(by_pair[by_pair$al.year==4,],aes(x=MRank,y=MO.AI.cor, group=MRank)) + geom_boxplot() +
  geom_smooth(aes(group=1),color = "blue", alpha = 0.5, method = "lm", se = T) + 
  xlab("Maternal rank") + ylab("Corrleation of mother and offspring AIs") + 
  theme(strip.background= element_rect(colour = "black", fill = "white")) + ggtitle("4th year")

Fig.2E = ggplot(by_pair[by_pair$al.year==5,],aes(x=MRank,y=MO.AI.cor, group=MRank)) + geom_boxplot() +
  geom_smooth(aes(group=1),color = "blue", alpha = 0.5, method = "lm", se = T) + 
  xlab("Maternal rank") + ylab("Corrleation of mother and offspring AIs") + 
  theme(strip.background= element_rect(colour = "black", fill = "white")) + ggtitle("5th year")

Fig.2F = ggplot(y1.corr[!is.na(y1.corr$MRank),],aes(x=MRank, group=MRank,y=mean.diff))+geom_boxplot()+
  geom_smooth(method = "lm")+xlab("Maternal rank")+ylab("Mean difference of MO AIs") +
  geom_smooth(aes(group=1),color = "dark red", alpha = 0.5, method = "lm", se = T)

Fig.2 = plot_grid(Fig.2A, Fig.2B, Fig.2C,Fig.2D,Fig.2E,Fig.2F, labels=c("A","B","C","D","E","F"))
save_plot("Fig_2.pdf", Fig.2, base_width=4, base_height=4, ncol=3,nrow=2)

# Social inheritance and maternal rank in the 1st year
mod.rank2 = lmer(MO.AI.cor ~ MO.AI+ MRank + sex + (1|Mother) + (1|Year), data = by_pair[by_pair$al.year==1,])
summary(mod.rank2)

# Social inheritance and maternal rank in subsequent years
mod.rank3 = lmer(MO.AI.cor ~ MO.AI+ MRank + sex + (1|Mother/Offspring) + (1|Year), 
                 data = by_pair[by_pair$al.year!=1,])
summary(mod.rank3)

# Effect of maternal rank on offspring bonds
y1.corr=by_pair[by_pair$al.year==1,]

ggplot(y1.corr[!is.na(y1.corr$MRank),],aes(x=MRank, group=MRank,y=mean.diff))+geom_boxplot()+
  geom_smooth(method = "lm")+xlab("Rank of mother")+ylab("Mean difference of MO AIs") +
  geom_smooth(aes(group=1),color = "blue", alpha = 0.5, method = "lm", se = T)

diff.cor=rep(0,length(y1.corr$Offspring))
o.m.diff=rep(0,length(y1.corr$Offspring))

for (i in 1:length(y1.corr$Offspring)){
  y = y1.corr$Year[i]
  o = y1.corr$Offspring[i]
  m = y1.corr$Mother[i]
  net = all.wave[[y]]
  net.size=dim(net)[1]
  y.ranks = ranks[ranks$year==y+1988,] # get the ranks of this year
  ai.m = net[which(colnames(net)==m),]
  ai.o = net[which(colnames(net)==o),]
  ai.diff = ai.o - ai.m
  ai.diff = ai.diff[names(ai.diff)!=m]
  ai.diff = ai.diff[names(ai.diff)!=o]
  ai.m = ai.m[names(ai.m)!=o]
  ai.o = ai.o[names(ai.o)!=m]
  ai.m = ai.m[names(ai.m)!=m]
  ai.o = ai.o[names(ai.o)!=o]
  diff.df = data.frame(name=names(ai.diff),diff=ai.diff,ai.m=ai.m,ai.o=ai.o)
  dm = merge(diff.df,y.ranks,by.x="name",by.y="id",all.x=F)
  diff.cor[i]=cor(dm$rank,dm$diff)
  o.m.diff[i]=mean(ai.diff)
}
y1.corr$diff.cor=diff.cor
y1.corr$mean.diff=o.m.diff

diff.rank.mod = lmer(mean.diff ~ MRank + (1|Mother) + (1|Year), data = y1.corr)
summary(diff.rank.mod)
