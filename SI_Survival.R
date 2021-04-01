# Survival of offspring and mothers

library(survival)
library(broom)
library(viridis)
library(nlme)
library("piecewiseSEM")
library(qpcR)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(lmtest)


f1 = cent_by_pair %>% filter(sex=="f",al.year==1,last.age>365,!is.na(MRank))
f1 = f1 %>% ungroup()


cox.test = coxph(Surv(last.age) ~ MO.AI.cor * MRank, data = f1)
summary(cox.test)

cox.pred = expand.grid(MO.AI.cor = seq(-0.5,1,length.out = 20), MRank = seq.int(1,50))
cox.pred.lp = augment(cox.test, newdata = as_tibble(cox.pred))

ggplot(cox.pred.lp, aes(x=MRank, y= MO.AI.cor)) + geom_raster(aes(fill=.fitted),interpolate = T) + 
  scale_fill_gradientn(colors=viridis(10,alpha = 0.8)) + geom_point(data = f1,alpha=0.6) + 
  labs(x="Maternal rank",y="Mother-offspring AI correlation",fill = "fitted\nhazard\nratio") + theme_cowplot()

cox.pred2 = expand.grid(MO.AI.cor = seq(0,1,length.out = 20), MRank = c(1,30))
x <- survfit(cox.test, newdata = cox.pred2)
cox.pred2$Results <- summary(x)$table[,"median"]

cox.pred2 %>% filter(MO.AI.cor==1,MRank==1) %>% select(Results) - cox.pred2 %>% filter(MO.AI.cor==0,MRank==1) %>% 
  select(Results) # 3070 days difference in lifespan
cox.pred2 %>% filter(MO.AI.cor==1,MRank==30) %>% select(Results) - cox.pred2 %>% 
  filter(MO.AI.cor==0,MRank==30) %>% select(Results) # -867 days difference in lifespan

mod_si_ai = lmer(MO.AI.cor ~ log(MO.AI) + (1|Mother) + (1|Year), data = cent_by_pair[cent_by_pair$al.year=="1" & cent_by_pair$MO.AI>0,])

pred_mod_si_ai = cent_by_pair[cent_by_pair$al.year=="1" & cent_by_pair$MO.AI>0,]
pred_mod_si_ai$pred = as.numeric(predict(mod_si_ai,re.form = NA))
ggplot(pred_mod_si_ai, aes(x=MO.AI,y=MO.AI.cor,color=sex)) + geom_point(size=1.5, alpha = 0.7) +
  geom_line(aes(y=pred),color="black") + labs(x="Mother-offspring AI",y="Mother-offspring AI correlation") + 
  scale_color_brewer(palette = "Dark2") + theme_cowplot()

summary(mod_si_ai)

# Structural equation modeling
lav.dat2 = f1[complete.cases(f1),]

modlist1a= psem (
  lme (MO.AI.cor ~ MO.AI, random = ~1|Mother, data = lav.dat2, na.action = na.omit),
  lme (MO.AI ~ MRank, random = ~1|Mother, data = lav.dat2, na.action = na.omit),
  lme (last.age ~ MO.AI.cor + MO.AI + MRank, random = ~1|Mother,data = lav.dat2, na.action = na.omit )
)
summary(modlist1a)

modlist1c= psem (
  lme (MO.AI.cor ~ MO.AI, random = ~1|Mother, data = lav.dat2, na.action = na.omit),
  lme (MO.AI ~ MRank, random = ~1|Mother, data = lav.dat2, na.action = na.omit),
  glmer (last.age ~ MO.AI + MRank + (1|Mother), family=Gamma(link = "log"),data = lav.dat2, na.action = na.omit )
)
summary(modlist1c)

akaike.weights(c(30.09,30.67))

# coefficients from modlist1a
sem1 <- grViz("
digraph SEM {

graph [       overlap = true,
       outputorder = edgesfirst]

node [shape = rectangle]

a [pos = '0,1!', label = 'Mother-offspring\nassociation strength']
b [pos = '3,1!', label = 'Mother-offspring\nsocial similarity']
c [pos = '-3,1!', label = 'Maternal rank']
d [pos = '0,-1!', label = 'Offspring survival']

a->b [label = '0.67', penwidth = 2.01, arrowsize = 1.34]
a->d [label = '0.41', penwidth = 1.23, arrowsize = 0.82]
c->a [label = '-0.03', penwidth = 0.09, arrowsize = 0.15]
c->d [label = '-0.25', penwidth = 0.75, arrowsize = 0.50]

}
")

sem1 %>% export_svg %>% charToRaw %>% rsvg_pdf("sem1_2021.pdf")

# Social inheritance and maternal survival

# finding the last year of each hyena:
last.year.of = data.frame(hyena=character(), year=integer())
current.year=1989
for (i in 1:ny){
  net=all.wave[[i]]
  for (j in 1:dim(net)[1]){
    loc = match(colnames(net)[j],last.year.of$hyena, nomatch = 0)
    if (loc == 0){ # hyena not there yet
      new = data.frame(hyena=colnames(net)[j], year=current.year) 
      last.year.of = rbind(last.year.of,new)
    }
    else last.year.of$year[loc] = current.year
  }
  current.year = current.year + 1
}
last.year.of=mutate(last.year.of,before.last=year-1)

max.year=dplyr::slice(group_by(by_pair,pair),which.max(al.year))
before.max=dplyr::slice(arrange(group_by(by_pair,pair),al.year),which.max(al.year)-1)
before.max=mutate(before.max,which=1)
max.year=mutate(max.year,which=0)
max.year$year=max.year$Year+1988
before.max$year=before.max$Year+1988

max.year3 = semi_join(max.year,last.year.of, by = c("Mother"="hyena","year"))
before.max3 = semi_join(before.max,last.year.of, by = c("Mother"="hyena","year"="before.last"))
all.max3=bind_rows(max.year3,before.max3)

all.max3=mutate(all.max3, Survived=which)
all.max3 = all.max3 %>% mutate(Survived = ifelse(which==0,"No","Yes"))

ggplot(all.max3[all.max3$al.year<7,], aes(x=factor(al.year),y=MO.AI.cor,fill=Survived)) + geom_boxplot() + 
  scale_fill_brewer(name="Mother\nsurvived",palette="Set1") + ylab("Mother-offspring AI correlation") + xlab("Years after leaving the den")  + theme_cowplot()
theme(legend.justification=c(1,1),legend.position=c(0.96,0.96))

log_m=(glm(which~al.year+MO.AI.cor+MRank, family=binomial(link='logit'),data=all.max3[all.max3$al.year<7 & all.max3$sex=="f",]))
summary(log_m)
lrtest(log_m,2)

