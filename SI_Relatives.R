# Relatives

new.si$tree.size3=NA
family.size = 3
for (i in 1:length(new.si$Offspring))
{
  offspring = new.si$Offspring[i]
  mom = new.si$Mother[i]
  year = new.si$Year[i]
  comp = components(talek)$membership[names(components(talek)$membership) == offspring] #finding the relevant tree
  net = all.wave[[year]]
  net.size = dim(net)[1]
  curr.tree = induced_subgraph(talek,subcomponent(talek,offspring))
  rele.tree = induced_subgraph(curr.tree, which(V(curr.tree)$name %in% row.names(net)))
  dists = distances(as.undirected(rele.tree))
  off.dist = sum(dists[row.names(dists)==offspring] < family.size) - 1
  new.si$tree.size3[i] = off.dist
}

si.mod.3 = lmer(MO.AI.cor ~ scale(MO.AI) + scale(MRank) + sex + scale(tree.size3) + (1|Mother) + (1|Year), 
                data = new.si ,REML = T)
summary(si.mod.3)

new.si$tree.size5=NA
family.size = 5
for (i in 1:length(new.si$Offspring))
{
  offspring = new.si$Offspring[i]
  mom = new.si$Mother[i]
  year = new.si$Year[i]
  comp = components(talek)$membership[names(components(talek)$membership) == offspring] #finding the relevant tree
  net = all.wave[[year]]
  net.size = dim(net)[1]
  curr.tree = induced_subgraph(talek,subcomponent(talek,offspring))
  rele.tree = induced_subgraph(curr.tree, which(V(curr.tree)$name %in% row.names(net)))
  dists = distances(as.undirected(rele.tree))
  off.dist = sum(dists[row.names(dists)==offspring] < family.size) - 1
  new.si$tree.size5[i] = off.dist
}

si.mod.5 = lmer(MO.AI.cor ~ scale(MO.AI) + scale(MRank) + sex + scale(tree.size5) + (1|Mother) + (1|Year), data = new.si ,REML = T)
summary(si.mod.5)

relate.si.df = data.frame(Offspring=y1.corr$Offspring,Mother=y1.corr$Mother,relate.ai.cor=rep(0,length(y1.corr$Offspring)),AB.AI.cor=rep(0,length(y1.corr$Offspring)))
for (i in 1:length(y1.corr$Offspring))
{
  x = y1.corr$Offspring[i]
  net = all.wave[[(y1.corr$Year[i])]]
  x.ai = net[colnames(net)==x,]
  x.df = data.frame(name=names(x.ai),ai=x.ai)
  x.df = x.df %>% filter(name != x)
  x.dist = mul.si.df[mul.si.df$A==x | mul.si.df$B==x,] %>% select(A,B,AB.AI,AB.AI.cor,Tree.dist,Year,Pair) %>% group_by(Pair) %>% filter(Year==min(Year))
  tmp.A = x.df %>% left_join(x.dist, by=c("name"="A")) %>% filter(!is.na(B))
  tmp.B = x.df %>% left_join(x.dist, by=c("name"="B")) %>% filter(!is.na(A))
  x.df = bind_rows(tmp.A %>% select(-B), tmp.B %>% select(-A))
  #x.df = x.df %>% filter(Tree.dist != "Inf")
  x.df = x.df %>% mutate(Tree.dist=ifelse(Tree.dist=="Inf",15,Tree.dist))
  relate.si.df$relate.ai.cor[i] = cor(1/x.df$Tree.dist,x.df$AB.AI)
  relate.si.df$AB.AI.cor[i] = y1.corr$MO.AI.cor[i]
}
relate.si.df %>% ggplot(aes(relate.ai.cor,AB.AI.cor)) + geom_point() + labs(x="Correlation between genetic relatdness and AI",y="Mother-offspring association similarity") + theme_cowplot() 
relate.si.df %>% select(Mother) %>% distinct() # There are 82 different mothers
summary(lmer(AB.AI.cor ~ relate.ai.cor + (1|Mother), data = relate.si.df))
