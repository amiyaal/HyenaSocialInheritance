# Theoretical model

library(bbmle)

a.mx = alt.ai.df$X.Mother.AI
real.a.ox = alt.ai.df$Off.X.AI

si_proc_9 = function(pn, beta, beta.r, a.mx.i, a.zero.i, real.a.ox.i){
  real.a.ox.i = real.a.ox.i + 0.0001
  alpha.n.i = (beta*a.mx.i)/(1-a.mx.i) + 0.0001
  R1 = suppressWarnings(dbeta(real.a.ox.i,alpha.n.i,beta))
  alpha.r.i = (beta.r*a.zero.i)/(1-a.zero.i) + 0.0001
  R2 = suppressWarnings(dbeta(real.a.ox.i,alpha.r.i,beta.r))
  R = pn*R1 + (1-pn)*R2
  R
}

si_like_9 = function(pn, beta, beta.r, a.zero){
  res = rep(0,length(a.mx))
  for (i in 1:length(a.mx)){
    res[i] = si_proc_9(pn, beta, beta.r, a.mx[i], a.zero, real.a.ox[i])
  }
  -sum(log(res))
}

est = mle2(minuslogl = si_like_9, start = list(pn = 0.9, beta = 5, beta.r = 5, a.zero = 0.01),
           lower = c(0, 0.1,0.1,0.0001),upper = c(1, 200, 20, 0.5),method = "L-BFGS-B")
summary(est)

# Validating the model

beta.n = 98.06
beta.r = 10.6
a.zero.valid = 0.037
validation = data.frame(a.mx = rbeta(10000,beta.r,beta.n),a.ox = rep(0,10000))
for (i in 1:10000){
  if(sample(c(1,0),1,prob = c(pn,1-pn))) 
    validation$a.ox[i] = rbeta(1,(beta.n*validation$a.mx[i])/(1-validation$a.mx[i]),beta.n)
  else validation$a.ox[i] = rbeta(1,(beta.r*a.zero.valid)/(1-a.zero.valid),beta.r)
}                         
validation$a.zero = rep(a.zero.valid,10000)

si_proc_v = function(pn, beta, beta.r, a.mx.i, a.zero.i, real.a.ox.i){
  real.a.ox.i = real.a.ox.i + 0.0001
  a.zero.i = a.zero.i + 0.0001
  alpha.n.i = (beta*a.mx.i)/(1-a.mx.i) + 0.0001
  R1 = suppressWarnings(dbeta(real.a.ox.i,alpha.n.i,beta))
  alpha.r.i = (beta.r*a.zero.i)/(1-a.zero.i) + 0.0001
  R2 = suppressWarnings(dbeta(real.a.ox.i,alpha.r.i,beta.r))
  R = pn*R1 + (1-pn)*R2
  R
}

si_like_v = function(pn, beta, beta.r, a.zero){
  res = rep(0,10000)
  for (i in 1:10000){
    res[i] = si_proc_v(pn, beta, beta.r, validation$a.mx[i], validation$a.zero[i], validation$a.ox[i])
  }
  -sum(log(res))
}

pn.vals = seq(0.1,0.9, by=0.1)
valid.summary.df = data.frame(pn=rep(0,length(pn.vals)*5),est.pn=rep(0,length(pn.vals)*5))
a=1
for (pn.i in pn.vals){
  for (c in 1:5){
    validation = data.frame(a.mx = a.mx,a.ox = rep(0,length(a.mx)))
    validation$a.zero = a.zero.valid
    for (i in 1:length(a.mx)){
      if(sample(c(1,0),1,prob = c(pn.i,1-pn.i))) 
        validation$a.ox[i] = rbeta(1,(beta.n*validation$a.mx[i])/(1-validation$a.mx[i]),beta.n)
      else validation$a.ox[i] = rbeta(1,(beta.r*validation$a.zero[i])/(1-validation$a.zero[i]),beta.r)
    }                         
    est.valid = mle2(minuslogl = si_like_v, 
                     start = list(pn = pn.i, beta = beta.n, beta.r = beta.r, a.zero = a.zero.valid),
                     lower = c(0, 0.1,0.1,0.01),upper = c(1, 200,50,0.2),method = "L-BFGS-B")
    valid.summary.df$pn[a] = pn.i
    valid.summary.df$est.pn[a] = coef(est.valid)[1]
    a = a + 1
    cat(a)
  }
}

fig.mle.valid = ggplot(valid.summary.df, aes(pn, est.pn)) + geom_jitter(width=0.02) + theme_cowplot() + 
  scale_x_continuous(breaks = pn.vals,labels = pn.vals) + scale_y_continuous(breaks = pn.vals,labels = pn.vals) + 
  labs(x="pn",y="Inferred pn")
ggsave("MLE_valid2021.pdf",fig.mle.valid)
