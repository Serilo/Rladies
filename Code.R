# Libro de referencia: Bayesian Computation with R, Albert Jim -----------------------

library(LearnBayes)
library(ggplot2)
library(knitr)
library(tidyr)
library(dplyr)
library(LaplacesDemon)



# Capítulo 2: Introducción al pensamiento Bayesiano -----------------------

Param_previa<-  function(q,p,init, dist){
  previa_optim <- function(param) {
    (dist(p[1],param[1],param[2])-q[1])^2 + (dist(p[2],param[1],param[2])-q[2])^2
  }
  r <-  optim(init,previa_optim)
  r <-  unlist(r)
  r <-  c(r[1],r[2])
  return(data.frame(r))
}

res <- Param_previa(c(0.3,0.5), c(0.5,0.9),c(1,1), qbeta)

p <-  seq(0, 1, length = 500)
a <-  res$r[1]
b <-  res$r[2]
s = 11
f = 16

df <- data.frame(p = seq(0, 1, length = 500),Previa= dbeta(p,a,b),Verosimilitud= dbeta(p,s+1,f+1),Posterior= dbeta(p,a+s,b+f))
df <- gather(df,key = "Dist", value = "Densidad", -p)

ggplot(data=df, aes(x=p, y=Densidad, colour=Dist))+
  geom_line(lwd=1)+
  theme_bw()+
  ggtitle("Comparación")

# Distribución predictiva

ab=c(a, b)
m=20; xs=0:20
pred=pbetap(ab, m, xs)

p=rbeta(10000,a, b)
x = rbinom(10000, 20, p)
table(x)

freq=table(x)
xs=c(0:max(x))
predprob=freq/sum(freq)


df <- data.frame(x=0:(length(predprob)-1),Pred_sim=as.vector(predprob), Pred_teo=pred[0:length(predprob)])

df <- gather(df,"Densidad", "Probabilidad_predictiva",-x)

ggplot(data=df, aes(x=x, y=Probabilidad_predictiva, color=Densidad))+
  geom_point()+
  labs(title=("Probabilidad predictiva del número de dormilones"))+
  geom_vline(xintercept = 5, colour="red")+
  facet_wrap(~ Densidad,dir = "v")+
  theme_bw()+
  theme(legend.position="none")


# Capítulo 3: Modelos Univariados -----------------------------------------

## Distribución Normal con media conocida y varianza desconocida

data(footballscores)
attach(footballscores)
d <- favorite - underdog - spread
n <- length(d)
v <- sum(d^2)

P <- rchisq(1000, n)/v
s <- sqrt(1/P)

quantile(s, probs = c(0.025, 0.5, 0.975))

ggplot(data.frame(s), aes(x=s))+ 
geom_histogram(aes(y=..density..),color="darkblue", fill="white")+
geom_density(alpha=.2, fill="#FF6666")+ 
theme_bw()+
labs(x="s", y="Frequencia", title="Histograma de las simulaciones de la desviación estándar")




## Estimación de la tasa de mortalidad en trasplantes de corazón

alpha<-16;beta<-15174
yobs<-1; ex<-66
y<-0:10
lam<-alpha/beta
py<-dpois(y, lam*ex)*dgamma(lam, shape = alpha,rate = beta)/dgamma(lam, shape= alpha + y,rate = beta + ex)
cbind(y, round(py, 3))

lambdaA <- rgamma(1000, shape = alpha + yobs, rate = beta + ex)

ex <- 1767; yobs<-4
y <- 0:10
py <- dpois(y, lam * ex) * dgamma(lam, shape = alpha,rate = beta)/dgamma(lam, shape = alpha + y, rate = beta + ex)
cbind(y, round(py, 3))

lambdaB <- rgamma(1000, shape = alpha + yobs, rate = beta + ex)
lambda <- data.frame(lambdaA, lambdaB)
lambda <- gather(lambda) 
 
grid <- with(lambda, seq(min(value), max(value), length = 1000))
 
gammad <- lambda %>%group_by(key) %>% mutate(
     value = grid,
     density = dgamma(grid, alpha, beta))


ggplot(lambda, aes(x=value))  + 
  geom_histogram(aes(y=..density..), color="darkblue", fill="lightblue") + 
  geom_line(aes(y = density),data=gammad, colour = "red") +
  facet_wrap(~ key,dir = "v")+
  theme_bw()


## Una ilustración de la robustez bayesiana

res <- Param_previa(c(80,120),c(0.05,0.95),c(1,1),qnorm)
res

mu <- res$r[1]; tau <- res$r[2] 
sigma <- 15 
n<-4
se <- sigma/sqrt(4)
ybar <- c(110, 125, 140)
tau1 <- 1/sqrt(1/se^2 + 1/tau^2)
mu1 <- (ybar/se^2 + mu/tau^2)*tau1^2
summ1 <- cbind(ybar, mu1, tau1)
summ1

Param_previa_ts <-  function(q,p,init) {
  previa_optim <- function(param) {
    (qst(p=p[1],mu=param[1],sigma=param[2],nu=2)-q[1])^2 + (qst(p=p[2],mu=param[1],sigma=param[2],nu=2)-q[2])^2
  }
  r <-  optim(init,previa_optim)
  r <-  unlist(r)
  r <-  c(r[1],r[2])
  return(data.frame(r))
}

rest <- Param_previa_ts(c(80,120),c(0.05,0.95),c(1,1))
rest

theta <-  seq(60, 140, length = 200)

mu <- rest$r[1]
tscale <-  rest$r[2]

df <- data.frame(Theta=seq(60, 140, length = 200),
                 tstudent=1/tscale*dt((theta-mu)/tscale,2),
                 normal=1/10*dnorm((theta-mu)/tau))
df <- gather(df,key = "Dist",
             value = "Densidad",
             -Theta)

ggplot(data=df, aes(x=Theta, y=Densidad, colour=Dist))+
  geom_line()+
  theme_bw(base_size=18)+
  ggtitle("Comparación de previas")+
  theme(legend.position = "bottom")

summ2 = c()
postts <- matrix(nrow=500,ncol=3)
for(i in 1:3){
  theta = seq(60, 180, length = 500)
  like = dnorm((theta - ybar[i])/7.5)
  prior = dt((theta - mu)/tscale, 2)
  post = prior * like
  post = post/sum(post)
  postts[,i] <-post 
  m = sum(theta * post)
  s = sqrt(sum(theta^2 * post) - m^2)
  summ2 = rbind(summ2, c(ybar[i], m, s))
}

postts <- cbind(c(rep(110,500),rep(125,500),rep(140,500)),c(postts[,1:3]))

normpost = matrix(nrow=500,ncol=3)
normpost = sapply(1:3, function(i) normpost[,i] = dnorm(theta, mu1[i], tau1)/sum(dnorm(theta, mu1[i], tau1)))
normpost <- cbind(c(rep(110,500),rep(125,500),rep(140,500)),c(normpost[,1:3]))

moments <- data.frame(round(summ1,2),round(summ2[,-1],2))
names(moments) <- c("y","E_n", "Sd_n", "E_ts","Sd_ts")
moments

df <- data.frame(Theta=rep(seq(60, 180, length = 500),6), Tstudent_post=postts[,2],Normal_post=normpost[,2], y=c(rep(110,500),rep(125,500),rep(140,500)))
df <- gather(df, "Dist", "Densidad_posterior", -Theta,-y)

ggplot(data=df, aes(x=Theta, y=Densidad_posterior, colour=Dist))+
  geom_line()+
  theme_bw(base_size=14)+
  ggtitle("Comparación de posteriores")+
  theme(legend.position = "bottom")+
  facet_wrap(~ y,dir = "h")


