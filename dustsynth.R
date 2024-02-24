library(bayestestR)
betas = read.csv("bpost_column.csv")
n = 40
amt = 50

mean(betas[5001:45000,1])
#beta = map_estimate(betas[,1], precision = 2^10)
beta = .139

tot.dust = rep(NA,2)
mean.inf = c()

for (l in 1:100)
{
days = c(-9:10)
#beta = sample(betas[5001:45000,1],1)

###4 individuals enter iso on days -9, -8, . . ., 10
t.iso = rep(0, n)#rep(c(-9:10),n/10)
m = length(t.iso)

###test positive one day before isolation
#t.pos = t.iso-1

###draw total length of infection and placement of peak using Kissler et al
time.topeak = rnorm(m,4.2,0.5)
inf.time = time.topeak + rnorm(m,7.3,0.6)

###draw iso times relative to infection times
U.shift = runif(m,rep(0,m),inf.time)
#iso.start = pos.loc + 1

###individuals leave isolation after 10 days
#iso.end = iso.start + 10

###fix infection times in absolute time
inf.start = - U.shift - 1
inf.peak = inf.start + time.topeak
inf.end = inf.start + inf.time

mean.inf = c(mean.inf, inf.start)
###draw theta for each individual
#theta = rexp(80, median(betas[,2]))
theta = rexp(m, beta)
dust = rep(NA, m)

###deacy rate
zeta = log(2)/7

###generate dust for each individual
for (i in 1:m)
  {
  left.survive = 0
  left.arrivals = 0
  right.survive = 0
  right.arrivals = 0
  
  ###generate dust particles from the pre-peak portion of shedding if any of it falls in [0,T] \cap [t.iso+t.iso+10]
  if (inf.peak[i] > max(0,t.iso[i])) {
    left.end = max(c(inf.start[i], 0, t.iso[i]))
    right.end = min(c(inf.peak[i],10,t.iso[i]+10))
    C.left = log(theta[i] + 1)/(inf.peak[i]-inf.start[i])
    left.lambda = (exp(-C.left*inf.start[i])/C.left*exp(C.left*right.end) - right.end)-(exp(-C.left*inf.start[i])/C.left*exp(C.left*left.end)- left.end)
    left.arrivals = rpois(1, amt*left.lambda)
    
    left.times = rep(NA, left.arrivals)
    k = 1
    
    ###If there are any arrivals, draw their arrival times with rejection sampling
    while (k <= left.arrivals)
      {
      proposal.x = runif(1,left.end, right.end)
      proposal.y = runif(1, 0, exp(log(theta[i]+1))-1)
      if (proposal.y < exp(log(theta[i]+1)*(proposal.x-inf.start[i])/(inf.peak[i]-inf.start[i]))-1)
        {
        left.times[k] = proposal.x
        k = k + 1
        }
      }
    }
  
    ###If there are any arrivals, compute how many survive using their arrival times and zeta  
    if (left.arrivals > 0)
      {
      left.flips = exp(-zeta*(10-left.times))
      left.survive = rbinom(1, left.arrivals, left.flips)
      }
  
  ###generate dust particles from post-peak shedding if any of it falls in [0,T] \cap [t.iso,t.iso+10]
  if (inf.peak[i] < min(10,t.iso[i]+10) & inf.end[i] > max(0,t.iso[i])) {
    left.end = max(c(0, inf.peak[i],t.iso[i]))
    right.end = min(c(inf.end[i],10,t.iso[i]+10))
    C.right = log(theta[i] + 1)/(inf.end[i]-inf.peak[i])
    right.lambda = (-exp(C.right*inf.end[i])/C.right*exp(-C.right*right.end) - right.end)-(-exp(C.right*inf.end[i])/C.right*exp(-C.right*left.end) - left.end)
    right.arrivals = rpois(1,amt*right.lambda)
    
    right.times = rep(NA, right.arrivals)
    h = 1
    
    ###If there are any arrivals, draw their arrival times with rejection sampling
    while(h <= right.arrivals)
      {
      prop.x = runif(1, left.end, right.end)
      prop.y = runif(1, 0, exp(log(theta[i]+1))-1)
      if (prop.y < exp(log(theta[i]+1)*(inf.end[i] - prop.x)/(inf.end[i] - inf.peak[i]))-1)
        {
        right.times[h] = prop.x
        h = h+1
        }
      }
    
    ###If there are any arrivals, compute how many survive using their arrival times and zeta  
    if (right.arrivals > 0)
      {
      right.flips = exp(-zeta*(10-right.times))
      right.survive = rbinom(1, right.arrivals, right.flips)
      }
  }
  
  ###Compute total amount of surviving dust
  dust[i] = left.survive + right.survive
}

tot.dust[l] = sum(dust)
}

hist(tot.dust, main = "Histogram of virus particles from 10 individuals with 10mp dust", xlab = "dust")
print(tot.dust)
print(mean(tot.dust))

dust.out = tot.dust

write.csv(dust.out, file="dust_40.csv")
