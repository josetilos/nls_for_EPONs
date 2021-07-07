
# Chaotic map model
# Jose Alberto Hernandez
# July 2021



# Function

getChaoticMap <- function (m1=2,m2=2,d,Nsim) {
  x0 = runif(1); x = x0*(1:Nsim); y = as.numeric(x0<d)*(1:Nsim)
  for (kk in 2:Nsim) {
    x[kk] = ifelse(x[kk-1]<d,x[kk-1]+(1-d)*(x[kk-1]/d)^m1,x[kk-1]-d*((1-x[kk-1])/(1-d))^m2)
    y[kk] = ifelse(x[kk]<=d,0,1)
  }
  return(y);
}



# Model settings


m1 = 2; m2 = 2



UserRate = 100e6 # 100 Mb/s
deltat = 0.5e-6 # each time-slot is 0.5 us
Cuspon = 1.25e-9

# aa gives ON/OFF periods, aadur gives the duration of each ON periods

load = 0.1
aa = getChaoticMap(m1,m2,1-load,1e4) #  ON/OFF periods, multiply by *(UserRate*deltat) to get bits 





# Generic parameters

Nonu = 16
Guard = 1e-6
Load = 0.6
Load_peronu = Load/Nonu

Upstream_capacity = 1e9
PacketSize = 1250 # 1250 Bytes per packet
time_slot = PacketSize/1e9 # 1 time-slot is the size of a packet, i.e. 10us at 1G-EPON (or 1us at 10G-EPON)

# Poisson simulator (arrivals not correlated, just random)

Nslots = 1e5
Poisson_packets = 0*(1:Nslots)
Poisson_packets[sample((1:Nslots),round(Load_peronu*Nslots),replace=F)]=1

acf(Poisson_packets)


# Chaotic map model:

Chaotic_packets = getChaoticMap(2,2,1-Load_peronu,1e5)
acf(Chaotic_packets)






