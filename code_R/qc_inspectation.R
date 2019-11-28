qcs <- read.csv("QCs.csv")
# table(qcs$days)
outs == 0
for(i in 2:ncol(qc)){
    if(zs[i,]/40>=0.1){
      print("please drop me!")
      print(rownames(zs)[i])
    }else{
      # print("nothing to print")
    }
  
}


drops_20 <- c("mz693.1_rt0.2653", "mz209.8_rt0.3113","mz323.1_rt0.7153","mz911.5_rt10.49")
