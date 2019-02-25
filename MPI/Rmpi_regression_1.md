
# Rmpi regression example

### 1. R code : regression


```R
library(Rmpi)
start_time = Sys.time()
mpi.spawn.Rslaves(nslaves = mpi.universe.size() - 1)

.Last <- function(){
  if (is.loaded("mpi_initialize")){
    if (mpi.comm.size(1) > 0){
      print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
}

beta = rep(0,10)
gamma = 1

slavefunction_reg = function(data, beta)
{
  x = data[,1:10]
  y = data[,11]
  for(i in 1:10)
  {
    gradient = t(x) %*% (y - x%*%beta)/length(y)
    beta = beta + gamma * gradient
  }
  return(beta)
}

mpi.bcast.Robj2slave(slavefunction_reg)
mpi.bcast.Robj2slave(beta)
mpi.bcast.Robj2slave(gamma)

mpi.bcast.cmd(cmd = assign(paste("new_testdata", mpi.comm.rank(), sep=""), as.matrix(read.csv(file = paste("/home/user/new_data", mpi.comm.rank(), "_mpi.csv", sep=""), header = T))))
i = 1
while(i <= 100)
{
  beta_list = mpi.remote.exec(cmd = slavefunction_reg(get(paste("new_testdata", mpi.comm.rank(), sep="")), beta))
  beta_vec = matrix(0, 10, 1)
  for(k in 1 : (mpi.comm.size() - 1)){
    beta_vec = beta_vec + as.matrix(beta_list[[k]])
  }
  beta = beta_vec / (mpi.comm.size() - 1)
  # print(beta)
  mpi.bcast.Robj2slave(beta)
  i = i + 1
}

write.csv(beta, file = "/home/user/beta_result.csv",row.names = F)

mpi.close.Rslaves()
end_time = Sys.time()
print(end_time - start_time)
mpi.quit()
```

### 2. sample data generation code


```R
# 각 slave의 rank에 맞게 for문 안의 제어문을 수정하여 데이터를 생성한다
for(i in 1:4){
  x = matrix(rnorm(1000000), 100000, 10)
  y = x %*% c(1,0,0,0,0,0,0,0,0,0) + rnorm(100000)
  data = cbind(x, y)
  write.csv(data, file = paste("/home/user/new_data",i,"_mpi.csv", sep=""),row.names = F)
}
```

### 3. console 결과 창

![console image](https://github.com/SeungHwan-AN/ML_study/blob/master/MPI/rmpi_regression_console.png)

* 1회 통신(send or receive)에 0.15초 정도 걸리는 것을 확인할 수 있다
    - R code에서 통신이 이루어지는 while문이 10번 반복하기 때문

### 4. regression result csv file

![csv file](https://github.com/SeungHwan-AN/ML_study/blob/master/MPI/rmpi_regression_result.PNG)
