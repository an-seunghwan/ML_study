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

beta = rep(0,3)
gamma = 1

slavefunction_reg = function(data, beta)
{
  x = data[,1:3]
  y = data[,4]
  for(i in 1:1000)
  {
    gradient = t(x) %*% (y - x%*%beta)/length(y)
    beta = beta + gamma * gradient
  }
  return(beta)
}

mpi.bcast.Robj2slave(slavefunction_reg)
mpi.bcast.Robj2slave(beta)
mpi.bcast.Robj2slave(gamma)

mpi.bcast.cmd(cmd = assign(paste("testdata", mpi.comm.rank(), sep=""), as.matrix(read.csv(file = paste("/home/user/data", mpi.comm.rank(), "_mpi.csv", sep=""), header = T))))
i = 1
while(i <= 10)
{
  beta_list = mpi.remote.exec(cmd = slavefunction_reg(get(paste("testdata", mpi.comm.rank(), sep="")), beta))
  beta_vec = matrix(0, 3, 1)
  for(k in 1 : (mpi.comm.size() - 1)){
    beta_vec = beta_vec + as.matrix(beta_list[[k]])
  }
  beta = beta_vec / (mpi.comm.size() - 1)
  # print(beta)
  mpi.bcast.Robj2slave(beta)
  i = i + 1
}

write.csv(beta, file = "/home/user/beta_result2.csv",row.names = F)

mpi.close.Rslaves()
end_time = Sys.time()
print(end_time - start_time)
mpi.quit()