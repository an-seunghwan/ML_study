
# Rmpi code upgrade 1

### 1. spawn slaves


```R
library(Rmpi)
# 전체 mpi 과정의 걸리는 시간을 확인하기 위함
start_time = Sys.time()
mpi.spawn.Rslaves(nslaves = mpi.universe.size() - 1)
```

### 2. define global parameters(variables)


```R
# global parameters
beta = rep(0,10)
gamma = 1

# send global variables at once(using list data type)
global_obj = list(beta, gamma)
mpi.bcast.Robj2slave(global_obj)
```

* 이때 global parameters를 일일이 broadcasting하지 않고, list 타입으로 묶어 한번에 broadcasting을 한다.
* 그리고 나중에 list indexing을 이용해 slave들이 parameter들을 읽을 수 있도록 한다.

### 3. broadcast commands


```R
# 각 slave마다 각 server에 저장되어 있는 자신의 rank에 맞는 data를 읽도록 명령 
mpi.bcast.cmd(cmd = assign(paste("new_testdata", mpi.comm.rank(), sep=""), as.matrix(read.csv(file = paste("/home/user/new_data", mpi.comm.rank(), "_mpi.csv", sep=""), header = T))))
mpi.bcast.cmd(cmd = source("reg_fun.R"))
```

* ```mpi.bcast.cmd(cmd = source("reg_fun.R"))```의 역할
    + 각 slave마다 다른 함수를 사용해야 하는 경우가 생길 수 있다.
    + 따라서 그런 경우에 대비하여 각 server마다 미리 slavefunction_reg을 정의한 reg_fun.R이라는 파일을 source 함수를 이용해 읽도록 함으로써 각 server의 slave들이 각자의 function을 수행할 수 있게 된다.
    + 또한 reg_fun.R 파일에 미리 필요한 변수나 함수를 정의하여 master에서 broadcasting을 굳이 여러번 하지 않도록 한다.
    + reg_fun.R 파일에 저장된 코드는 다음 코드와 같다.


```R
slavefunction_reg = function(data, beta)
{
  # runtime manage(less than 3sec)
  slave_start_time = as.numeric(format(Sys.time(), "%s"))
  x = data[,1:10]
  y = data[,11]
  for(i in 1:10)
  {
    gradient = t(x) %*% (y - x%*%beta)/length(y)
    beta = beta + global_obj[[2]] * gradient
    # runtime이 3초를 초과하면 바로 beta값을 return하도록 한다
    if((as.numeric(format(Sys.time(), "%s")) - slave_start_time) > 3) 
      return(beta)
  }
  return(beta)
}
```

### 4. communication and parameter update part


```R
i = 1
while(i <= 100)
{
  # 각 server마다 별도로 정의되어 있는 slavefunction_reg를 수행하도록 함
  beta_list = mpi.remote.exec(cmd = slavefunction_reg(get(paste("new_testdata", mpi.comm.rank(), sep="")), global_obj[[1]]))
  beta_vec = matrix(0, 10, 1)
  for(k in 1 : (mpi.comm.size() - 1)){
    beta_vec = beta_vec + as.matrix(beta_list[[k]])
  }
  beta = beta_vec / (mpi.comm.size() - 1)
  # update global parameters
  global_obj[[1]] = beta
  # update된 global parameter를 다시 broadcasting
  mpi.bcast.Robj2slave(global_obj)
  i = i + 1
}
```

### 5. save the result and mpi.quit


```R
# 결과를 저장
write.csv(beta, file = "/home/user/beta_result3.csv",row.names = F)
# slave를 close하고 mpi를 종료
mpi.close.Rslaves()
end_time = Sys.time()
print(end_time - start_time)
mpi.quit()
```
