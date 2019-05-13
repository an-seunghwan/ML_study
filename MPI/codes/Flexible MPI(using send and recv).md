
# Flexible MPI
## using send and recv


```R
if(!require(Rmpi)) install.packages('Rmpi')
require(Rmpi)

mpi.spawn.Rslaves(nslaves = 3)
```

    Loading required package: Rmpi
    Warning message:
    "package 'Rmpi' was built under R version 3.5.2"


    Error in mpi.comm.spawn(slave = system.file("Rslaves64.cmd", package = "Rmpi"), : Internal MPI error!, error stack:
    MPI_Comm_spawn(cmd="C:/Users/dpelt/R_personal/Rmpi/Rslaves64.cmd", argv=0x00000000197C83E0, maxprocs=3, MPI_INFO_NULL, root=0, MPI_COMM_SELF, intercomm=0x0000000016E58758, errors=0x000000001763E838) failed
    Internal MPI error!  FAILspawn not supported without process manager
    Traceback:
    

    1. mpi.spawn.Rslaves(nslaves = 3)

    2. mpi.comm.spawn(slave = system.file("Rslaves64.cmd", package = "Rmpi"), 
     .     slavearg = arg, nslaves = nslaves, info = 0, root = root, 
     .     intercomm = intercomm, quiet = quiet)



```R
### comparison of each slave's message arriving time by giving different tasks ###
# (computation cost)
slavefunction = function() {
  message = mpi.recv.Robj(source = 0, tag = mpi.any.tag())
  message_info = mpi.get.sourcetag()
  tag = message_info[2]
  if(tag == 1) {
    temp = sum(task[1:10])
    mpi.send.Robj(temp, dest = 0, tag = 1)
  } else if(tag == 2) {
    temp = sum(task[1:1000])
    mpi.send.Robj(temp, dest = 0, tag = 1)
  } else if(tag == 3) {
    max = task[1]
    for(i in 2:n) {
        if(task[i] > max) max = task[i]
    }
    mpi.send.Robj(max, dest = 0, tag = 1)
  }
}

mpi.bcast.Robj2slave(slavefunction)
```


```R
### comparison of each slave's message arriving time using time consuming while loop ###
### -> result : order of slave 3 -> slave 2 -> slave 1 ###
# slavefunction = function() {
#   message = mpi.recv.Robj(source = 0, tag = mpi.any.tag())
#   message_info = mpi.get.sourcetag()
#   tag = message_info[2]
#   if(tag == 1) {
#     temp = sum(task[1:n])
#     time <- Sys.time()
#     while((as.numeric(Sys.time()) - as.numeric(time)) < 10^5) {}
#     mpi.send.Robj(temp, dest = 0, tag = 1)
#   } else if(tag == 2) {
#     temp = sum(task[1:n])
#     time <- Sys.time()
#     while((as.numeric(Sys.time()) - as.numeric(time)) < 10) {}
#     mpi.send.Robj(temp, dest = 0, tag = 1)
#   } else if(tag == 3) {
#     temp = sum(task[1:n])
#     time <- Sys.time()
#     while((as.numeric(Sys.time()) - as.numeric(time)) < 2) {}
#     mpi.send.Robj(max, dest = 0, tag = 1)
#   }
# }
```


```R
set.seed(1)
n = 10^4
task = rpois(n, lambda = 10)
mpi.bcast.Robj2slave(n)
mpi.bcast.Robj2slave(task)

# wake slaves up
junk = 0
# slave들에게 Robj를 동시에 보낼 수 있는 방법 필요*****
mpi.send.Robj(junk, 1, 1); mpi.send.Robj(junk, 2, 2); mpi.send.Robj(junk, 3, 3)
```


```R
s = 2 # set partial barrier
message = c()
message_info = matrix(numeric(0), 0, 2)

mpi.bcast.cmd(slavefunction()) # get message from slaves

for(i in 1:s) # the number of message = partial barrier(s)
{
  temp_message = mpi.recv.Robj(mpi.any.source(), mpi.any.tag())
  message = c(message, temp_message)
  info = mpi.get.sourcetag()
  message_info = rbind(message_info, c(info[1], info[2])) # save slave rank & tag
}
# 이때 for문을 빠져나오면 받지않은 나머지 1개의 message를 잃어버리는 문제 발생*****
(message)
(message_info)
```


```R
mpi.close.Rslaves()
```


```R
### using message length -> not working ###
# message = c()
# message_info = matrix(numeric(0), s)
# while(length(message) < s)
# {
#   temp_message = mpi.recv.Robj(mpi.any.source(), mpi.any.tag())
#   message = c(message, temp_message)
#   info = mpi.get.sourcetag()
#   message_info = cbind(message_info, c(info[1], info[2])) # slave rank & tag
# }
# (message)
# (message_info)
```
