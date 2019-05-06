
# Flexible MPI using send and recv

* 반드시 recv를 실행하기 전, send를 실행하여 message passing을 활성화 시킨다.
    - 만약 그렇지 않으면, recv가 무한 루프에 빠진 것처럼 R의 정상적 종료가 불가능해진다.


```R
if(!require(Rmpi)) install.packages("Rmpi")
require(Rmpi)

mpi.spawn.Rslaves(nslaves = 2)

# (tag별로 수행해야 할 동작을 다르게  설정하면 더욱 유연한 slave function 가능)
slavefunction = function()
{
  task = mpi.recv.Robj(0, mpi.any.tag()) # master로 부터 임의의 tag가 달린 task를 받는다
  temp = task * 2
  mpi.send.Robj(temp, 0, 1) # master(rank = 0) 에게 tag가 1이달린 temp 결과를 send
}

mpi.bcast.Robj2slave(slavefunction)
mpi.bcast.cmd(slavefunction())

# 처리해야 할 데이터를 list 형식으로 만든다
task = list()
for(i in 1:10) {
  task[[i]] = i
}

results = vector("list", 2)
# task 데이터를 역순으로 slave 1, 2에게 먼저 send
# (꼭 역순으로 하지 않는 방법 필요)
# (이를 한번에 실행시킬 수 있는 방법 필요 - for, while, ...)
mpi.send.Robj(task[[length(task)]], 1, 1)
mpi.send.Robj(task[[length(task)]], 2, 1)

while(length(task) >= 1) # 처리해야 할 task가 남지 않을 때까지 실행
{
  task[[length(task)]] = NULL # 처리한 task의 부분은 NULL 처리
  mpi.bcast.cmd(slavefunction())
        # 앞에서 send를 실행하여 slave에게 task 데이터를 보냈으므로
        # slave function을 수행하도록 하여 앞에서 task 데이터의 처리를 명령
  message = mpi.recv.Robj(mpi.any.source(), mpi.any.tag())
        # slave로 부터 처리된 task 데이터를 받는다
        # 임의의 slave와 임의의 tag
        # (이때 slave의 message가 겹치면 둘 중 하나가 덮어씌어지는 것처럼 무시되어 사라진다. - 해결 필요*****)
        # (각각의 message를 slave 별로 저장할 수 있도록 하는 방법 필요)
  message_info = mpi.get.sourcetag() # message의 정보를 얻는 과정
  slave_id = message_info[1]
  # tag = message_info[2]
  results[[slave_id]] = c(results[[slave_id]], message)
  if(length(task) != 0) mpi.send.Robj(task[[length(task)]], slave_id, 1)
        # 처리해야 할 task가 남았으면 다시 slave 들에게 task 데이터 send 진행
}
print(results)
mpi.close.Rslaves()
mpi.quit()
```

    Loading required package: Rmpi
    Warning message:
    "package 'Rmpi' was built under R version 3.5.2"


    Error in mpi.comm.spawn(slave = system.file("Rslaves64.cmd", package = "Rmpi"), : Internal MPI error!, error stack:
    MPI_Comm_spawn(cmd="C:/Users/dpelt/R_personal/Rmpi/Rslaves64.cmd", argv=0x0000000017E2B690, maxprocs=2, MPI_INFO_NULL, root=0, MPI_COMM_SELF, intercomm=0x000000001744FAD8, errors=0x0000000019EF4E78) failed
    Internal MPI error!  FAILspawn not supported without process manager
    Traceback:
    

    1. mpi.spawn.Rslaves(nslaves = 2)

    2. mpi.comm.spawn(slave = system.file("Rslaves64.cmd", package = "Rmpi"), 
     .     slavearg = arg, nslaves = nslaves, info = 0, root = root, 
     .     intercomm = intercomm, quiet = quiet)

