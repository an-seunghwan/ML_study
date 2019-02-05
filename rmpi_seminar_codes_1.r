rm(list = ls())
gc()

library(Rmpi)

##### 1
# Clean up if R quits unexpectedly
.Last = function(){
  if(is.loaded("mpi_initialize")){
    if(mpi.comm.size(1) > 0){
      print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
}

##### 2. Spawn slaves
# Spawn as many slaves as possible
mpi.spawn.Rslaves()
# Spawn two slaves using nslaves
mpi.spawn.Rslaves(nslaves = 2)
# Let each slave return their ranks
mpi.remote.exec(paste("I am", mpi.comm.rank(),"of", mpi.comm.size()))

# $slave1
# [1] "I am 1 of 3"
# $slave2
# [1] "I am 2 of 3"

mpi.close.Rslaves()

##### 2.1
mpi.spawn.Rslaves()
mpi.remote.exec(sum(1:mpi.comm.rank()))

##### 2.2 Measure time to compute
mpi.spawn.Rslaves(nslaves=4) # changing nslaves
ptm<-proc.time() 
mpi.iparReplicate(n = 400, expr = mean(rnorm(1000000)))
print(proc.time() - ptm)
mpi.close.Rslaves()

##### 3. Sending Data

# 1. Sending Identical Data to Slave CPUs

# ex) Send a character string to 3 slaves
mpi.spawn.Rslaves(nslaves = 3)
mpi.bcast.cmd(x <- mpi.bcast.Robj()) # to send an object to each slave(use <-)
# slaves will execute some functions with 'x' as master variable(?)
# cmd is not clarified
x = c("This is a test")
mpi.bcast.Robj(x)
mpi.remote.exec(x)
# $slave1
# [1] "This is a test."
# $slave2
# [1] "This is a test."
# $slave3
# [1] "This is a test."

mpi.close.Rslaves()

# 2. Sending Non-identical Data to Slave CPUs
# This is useful when we have a large set of data.
# First, divide the data and put them in a list.
# Second, send the list to each slave.

# We need to have exactly the same number of object portions as the number of the CPUs.
# If not, we will get an error message.
# For example, if you spawn 3 slave CPUs, and your object to send is a list of 4,
# mpi.scatter.Robj2slave() will not work because you don't have the equal numbers.
# However, mpi.scatter.Rbj() will work because master plus slaves equal to the number of list in the object.

# ex1) Split and send an object of list of 4 to master and slaves
mpi.spawn.Rslaves(nslaves = 3)
x = list("master", rnorm(3), letters[2], 1:10) # has 4 portions
mpi.bcast.cmd(x <- mpi.scatter.Robj())
mpi.scatter.Robj(x)
mpi.remote.exec(x)
# return(print) a splitted list x
# the master and each slaves print one portion of list x respectively

mpi.close.Rslaves()
# Note that master receive object x but it does not overwrite existing x.

# ex2) Divide 8*4 matrix into 4 blocks and send to slaves.
# 1. Spawn 4 slaves CPUs, and create a 8*4 matrx with random numbers.
mpi.spawn.Rslaves(nslaves = 4)
mat = matrix(rnorm(32), 8)
# 2. Split the matrix into 4 of 2*4 matrices
smat = lapply(.splitIndices(nrow(mat), 4), function(i) mat[i,])
# 3. Send each matrix to slave CPUs
mpi.scatter.Robj2slave(smat) # 특별히 수행할 cmd가 없으면 mpi.bcast.cmd 필요 없음...?
mpi.remote.exec(smat)
# return(print) a splitted matrix smat
# the master and each slaves print one portion of matrix smat respectively

mpi.close.Rslaves()

##### 4. Getting Data Back From Slaves -> 예제 체크 필요!!!

# 1. mpi.gather.Robj()
# This will retrieve data from each slave and put them together.

# ex) Getting slave number information from each slave
mpi.spawn.Rslaves(nslaves = 3)
mpi.bcast.cmd(id <- mpi.comm.rank())
mpi.bcast.cmd(x <- paste("I am no.", id))
mpi.bcast.cmd(mpi.gather.Robj(x))
x = "I am a master"
mpi.gather.Robj(x)
# [1] "I am a master"    "I am slave no. 1" 
# "I am slave no. 2" "I am slave no. 3"

mpi.remote.exec(x)
# $slave1
# [1] "I am slave no. 1"
# $slave2
# [1] "I am slave no. 2"
# $slave3
# [1] "I am slave no. 3"

mpi.close.Rslaves()

# 2. mpi.allgather.Robj()
# This will retrieve data from each slave and give the whole data to all slaves.

# ex) Retrieve a string data from each slave and send the all data to master and all slaves.
mpi.spawn.Rslaves(nslaves = 3)
x = c("fruits", "apple", "banana", "orange")
mpi.bcast.cmd(x <- mpi.scatter.Robj())
mpi.scatter.Robj(x)
mpi.remote.exec(x)
# return resuls

mpi.bcast.cmd(x <- mpi.allgather.Robj(x))
mpi.allgather.Robj(x)
# [1] "fruits" "apple" "banana" "orange" 

mpi.remote.exec(x)
# $slave1 
# [1] "fruits" "apple" "banana" "orange" 
# $slave2 
# [1] "fruits" "apple" "banana" "orange" 
# $slave3 
# [1] "fruits" "apple" "banana" "orange"

mpi.close.Rslaves()

# 3. mpi.reduce(), mpi.allreduce()
# : Reduce data by simple operation

# ex) set a value of x to 1 in the master, 2, 3, 4 in slave 1, 2, 3
# Then using mpi.reduce to return a sum of all x
mpi.spawn.Rslaves(nslaves = 3)
# Define function for reduction
red = function(option = "sum") {
  mpi.reduce(x, type=2, op=option)
}
mpi.bcast.Robj2slave(red)
x = c(1,2,3,4)
mpi.bcast.cmd(x <- mpi.scatter.Robj())
mpi.scatter.Robj(x)
mpi.remote.exec(x)
# return results

# call the function in slaves
mpi.remote.exec(red("sum"))
# return results

mpi.reduce(x, 2, "sum")
# output is 10

mpi.close.Rslaves()

# Two more options
# : minloc and maxloc
mpi.reduce(x,2,"maxloc")
mpi.reduce(x,2,"minloc")
# the command will return two values, the value resulting from the operation and the location of the value.
# This can be useful to find slave which provide the value.

##### 5. Example with Fibonacci Sequence

# 1. Load Rmpi and spawn the slaves.
library(Rmpi)
mpi.spawn.Rslaves(nslaves = 2)

# 2. Write any necessary functions.
fib = function(n)
{
  if(n < 1)
    stop("Input must be an integer >= 1")
  if(n == 1 | n == 2)
    1
  else 
    fib(n-1) + fib(n-2)
}

# 3. Wrap all code to make sure that the slaves should run into a worker function.
# worker function
slavefunction = function() # worker 입장에서 수행해야할 함수
{
  # protocol for sending message
  # 1 = ready to process a task, 2 = task completed, 3 = slave terminated
  # protocol for receiving message(= tag???)
  # 1 = task to do, 2 = task completed
  junk = 0 # junk R object(???)
  done = 0
  while(done != 1) {
    # send a message that ready to process a task
    mpi.send.Robj(junk, 0, 1) # slave -> master(number is 0)
    # receive the task
    task = mpi.recv.Robj(mpi.any.source(), mpi.any.tag())
    task_info = mpi.get.sourcetag()
    tag = task_info[2]
    
    if(tag == 1) {
      # task processing and send process results
      # task : calculate the fibonacci sequence
      n = task$n
      val = ifelse(n == 1 | n == 2, 1, fib(n-1) + fib(n-2))
      print(val)
      res = list(result = val, n = n)
      # send results of task and a message that the task has completed.
      mpi.send.Robj(res, 0, 2)
    }
    else if(tag == 2) {
      # master said that all tasks are terminated
      done = 1
    }
    # ignore other messages or print error log
  }
  # a message that slave was terminated was send to master
  mpi.send.Robj(junk, 0, 3)
}

# broadcast slave function and fibonacci function to slaves
# use Robj2slave instead of cmd to broadcast R object(= function)(???)
mpi.bcast.Robj2slave(slavefunction)
mpi.bcast.Robj2slave(fib)

# 4. Tell the slaves to execute their worker function.
mpi.bcast.cmd(slavefunction())

# 5. Initialize data.(Generate task)
tasks = vector("list")
for(i in 1: 25)
{
  tasks[[i]] = list(n = i)
}

junk = 0
closed_slaves = 0
n_slaves = mpi.comm.size() - 1 # one less than the maximum size
results = rep(0, 25)

# 6. Communicate with the slaves to perform the computation.
# repeat till all tasks is finished and slaves terminate
while(closed_slaves < n_slaves) { # master 입장에서 수행하는 반복문
  # receive message from slave
  message = mpi.recv.Robj(mpi.any.source(), mpi.any.tag()) # what is result of this????? -> print this
  message_info = mpi.get.sourcetag()
  slave_id = message_info[1]
  tag = message_info[2]
  
  if(tag == 1) { # when a slave is ready to operate the task
    # send one of the sending message protocol
    if(length(task) > 0) { # there are left tasks
      # send a task, and then remove it from the task list
      mpi.send.Robj(tasks[[1]], slave_id, 1)
	  # slave performs at this moment????? 
      tasks[[1]] = NULL
    }
    else { # there are no tasks anymore
      mpi.send.Robj(junk, slave_id, 2)
    }
  }
  else if(tag == 2) { # when the task is done
    # message includes result of task processing
    # save at result vector
    res = message$result
    result[message$n] = res
  }
  else if(tag == 3) {
    # receive a message that a slave process is destoryed
    closed_slaves = closed_slaves + 1
  }
}

# 7. Operate on the results in list datatype.
# get next processing with received results
#...

# 8. Halt the slaves and quit.
mpi.close.Rslaves()
mpi.quit()