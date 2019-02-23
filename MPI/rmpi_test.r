library(Rmpi)

##### 1
# Clean up if R quits unexpectedly
.Last = function(){
  if(is.loaded("mpi_initialize")){
    if(mpi.comm.size(1) > 0){
      ;
      print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
}
mpi.spawn.Rslaves(nslaves = mpi.universe.size() - 1)
mpi.remote.exec(paste("I am", mpi.comm.rank(),"of", mpi.comm.size(), Sys.info()["nodename"]))
mpi.close.Rslaves()
mpi.quit()
