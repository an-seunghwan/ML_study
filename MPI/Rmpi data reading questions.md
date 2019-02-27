
# Rmpi의 데이터 broadcasting 방식

* Rmpi 코드에서 다음과 같은 부분이 있다


```R
mpi.bcast.cmd(cmd = assign(paste("new_testdata", mpi.comm.rank(), sep=""), as.matrix(read.csv(file = paste("/home/user/new_data", mpi.comm.rank(), "_mpi.csv", sep=""), header = T))))
```

* 질문 : server마다 필요한 데이터가 같은 경로에 저장되어 있으므로, master에서 이 전부를 읽고, 각 slave별로 broadcasting해주는 것이 아닐까?
    1. 기본적인 setting은 server1, 2, 3에 각각 1 ~ 4번, 5 ~ 8번, 9 ~ 12번의 데이터 파일이 저장되어 있다.
    2. 모든 데이터 파일의 경로는 /home/user로 동일하다
    
* 이때 만약 master가 데이터 파일 전부를 읽고 broadcasting을 하는 방식이라면, single machine으로 읽고 처리하기 힘든 매우 큰 데이터를 분산하여 처리하려는 우리의 목적과 맞지 않게 된다!!!

* 결론 : 그렇지 않다!
    1. 우선 ```mpi.bcast.cmd```에서 ```cmd```옵션만을 이용해 command만 broadcasting해주었다.
    2. 데이터를 읽는 시간의 차이가 발생한다.
        + 1~12번을 모두 합친 데이터를 한번에 읽는 경우(single machine) : Time difference of 52.67144 secs
        + ```mpi.bcast.cmd```를 이용해 각 slave별로 데이터 읽기를 수행하는 경우 : Time difference of 0.4392927 secs
    3. server1에서 server3에 저장된 데이터를 읽으려고 한 경우 : error 발생
        + ```read.csv(file = paste("/home/user/new_data", 10, "_mpi.csv", sep=""), header = T)```
        + Error in file(file, "rt") : cannot open the connection
          In addition: Warning message:
          In file(file, "rt") :
          cannot open file '/home/user/new_data10_mpi.csv': No such file or directory
