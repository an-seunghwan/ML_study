# openMPI with R

	본 설치는 ubuntu 16.04 버전에서 설치되었음
    또한 sudo 계정을 이용한 설치를 권장함(Rmpi를 사용하고자 하는 모든 컴퓨터의 local 아이디와 비밀번호 통일)
    ※ 특별한 경우가 아니면 모든 cluster내의 컴퓨터의 설치, 파일 저장, 라이브러리 경로를 동일하게 맞출 것

### 1. openMPI 설치
1. www.open-mpi.org 사이트에서 최신 release를 다운로드
	+ 저는 openmpi-4.0.0.tar.gz 파일을 받았습니다
2. 이를 가장 기본 경로에 저장
	+ /home/user/ 경로에 저장하였습니다
3. 파일의 압축을 풉니다
```
tar -xzf openmpi-4.0.0.tar.gz
```
4. /home/user/openmpi-4.0.0 와 같은 경로가 생성
5. 경로를 바꿉니다
```
cd openmpi-4.0.0
```
6. 설치의 설정을 다음을 이용해 configure 합니다. home 경로에 설치하시려면 sudo를 이용하지 않아도 됩니다
	+ ※ 이때 반드시 설치경로인 prefix의 경로는 cluster 내의 모든 컴퓨터가 동일해야 합니다
```
sudo ./configure --prefix="/usr/local/openmpi"
```
7. compile 과정
```
sudo make
```
8. install 과정
```
sudo make install
```
9. 경로 설정
```
export PATH="$PATH:/usr/local/openmpi/bin"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/openmpi/lib"
```
10. 필요 파일 설치
```
sudo apt-get install openmpi-bin openmpi-common openssh-client openssh-server libopenmpi1.10 libopenmpi-dev
```
11. openMPI 설치가 제대로 되었는지 확인한다
	+ hello.c 파일은 root 경로에 저장되어있다고 가정
	+ hello.c 파일은 본 github에서 다운 가능
```
mpicc -o hello hello.c
mpirun -np 2 ./hello
```
```
(결과)
user@statushdfs04:~$ mpicc -o hello hello.c                                     
user@statushdfs04:~$ mpirun -np 2 ./hello
Hello world from processor statushdfs04, rank 0 out of 2 processors
Hello world from processor statushdfs04, rank 1 out of 2 processors
```

### 2. cluster 환경 설정
	※ host name을 설정하는 과정에서 이름을 직접 설정해 주어야한다. 
1. cluster 내의 컴퓨터들의 IP주소 설정
```
sudo vi /etc/hosts
```
* 다음과 같이 수정한다
	+ 127.0.0.1		localhost
	+ 127.0.1.1		com1
	+ [IP주소]      com2
	+ ...
2. sudo vi host_file    
	* 각 줄마다 cluster내의 컴퓨터의 이름들을 입력한다(node 컴퓨터)
    	+ com1 
    	+ com2
    	+ com3
    	+ ...
    * master 컴퓨터의 경우 다음과 같이 설정한다
    	+ com1 slots = 5
    	+ com2 slots = 4
    	+ com3 slots = 4
    	+ ...
    		- 이때 slots 옵션은 각 컴퓨터 node마다 최대로 spawn할 slave의 개수를 말한다
    		- master가 작동하는 com1에 대해서는 master용 cpu 1개를 남기기위해 5개로 설정
3. ```touch .mpd.conf```
4. ```chmod 600 .mpd.conf```

### 3. SSH configuration 설정
	참고(youtube : https://www.youtube.com/watch?v=VcLEdqrvPKI)
1. [master 서버]
	* ```sudo apt-get install openssh-server```
	* ```ssh (slave IP 주소)``` - yes
		- 이때 slave IP주소로 이동한 것을 볼 수 있다
	* ```exit```
	* ```cd .ssh```
	* ```ls```
		- known_hosts 파일이 있다면 성공
2. [slave 서버]
	* ```ssh (master IP 주소)``` - yes
		- 이때 master IP주소로 이동한 것을 볼 수 있다
	* ```cd .ssh```
	* ```ls```
		- 이때 known_hosts 파일이 있다면 성공
3. [master 서버]
	* ```ssh (master 서버이름)```
	* ```cd .ssh```
	* ```ssh-keygen -t rsa```
		+ master 컴퓨터에서 한 번만 수행
	* ```enter 키 3번 누른다```
	* ```ls```
		- id_rsa id_rsa.pub known_hosts 가 생성되었으면 성공
	* ```scp id_rsa.pub user@(slave IP 주소):~/.ssh```
	* ```exit```
	* ```ssh (slave IP주소)```
	* ```cd .ssh```
	* ```ls```
		- id_rsa.pub known_hosts 가 생성되었으면 성공
	* ```cat id_rsa.pub >> authorized_keys```
	* ```ls```
		- authorized_keys id_rsa.pub known_hosts 가 생성되었으면 성공
	* ```exit```
4. [최종확인]
	* ```ssh (slave IP주소)``` 를 실행하였을 때 비밀번호 입력없이 서버를 이동하였다면 ssh 설정 완료

### 4. openMPI의 cluster 작동 여부 확인

1. openMPI 설치가 제대로 되었는지 확인한다
	+ master 컴퓨터에서 작동
```
user@stathdfs05:~$ mpicc -o hello hello.c
user@stathdfs05:~$ mpirun -n 13 -hostfile host_file ./hello
```
```
(결과)
Hello world from processor stathdfs05, rank 0 out of 13 processors
Hello world from processor stathdfs05, rank 1 out of 13 processors
Hello world from processor stathdfs05, rank 2 out of 13 processors
Hello world from processor stathdfs05, rank 3 out of 13 processors
Hello world from processor stathdfs05, rank 4 out of 13 processors
Hello world from processor stat04, rank 5 out of 13 processors
Hello world from processor stat04, rank 7 out of 13 processors
Hello world from processor stat04, rank 6 out of 13 processors
Hello world from processor stat04, rank 8 out of 13 processors
Hello world from processor statushdfs04, rank 9 out of 13 processors
Hello world from processor statushdfs04, rank 10 out of 13 processors
Hello world from processor statushdfs04, rank 11 out of 13 processors
Hello world from processor statushdfs04, rank 12 out of 13 processors
```

### 5. R을 이용한 MPI 활용

1. R을 3.5.2 버전으로 반드시 업데이트
2. ```.libPaths()```를 이용해 모든 컴퓨터 node의 library 경로가 같은지 확인
3. Rmpi package를 설치
```
install.packages("Rmpi")
```
4. R을 이용한 MPI 정상작동 확인
	+ rmpi_test.R 파일은 모든 컴퓨터의 root 경로에 저장되어있다고 가정
	+ rmpi_test.R 파일은 본 github에서 다운 가능
```
user@stathdfs05:~$ mpirun -n 1 -hostfile host_file R --slave -f /home/user/rmpi_test.R
        12 slaves are spawned successfully. 0 failed.
master  (rank 0 , comm 1) of size 13 is running on: stathdfs05
slave1  (rank 1 , comm 1) of size 13 is running on: stathdfs05
slave2  (rank 2 , comm 1) of size 13 is running on: stathdfs05
slave3  (rank 3 , comm 1) of size 13 is running on: stathdfs05
... ... ...
slave11 (rank 11, comm 1) of size 13 is running on: statushdfs04
slave12 (rank 12, comm 1) of size 13 is running on: statushdfs04
$slave1
[1] "I am 1 of 13 stathdfs05"

$slave2
[1] "I am 2 of 13 stathdfs05"

$slave3
[1] "I am 3 of 13 stathdfs05"

$slave4
[1] "I am 4 of 13 stathdfs05"

$slave5
[1] "I am 5 of 13 stat04"

$slave6
[1] "I am 6 of 13 stat04"

$slave7
[1] "I am 7 of 13 stat04"

$slave8
[1] "I am 8 of 13 stat04"

$slave9
[1] "I am 9 of 13 statushdfs04"

$slave10
[1] "I am 10 of 13 statushdfs04"

$slave11
[1] "I am 11 of 13 statushdfs04"

$slave12
[1] "I am 12 of 13 statushdfs04"

[1] 1
```
5. R코드에서의 mpi.universe.size가 host_file에서 설정한 slots의 수
6. command에서 -n 1은 실행시키고자 하는 코드(프로그램)을 master & slave 그룹을 이용해 몇번 반복할 것인지 지정한다
	* 따라서 1로 값을 설정해주어 프로그램을  1번만 수행하고, 이때 그룹은 master 1명과 그외 slave들로 이루어진다

### 자주 발생하는 error 대처법

1. sudo : unable to resolve host [hostname]
	* hostname과 sudo vi /etc/hostname의 이름을 동일하게 수정해준다
		* hostname : hostname 확인 명령어
	* hostname, hosts, host_file 들의 이름과 IP주소가 다른 것이 있는지 확인한다
2. Error in library("Rmpi") : there is no package called ‘Rmpi’
	* export LD_LIBRARY_PATH=/home/user/R/x86_64-pc-linux-gnu-library/3.5:$LD_LIBRARY_PATH
	* 즉, 모든 컴퓨터에 library가 동일한 경로로 설치되어 있어야 한다
3. 특정 노드에서 slave spawn이 안됨
	* 2번처럼 모든 node의 Rmpi package가 설치된 library가 일치되어야 한다.
	* .libPath()를 이용해서 모든 node의 Rmpi library path를 모두 동일하게 설정

