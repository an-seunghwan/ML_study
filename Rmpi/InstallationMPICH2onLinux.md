# Installation MPICH2 on Linux

#### 별도의 설명이 없으면 아래 명령어는 모두 cmd command임 
#### 설치시 파일의 경로에 대한 설정이 매우 중요하다!

##### 1. MPICH 설치 과정
* MPICH 홈페이지에서 mpich2 tar file을 다운로드 (home/dpeltms79/MPI folder에 있다고 가정)
* cd /home/dpeltms79/MPI
	+ (mkdir MPI)
* tar xfz mpich-3.3.tar.gz
	+ MPI directory에 sub-directory mpich-3.3가 생성됨
* mkdir mpich-install
	+ 이는 installation directory를 설정
* mkdir mpich-build
	+ 이는 build directory를 설정
* export RSHCOMMAND=ssh
* /home/dpeltms79/MPI/mpich-3.3/configure --enable-debuginfo --enable-fast=03 --enable-shared --with-pm=hydra:gforker:remshell --prefix=/home/dpeltms79/MPI/mpich-install |& tee c.txt
	+ configure를 설정하는 과정
* cd /home/dpeltms79/MPI/mpich-build
	+ build directory로 이동
* make 2>&1 | tee m.txt
	+ build MPICH2
* cd
	+ go to root directory
* make install |& tee mi.txt
	+ install MPICH
* export PATH=/home/dpeltms79/MPI/mpich-install/bin:$PATH
	+ add the bin directory to your path
* which mpd & mpicc & mpiexec & mpirun
	+ 제대로 설치가 되었는지 확인

##### 2. cluster 환경 설정
	* [중요!!!] host name을 설정하는 과정에서 이름을 직접 변경해 주어야한다. 
	* 이때 필요한 directory를 직접 만들어서 변경을 한다
	* 또한 sudo를 이용하여 command를 작성한다
* sudo mkdir /etc/sysconfig
* sudo vi /etc/sysconfig/network
	1) insert모드로 전환 후 'localhost' 입력
	2) 다음 wq를 이용해 저장한다
* sudo vi /etc/hosts
	1) 다음과 같이 수정한다
		+ 127.0.0.1		localhost
		+ 127.0.1.1		stathdsf05
		+ 172.16.207.139  stat04
	2) wq를 이용해 저장
* sudo vi mpi.hosts    
	1) 각 줄마다 machine의 이름을 입력한다
    	+ localhost
    	+ stathdsf05
    	+ stat04
    2) wq를 이용해 저장
* touch .mpd.conf
* chmod 600 .mpd.conf
* vi .mpd.conf
	+ MPD_SECRETWORD=<123456qwer!> 입력
	+ wq를 이용해 저장

##### 3. R package 설치
* Rmpi, foreach, doMPI package를 차례로 설치해준다


	* cannot find mpi.h header file 에러가 발생하는 경우
		- sudo apt-get install libopenmpi-dev
		- 설치 후 다시 Rmpi package 설치


* ssh configuration을 수행한다(youtube : https://www.youtube.com/watch?v=VcLEdqrvPKI)

##### [현재상황]
* stathdsf05와 stat04 서버는 ssh 통신이 가능
* 그러나 doMPI를 활용한 clustering 작업이 불가능


   