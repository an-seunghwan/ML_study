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
	
##### 3. SSH configuration 설정
* 참고(youtube : https://www.youtube.com/watch?v=VcLEdqrvPKI)
1. [master 서버]
	* sudo apt-get install openssh-server
	* (password 입력)
	* ssh (slave IP 주소) -> yes
	* (password 입력)
	* 이때 slave IP주소로 이동한 것을 볼 수 있다
	* exit
	* cd .ssh
	* 이때 known_hosts 파일이 있다면 성공
2. [slave 서버]
	* ssh (master IP 주소) -> yes
	* 이때 master IP주소로 이동한 것을 볼 수 있다
	* cd .ssh
	* 이때 known_hosts 파일이 있다면 성공
3. [master 서버]
	* ssh (서버이름) -> yes
	* cd .ssh
	* ssh-keygen -t rsa
	* enter 키 3번 누른다
	* ls 
		- id_rsa id_rsa.pub known_hosts 가 생성되었으면 성공
	* scp id_rsa.pub user@(slave IP 주소):/.ssh
	* (password 입력)
	* exit
	* ssh (slave IP주소)
	* (password 입력)
	* cd .ssh
	* ls
		- id_rsa.pub known_hosts 가 생성되었으면 성공
	* cat id_rsa.pub >> authorized_keys
	* ls
		- authorized_keys id_rsa.pub known_hosts 가 생성되었으면 성공
	* exit
4. [최종확인]
	* ~/.ssh$ ssh (slave IP주소) 를 실행하였을 때 비밀번호 입력없이 서버를 이동하였다면 ssh 설정 완료

##### 4. R package 설치
* Rmpi, foreach, doMPI package를 차례로 설치해준다


	* cannot find mpi.h header file 에러가 발생하는 경우
		- sudo apt-get install libopenmpi-dev
		- 설치 후 다시 Rmpi package 설치


##### [현재상황]
* stathdsf05와 stat04 서버는 ssh 통신이 가능
* 그러나 doMPI를 활용한 clustering 작업이 불가능


   
