# 각 slave의 rank에 맞게 for문 안의 제어문을 수정하여 데이터를 생성한다
for(i in 1:4){
  x = matrix(rnorm(1000000), 100000, 10)
  y = x %*% c(1,0,0,0,0,0,0,0,0,0) + rnorm(100000)
  data = cbind(x, y)
  write.csv(data, file = paste("/home/user/new_data",i,"_mpi.csv", sep=""),row.names = F)
}
