rm(list=ls())
gc()
set.seed(520)

# isotonic regression via PAVA algorithm
# (Pool-Adjacent-Violators Algorithm)
# monotonic trend-filtering의 일종!

# 1. parameters
n = 100 # n = n1 + n2 + n3 + n4
n1 = 10
n2 = 40
n3 = 30
n4 = 20
y = c(runif(n1, 1, 2), runif(n2, 4, 7), runif(n3, 8, 10), runif(n4, 10, 12)) # 4 groups
z = c(rep(1, n1), rep(2, n2), rep(3, n3), rep(4, n4))
# (y, z) 변수쌍은 이미 정렬되어 있다
weight = rep(1, n) # optional

# 2. init
x = y
blocks = lapply(1:n, function(i) {
  list(block = x[i], nums = 1)
})

# 3. PAVA
while(T) {
  ### adjacent pooling
  block_num = 1
  new_blocks = list()
  count = 1
  while(block_num < length(blocks)) {
    i = 1
    new_blocks[count] = blocks[block_num]
    while((blocks[[block_num]]$block > blocks[[block_num + i]]$block)) {
      new_blocks[[count]]$block = c(new_blocks[[count]]$block, blocks[[block_num + i]]$block)
      new_blocks[[count]]$nums = new_blocks[[count]]$nums + blocks[[block_num + i]]$nums
      i = i + 1
      if(block_num + i > length(blocks)) break
    }  
    block_num = block_num + i
    count = count + 1
  }
  
  ### solve f(x) - update estimates
  blocks = lapply(new_blocks, function(e) {
            list(block = mean(e$block),
            nums = e$nums)
           })
  
  ### stopping rule check
  flag = 0
  for(i in 2:length(blocks)) {
    if(blocks[[i - 1]]$block > blocks[[i]]$block) flag = 1
  }
  if(flag == 0) break
}
### result
est = unlist(sapply(1:length(blocks), function(i) {
        rep(blocks[[i]]$block, blocks[[i]]$nums) 
      }))

# 4. plot
par(mfrow = c(1,1))
plot(y, main = "isotonic regression")
points(est, type = 'l', col = 'red', lwd = 2)
