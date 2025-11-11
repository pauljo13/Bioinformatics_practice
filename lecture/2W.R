# Rdata structure
## Data frame
m1 <- matrix(c(1:9),nrow = 3)
m1
d1 <- as.data.frame(m1)
d1
summary(d1)

d1 <- data.frame(c1=1:4,
                 c2=c('a','b','c','d'),
                 c3=c(T,T,F,F),
                 row.names=c('r1','r2','r3','r4'))
d1

# Quiz
myMat <- matrix(seq(1,by=2, length.out = 9), nrow = 3, byrow = T, dimnames = list(c('R1', 'R2', 'R3'), c('C1', 'C2','C3')))
myDf <- as.data.frame(myMat)
myDf$C4 <- paste('C4_R',1:3,sep = "")

myDf

## List
l1 <- list(3,c(T,T,F),matrix(1:4,nrow=2)); l1
l2 <- as.list(d1)
l2


# IF and Loop in R
## if
# Quiz
a <- sample(1:10, 1)
if (a%%2 == 1) {
  print(sprintf('a = %d ,a is an odd number', a))
} else {
  print(sprintf('a = %d ,a is an even number', a))
}

if (a%%2 == 1) b <- a+1 else b <- a; print(b)
b2 <- ifelse(a %% 2 == 1, a+1, a);b2
b3 <- if(a%%2 == 1) a+1 else a; print(b3)


# Quiz
book <-  10000
coffee <- 3000
candy <- 500
budget <- 4500
myitem <- c()
if (budget >= book) {
  myitem <- c(myitem, 'book')
} else if (budget < book & budget >= coffee){
  myitem <-  c(myitem, 'coffee')
} else if (budget < coffee & budget >= candy){
  myitem <-  c(myitem, 'candy')
} else {
  print('No mumry')
}
myitem



# While 
i <- 0
while(i < 5) {
  i = i+1
  j <- 0
  while (j < i) {
    j <- j + 1
    print(i, j)
    }
}

# Quiz
i <- 0
s <- 0
while (i < 100) {
  i <- i +1
  if(i %% 7 == 0){ 
    s <- s + i }
}; s


# Quiz
myVec <-  c('watermelon', 'banana', 'apply', 'tomato', 'peach')
vege <- c('tomato', 'carrot', 'radish')
i <- 0
while(i < length(myVec)) {
  i <- i + 1
  if(myVec[i] %in% vege) {
    next
  }
  print(myVec[i])
}

i <- 0
while(i < length(myVec)) {
  i <- i + 1
  if(myVec[i] %in% vege) {
    break
  }
  print(myVec[i])
}


# For
for (i in 1:10) {
  if (i%%2 == 1) {
    print(sprintf('i = %d ,a is an odd number', i))
  } else {
    print(sprintf('i = %d ,a is an even number', i))
  }
}

# Quiz
for (i in seq_along(myVec)) {
  x <- myVec[i]
  if (x %in% vege) {
    print(sprintf('index = %d, it is %s', i, x))
  }
}

myVec[which(myVec %in% vege)]

# repeat
multi_7 <- c()
repeat {
  i <- i +1
  if (i %% 7 == 0) {
    multi_7 <- c(multi_7, i)
  }
  if (i > 100){
    break
  }
}; multi_7


# Create a matrix using for loop
myMat <-  matrix(data = NA, nrow = 9, ncol = 9)
for (i in 1:ncol(myMat)) {
  val <- 1:9 * i
  myMat[,i] <- val
}

# Quiz
myMat2 <- matrix(c(1:81), nrow = 9, ncol = 9)
for (i in 1:ncol(myMat2)) {
  sum_val <- sum(myMat2[,i])
  print(sprintf('the sum of column %d: %s', i , sum_val))
}

# Quiz
mat <- c()
for (i in 1:3){
  mat <- rbind(mat, i * seq(2, 10 ,by = 2))
}; mat
apply(mat,1,mean)
