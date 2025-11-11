# R -> RStudio
print('hello, world')

# Quiz
a <- log2(128)
b <- log2(256)

a / b

'I love "R"'
"I love 'R'"

# paste()
a <- 'I'
b <- 'Love'
c <- 'R'

paste(a, b, c, sep = '')
paste(a, b, c, sep = '=')

paste(1:10,'a')
paste(1:10,'a', sep = '-')
paste(1:10,'a', sep = '-', collapse = '/')

# Quiz
idx.string <- paste('Index',1:5,sep = ':', collapse = ', ')
idx.string

# nchar()
a <- 'I love R'
length(a)
nchar(a)

# substring(a, start, stop) / substr(a, start, stop)
a[1]
substring(a,1,1)
substring(a,3,6)
substr(a,3,6)

#sprintf()
sprintf('I need to study Python%d' , 2)
sprintf('I need to study Python%s' , 2.7)
sprintf('I need to study Python%.2f' , 2.7)

# Quiz
circ <- 'pi'
a <- sprintf('The ratio of a circle\'s circumference to its diameter is called %s, and the value of %s is %.2f', circ,circ, pi)
a

# Vector
a <- c(1:3)
b <- c('abc', 1:3)
d <- c(a, c('1','2','3'))

1:-10

# seq()
seq(2,10)
seq(2,10, by=2)
seq(2,10, length.out = 3)
seq(2,10, along.with = c(1,2,2,3))

# rep()
rep(1:4)
rep(1:4, times = 2)
rep(1:4, length.out = 5)
rep(1:4, each = 2)
rep(1:4, each = 2, times = 3)
rep(1:4, c(3,1,2,1))

v1 <- c(1:4)
v2 <- c(1:2)
v1 + v2
v1 * v2
v1 - v2
v1 / v2

a = c(FALSE, TRUE)
!a


a <- paste('vector',1:5, sep='')
a
a[1]
a[c(1,3)]
a[1:3]
a[6]

a <- -3:3
a > 0
a[a>0]
a[!a>0]
a[a==1]
a[a %in% seq(-3,1)]

which(a > 0)
a[which(a > 0)]
a[which(!a > 0)]
which(!(a %in% c(1,0)))

names(a)

b <- paste('vector',seq_along(a), sep='')
names(a) <- b
a
names(a)
a['vector1']
a[c('vector1','vector3')]
a[7] <- 7;a
a[a < 0] <- 0;a
a['vector3'] = '3';a

#Quiz
myVec <- paste('vec',rep(1:3,each=2))
names(myVec)[2] <- 'rep.of.vec1'
myVec[5] <- '0'
myVec <- myVec[-1]
myVec

# factor
v1 <- c('C','C','T','C','A','G','C')
is.vector(v1)
is.factor(v1)
f1 <- as.factor(v1)
f1
is.factor(f1)
is.vector(f1)

f2 <- factor(f1, levels = c('A', 'T', 'C', 'G'))
f2
f2[1] <- 'N'
f2

#Quiz
myVec <- c('wetermelon', 'banana','apple','tomato','peach')
myFac <- factor(myVec, levels = c('wetermelon', 'banana','apple','tomato','peach'))
myFac


# matrix
matrix(c(1:8), nrow = 2)
matrix(c(1:8), nrow = 3, byrow = T)
matrix(c(1:8), nrow = 2, byrow = T)

# cbind, rbind
a1 <- 1:4
a2 <- 5:8
m1 <- cbind(a1, a2)
m1

m2 <- rbind(a1, a2);m2
dim(m1)
dim(m2)

matrix(c(1:8), nrow = 2, byrow = T, dimnames = list(c(paste0('r',1:2)), c(paste0('c', 1:4))))
rownames(m1) <- c(paste0('r',1:4))
colnames(m1) <- c('c1','c2')
m1
rownames(m1)[2] <- 'row2'
m1[4,1]
m1['row2',]

which(m1>3)
which(m1>3, arr.ind = T)

m1[m1 > 3] <-  10
m1

m1 <- matrix(1:9, nrow=3)
m2 <- matrix(1:9, nrow = 3, byrow = T)

#Quiz
myMat <- matrix(seq(1,by=2,length.out =9), nrow=3, byrow = T)
rownames(myMat) <- c(paste0('r',1:3))
colnames(myMat) <- c(paste0('c',1:3))
myMat
