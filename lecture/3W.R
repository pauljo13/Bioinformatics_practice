# User-defined function
hellowrld <-  function() {
  print('hello, world')
}
hellowrld()

echo <- function(what) {
  print(what)
}
echo('hi')

#quiz
sum_col <- function(x, y) {
  if(y > ncol(x)) {
    print(sprintf('%d is bigger than %s', y, ncol(x)))
  } else {
  return(sum(x[,y]))
  }
}
myMat2 <- matrix(1:81, nrow = 9, ncol = 9,byrow = TRUE)
sum_col(myMat2, 11)
sum_col(myMat2, 8)

#  Write & read
getwd()
setwd('/Users/bumsoojoe/Desktop/Study/Bioinformatics_practice/lecture/data/')

fileConn <- file('test.txt', open = 'wt')
writeLines(c('line1','line2'), con = fileConn)
close(fileConn)
file.show('test.txt')

# cat()
cat(c('line1','line2'), file = 'test3.txt', sep= '\n')
file.show('test3.txt')

# sink()
sink(file = 'text3.txt')
cat(c('line1','line2'), sep = '\n')
print('line1\nline2')
sink()
file.show('test3.txt')

sum_vals <- apply(myMat2, 2, sum); sum_vals
sumVec <- paste('the sum of column ',seq_along(sum_vals),': ',
                sum_vals, sep='')
sink(file='test_sink.txt')
cat(sumVec,sep="\n----------------\n")
sink()
file.show('test_sink.txt')

write.table()
myMat <- matrix(paste('v', 1:9, sep = ""), nrow = 3)
write.table(myMat, file = 'test_Mat1.txt', sep = '\t')
file.show('test_Mat1.txt')
write.table(myMat, file='test_Mat2.txt', sep=',', quote=F)
file.show('test_Mat2.txt')
write.table(myMat, file='test_Mat3.txt', sep='\t',row.names = F)
file.show('test_Mat3.txt')

# Quiz
colnames(mat) <- paste('myCol',1:ncol(mat),sep="_")
mat <- c()
for (i in c('a', 'b', 'c', 'd')) {
  mat <- cbind(mat,paste(i,1:5,sep=""))
}
write.table(mat,file='test_mat_quiz.csv',quote=F, sep=",", row.names = F)
                        

myDF <- read.table('test_mat3.txt', sep = '\t', header = T)
myDF <- read.table('test_mat2.txt', sep = ',', header = T)
myDF

fileConn <- file('test_mat4.txt', op ='wt')
writeLines(c('line1','line2','#end of line'), con = fileConn)
close(fileConn)

write.table(myMat, file='test_mat4.txt', sep='\t',row.names = F, col.names = F, append = T)
file.show('test_mat4.txt')

myDf <- read.table('test_mat4.txt')
myDf <- read.table('test_mat4.txt',header=F,sep="\t",skip=3)
myDf <- read.table('test_mat4.txt',header=F,sep="\t",skip=2,comment.char = "#")
myDf

# Quiz
myDf <- read.table('test_mat_quiz2.txt', header = T, row.names = 'rowname', sep = '\t')
#read.delim('test_mat_quiz2.txt',row.names = 1)
df_mean <- apply(myDf, myDf[,3:7], 2, mean)


# R package
.libPaths()
library(help = splines)
library(splines)
require(splines)
detach('package:splines')
search()

devtools::install_github("r-lib/devtools")
detach('package:r-lib')
