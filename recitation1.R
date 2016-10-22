# R Basics -> Recitation 1 (09/08/2014)

#----------------------------------------------------------------
# 1. Data types
#----------------------------------------------------------------
# Define vector
x = c(1,2,3,4);
y = c(3,4,5,6);
z = seq(1,9,1);  # additive sequence:  seq(start, end, stepsize)
w = (1:21);      # equivalent to seq(1,21,1)

# Define matrix
xmat = matrix(x, 2, 2); # form 2 by 2 matrix based on vector x
ymat = matrix(y, 2, 2); # same as above
zmat = matrix(z, 3, 3); # form 3 by 3 matrix based on vector z
omat = matrix(0, 3, 3); # form 3 by 3 matrix with all elements zeros
dmat = diag(1,3,3);     # form 3 by 3 diagonal matrix with diagonal elements one


# See what you have defined
ls();

# clear command window (console) : ctrl+L
# - this does not mean that you cleared variable from the memory

# clear (remove) variable from the memory
ls();  # let's check what we have
rm(w); # remove w from the memory
ls();  # there is no variable w

#----------------------------------------------------------------
# 2. Various operations
#----------------------------------------------------------------
# a) Indexing
print(x);      # recall x was a vector
x[1];          # pick first element of x
x[c(2,4)]; # pick second and fourth element of x
x[1] = 1000;   # change value in the first element of x
print(x);      # check the first element

# you can do samething with matrix and etc
print(xmat);  # recall x was a matrix
xmat[2,2];    # (2,2)th element of xmat
xmat[1:2, 2]; # take element row 1~2, column 2
xmat[1, ];    # take element row 1  , column all

# b) combining 
cbind(x,y); # column-wise binding
rbind(x,y); # row-wise binding

# c) matrix operation
print(xmat); # recall xmat
print(ymat); # recall ymat

t(xmat);     # Transpose

# Multiplication: two different multiplication
xmat*ymat;   # element-wise multiplication
xmat%*%ymat; # matrix multiplication

# Inverse of matrix
xinv = solve(xmat);   # inverse of matrix x
xmat%*%xinv;          # check xinv is whether inv(xmat)

#----------------------------------------------------------------
# Load Data
#----------------------------------------------------------------
# a) Manage working directory
getwd();             # see current working directory

# b) Load data
#mydata = read.table("dataset.txt");
#head(mydata); # show first 6 rows
#tail(mydata); # show last 6 rows

# c) Load CSV data
psdata = read.csv(file="caschool.csv",head=TRUE,sep=",")
head(psdata)

#----------------------------------------------------------------
# Looping
#----------------------------------------------------------------
# Looping (For) : for (i = (1:n)) {statement}
n        = 10;
cumulsum = matrix(1, n, 1); # create dummy matrix
for (i in (2:n)){
  print(i);
  cumulsum[i,1] = cumulsum[i-1,1] + i;
}
print(cumulsum);


# Looping (While) : while ( conditions ) {statement}
score705 = 10;
while (score705 < 90) {
  score705 = rnorm(1,70,10); #rnorm(n, m, s): generate n draws from N(b,c^2) 
  print(score705);
}

#----------------------------------------------------------------
# Create your own function
#----------------------------------------------------------------
rm(list= ls()); # let's delete all variables in the memory
ls(); # check again

# you can define a function (it is also object) in the script
# once it is defined, you can use it over and over

# function to calculate beta coefficient of the linear projection of Y onto space spanned by X
# systax: function name = function(arg1, arg2, ...) { specify what you want to do in this function}
myols = function(y,x) {
  myols = solve(t(x)%*%x) %*% (t(x) %*% y);
};
ls(); # see now our function is on the memory (that means that we can use it.

# Let's try to use this function
n = 100;
Y = rnorm(n, 0, 10);   # draw 100 random number from N(0,100)
X = rnorm(2*n, 0, 10); # draw 200 random number from N(0,100)
dim(X) = cbind(n,2);   # transform 200 by 1 vector into 100 by 2 matrix
# you could have done in other ways. For example, try X = matrix(c(X), 100, 2);
summary(X);            # print out summary statistics 

beta = myols(Y,X);
print(beta);

# Function with multiple outputs
myols2 = function(y,x){
  beta     = solve(t(x)%*%x) %*% (t(x) %*% y);
  residual = y-x%*%beta;
  output = list(b=beta, r=residual);
  return(output);
};
output = myols2(Y,X);
print(output);

#----------------------------------------------------------------
# scatter plot and save
#----------------------------------------------------------------
z = rnorm(n, 0, 10);
x = rnorm(n,0,10);
plot(z,x, main = "Scatter plot example", xlab = "Z", ylab = "X") # Draw plot
dev.copy(pdf, 'example1.pdf'); # save plot in pdf format
dev.off(); # close device for the previous figure.

#----------------------------------------------------------------
# 2013 PS1-Q5
#----------------------------------------------------------------

#function that computes E[fn(X)]

EfnX=function(n){
  cumsum=0;
  for(i in (1:4^{n})){
    cumsum=cumsum+pnorm(-sqrt(i/2^{n}));
  }
  output=1/2^{n-1}*cumsum;
};

# whether function works fine with n=1
output=EfnX(1);
compare=pnorm(-sqrt(1/2))+pnorm(-sqrt(2/2))+pnorm(-sqrt(3/2))+pnorm(-sqrt(4/2));

# Tabulate the sequence E[fn(X)] for n=1,...,8

k=8;
Eseq=matrix(0,k,1);
for(i in (1:k)){
  print(i);
  Eseq[i,1]=Eseq[i,1]+EfnX(i);
}
print(Eseq);


