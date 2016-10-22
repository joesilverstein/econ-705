sink("C:\\Users\\Joe\\Documents\\ps6_output.txt");

###############################
cat("# 2\n\n")

data = read.table("card.dat");

###############################
cat("# 2 (a)\n")

# dependent variable
log_wage = log(data$V26);

# endogenous variables
educ = data$V4; 

# exogenous variables
exper = data$V29;
exper_2 = exper^2;
south = data$V24;
black = data$V22;

# instrumental variables
near2 = data$V2;
near4 = data$V3; 
fatheduc = data$V6;
motheduc = data$V7;

n = length(log_wage); # number of observations
l = 8; # number of exogenous variables
Y = log_wage;

# Step 1: Compute preliminary consistent estimator via 2SLS
Z = cbind(rep(1,n), near2, near4, fatheduc, motheduc, exper, exper_2, south, black);
educ_hat = Z %*% solve(t(Z) %*% Z) %*% t(Z) %*% educ;
X_hat = cbind(rep(1,n), educ_hat, exper, exper_2, south, black);
beta_2SLS = solve(t(X_hat) %*% X_hat) %*% t(X_hat) %*% Y;

# Step 2: Compute the optimal weight matrix
X = cbind(rep(1,n), educ, exper, exper_2, south, black);
e_hat_2SLS = Y - X %*% beta_2SLS; # residual from preliminary estimator (should have dimensions nx1)

# Have to use scalar multiplication here because e_hat[i] is a scalar for each i.
# Also note that Z[i,] spits out a column vector for each i even though it's actually a row vector.
# I simultaneously compute the mean across observations here so I don't have to loop again.
g_hat = matrix(0, nrow=n, ncol = l+1);
sum = matrix(0,l+1,1);
for (i in 1:n) {
  g_hat[i,] = Z[i,] * e_hat_2SLS[i];
  sum = sum + g_hat[i,];
}
g_bar = (1/n)*sum; # Numbers are insignificant because the 2SLS estimator is consistent and we have a large sample size.

g_hat_star = matrix(0,n,l+1)
for (i in 1:n) {
  g_hat_star[i,] = g_hat[i,] - g_bar;
}

sum = 0;
for (i in 1:n) {
  sum = sum + g_hat_star[i,] %*% t(g_hat_star[i,]);
}
W = solve((1/n)*sum);

# Step 3: Compute the GMM estimator
beta_GMM = solve(t(X) %*% Z %*% W %*% t(Z) %*% X) %*% t(X) %*% Z %*% W %*% t(Z) %*% Y;
cat("\nGMM Coefficients:\n")
print(t(beta_GMM))

G_hat = -(1/n) * t(Z) %*% X;
cov = (1/n) * solve(t(G_hat) %*% W %*% G_hat); # REMEMBER THE 1/n TERM!!!!!!!!!!!!!!!!!!!!
stderr = sqrt(diag(cov));
cat("\nGMM Standard Errors:\n")
print(stderr);

###############################

cat("\n# 2 (b)\n")

# Calculate the J-Statistic for overidentification:

J = n * t(g_bar) %*% W %*% g_bar;
cat("\nJ Statistic:", J)

###############################

cat("\n\n# 3 (a)\n")

l = 3;
eta = 2; # number of regressors in the model (the columns of X)
beta = 1; 
PI = eta * matrix(1,l,1); # this is what the horrible notation means
n = 500; 
R = 5000;

# Homoskedastic Case:
# Takes ~1min due to slow matrix inversion when estimating GMM
betas_2SLS = matrix(0,R,1);
betas_GMM = matrix(0,R,1);
for (i in 1:R) {
  # Generate random sample:
  Z = matrix(rnorm(n*l, 0, 1), n, l); # each row is an observation
  E = rnorm(n,0,1); # because cov_u = I_2
  V = rnorm(n,0,1);
  X = Z %*% PI + V;
  Y = X %*% beta + E; 
  
  # Step 1: compute the 2SLS estimator and save it
  beta_2SLS = solve(t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)%*%t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y;
  betas_2SLS[i] = beta_2SLS;
  
  # Step 2: Compute the optimal weight matrix
  e_hat_2SLS = Y - X %*% beta_2SLS; # residual from preliminary estimator (should have dimensions nx1)
  
  g_hat = matrix(0, n, l);
  sum = matrix(0,l,1);
  for (j in 1:n) {
    g_hat[j,] = Z[j,] * e_hat_2SLS[j];
    sum = sum + g_hat[j,];
  }
  g_bar = (1/n)*sum; # Numbers are insignificant because the 2SLS estimator is consistent and we have a large sample size.
  
  g_hat_star = matrix(0,n,l)
  for (j in 1:n) {
    g_hat_star[j,] = g_hat[j,] - g_bar;
  }
  
  sum = 0;
  for (j in 1:n) {
    sum = sum + g_hat_star[j,] %*% t(g_hat_star[j,]);
  }
  W = solve((1/n)*sum);
  
  # Step 3: Compute the GMM estimator
  betas_GMM[i] = solve(t(X) %*% Z %*% W %*% t(Z) %*% X) %*% t(X) %*% Z %*% W %*% t(Z) %*% Y;
}
cat("\n2SLS finite sample variance (Homoskedastic Errors):", var(betas_2SLS)); # 0.0001664162
cat("\nGMM finite sample variance (Homoskedastic Errors):", var(betas_GMM)); # 0.000167618

###############################

cat("\n\n# 3 (b)\n")

# Heteroskedastic Case:
# Takes ~2min
betas_2SLS = matrix(0,R,1);
betas_GMM = matrix(0,R,1);
for (i in 1:R) {
  # Generate random sample with heteroskedastic errors as specified in the instructions:
  Z = matrix(0,n,l);
  E = matrix(0,n,1);
  V = matrix(0,n,1);
  for (j in 1:n) {
    z = matrix(rnorm(l, 0, 1), 1, l); # each row is an observation
    Z[j,] = z; # save as a new observation in the big matrix
    E[j] = rnorm(1,0,abs(z[1])); # Takes sd as argument -- NOT variance!!! 
    V[j] = rnorm(1,0,abs(z[2]));
  }
  X = Z %*% PI + V;
  Y = X %*% beta + E; 
  
  # Step 1: compute the 2SLS estimator and save it
  beta_2SLS = solve(t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)%*%t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y;
  betas_2SLS[i] = beta_2SLS;
  
  # Step 2: Compute the optimal weight matrix
  e_hat_2SLS = Y - X %*% beta_2SLS; # residual from preliminary estimator (should have dimensions nx1)
  
  g_hat = matrix(0, n, l);
  sum = matrix(0,l,1);
  for (j in 1:n) {
    g_hat[j,] = Z[j,] * e_hat_2SLS[j];
    sum = sum + g_hat[j,];
  }
  g_bar = (1/n)*sum; # Numbers are insignificant because the 2SLS estimator is consistent and we have a large sample size.
  
  g_hat_star = matrix(0,n,l)
  for (j in 1:n) {
    g_hat_star[j,] = g_hat[j,] - g_bar;
  }
  
  sum = 0;
  for (j in 1:n) {
    sum = sum + g_hat_star[j,] %*% t(g_hat_star[j,]);
  }
  W = solve((1/n)*sum);
  
  # Step 3: Compute the GMM estimator
  betas_GMM[i] = solve(t(X) %*% Z %*% W %*% t(Z) %*% X) %*% t(X) %*% Z %*% W %*% t(Z) %*% Y;
}
cat("\n2SLS finite sample variance (Heteroskedastic Errors, n=500, l=3):", var(betas_2SLS)); # 0.0002728276
cat("\nGMM finite sample variance (Heteroskedastic Errors, n=500, l=3):", var(betas_GMM)); # 0.0002117497

###############################

cat("\n\n# 3 (c)\n")

n = 1000; # this is the only difference from 3 (b)

# Takes ~4min
betas_2SLS = matrix(0,R,1);
betas_GMM = matrix(0,R,1);
for (i in 1:R) {
  # Generate random sample with heteroskedastic errors as specified in the instructions:
  Z = matrix(0,n,l);
  E = matrix(0,n,1);
  V = matrix(0,n,1);
  for (j in 1:n) {
    z = matrix(rnorm(l, 0, 1), 1, l); # each row is an observation
    Z[j,] = z; # save as a new observation in the big matrix
    E[j] = rnorm(1,0,abs(z[1])); # Takes sd as argument -- NOT variance!!! 
    V[j] = rnorm(1,0,abs(z[2]));
  }
  X = Z %*% PI + V;
  Y = X %*% beta + E; 
  
  # Step 1: compute the 2SLS estimator and save it
  beta_2SLS = solve(t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)%*%t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y;
  betas_2SLS[i] = beta_2SLS;
  
  # Step 2: Compute the optimal weight matrix
  e_hat_2SLS = Y - X %*% beta_2SLS; # residual from preliminary estimator (should have dimensions nx1)
  
  g_hat = matrix(0, n, l);
  sum = matrix(0,l,1);
  for (j in 1:n) {
    g_hat[j,] = Z[j,] * e_hat_2SLS[j];
    sum = sum + g_hat[j,];
  }
  g_bar = (1/n)*sum; # Numbers are insignificant because the 2SLS estimator is consistent and we have a large sample size.
  
  g_hat_star = matrix(0,n,l)
  for (j in 1:n) {
    g_hat_star[j,] = g_hat[j,] - g_bar;
  }
  
  sum = 0;
  for (j in 1:n) {
    sum = sum + g_hat_star[j,] %*% t(g_hat_star[j,]);
  }
  W = solve((1/n)*sum);
  
  # Step 3: Compute the GMM estimator
  betas_GMM[i] = solve(t(X) %*% Z %*% W %*% t(Z) %*% X) %*% t(X) %*% Z %*% W %*% t(Z) %*% Y;
}
cat("\n2SLS finite sample variance (Heteroskedastic Errors, n=1000, l=3):", var(betas_2SLS)); # 0.0001376627
cat("\nGMM finite sample variance (Heteroskedastic Errors, n=1000, l=3):", var(betas_GMM)); # 0.0001047472

###############################

cat("\n\n# 3 (d)\n")

n = 500; # set back to value from 3 (b)
l = 5; # the only difference from 3 (b)
PI = eta * matrix(1,l,1); # have to make PI longer

# Takes ~2min
time = proc.time();
betas_2SLS = matrix(0,R,1);
betas_GMM = matrix(0,R,1);
for (i in 1:R) {
  # Generate random sample with heteroskedastic errors as specified in the instructions:
  Z = matrix(0,n,l);
  E = matrix(0,n,1);
  V = matrix(0,n,1);
  for (j in 1:n) {
    z = matrix(rnorm(l, 0, 1), 1, l); # each row is an observation
    Z[j,] = z; # save as a new observation in the big matrix
    E[j] = rnorm(1,0,abs(z[1])); # Takes sd as argument -- NOT variance!!! 
    V[j] = rnorm(1,0,abs(z[2]));
  }
  X = Z %*% PI + V;
  Y = X %*% beta + E; 
  
  # Step 1: compute the 2SLS estimator and save it
  beta_2SLS = solve(t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)%*%t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y;
  betas_2SLS[i] = beta_2SLS;
  
  # Step 2: Compute the optimal weight matrix
  e_hat_2SLS = Y - X %*% beta_2SLS; # residual from preliminary estimator (should have dimensions nx1)
  
  g_hat = matrix(0, n, l);
  sum = matrix(0,l,1);
  for (j in 1:n) {
    g_hat[j,] = Z[j,] * e_hat_2SLS[j];
    sum = sum + g_hat[j,];
  }
  g_bar = (1/n)*sum; # Numbers are insignificant because the 2SLS estimator is consistent and we have a large sample size.
  
  g_hat_star = matrix(0,n,l)
  for (j in 1:n) {
    g_hat_star[j,] = g_hat[j,] - g_bar;
  }
  
  sum = 0;
  for (j in 1:n) {
    sum = sum + g_hat_star[j,] %*% t(g_hat_star[j,]);
  }
  W = solve((1/n)*sum);
  
  # Step 3: Compute the GMM estimator
  betas_GMM[i] = solve(t(X) %*% Z %*% W %*% t(Z) %*% X) %*% t(X) %*% Z %*% W %*% t(Z) %*% Y;
}
time = proc.time() - time;
cat("\n2SLS finite sample variance (Heteroskedastic Errors, n=500, l=5):", var(betas_2SLS)); # 0.0001409414
cat("\nGMM finite sample variance (Heteroskedastic Errors, n=500, l=5):", var(betas_GMM)); # 0.000113114

###############################

sink();