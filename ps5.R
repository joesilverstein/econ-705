sink("C:\\Users\\Joe\\Documents\\ps5_output.txt");

cat("# 3\n\n")

data = read.table("card.dat")

cat("# 3 (a)\n")

log_wage = log(data$V26);
educ = data$V4;
exper = data$V29;
exper_2 = exper^2;
south = data$V24;
black = data$V22;
olsspec_a = lm(log_wage ~ educ + exper + exper_2 + south + black);
summary(olsspec_a)

cat("# 3 (b)\n\n")

n = length(log_wage);

near4 = data$V3;
Z = cbind(rep(1,n), exper, exper_2, south, black, near4);
educ_hat = Z %*% solve(t(Z) %*% Z) %*% t(Z) %*% educ;
olsspec_b = lm(log_wage ~ educ_hat + exper + exper_2 + south + black);
out = summary(olsspec_b);

beta_2SLS = out$coefficients[,1];
X = cbind(rep(1,n), educ, exper, exper_2, south, black);
e_hat = log_wage - X %*% beta_2SLS; # residual from fitting model using 2SLS estimates (NOT the residual from the 2nd stage!)

cat("Coefficient Estimates:\n")
cat("Intercept:", out$coefficients[1,1], "\n")
cat("educ:", out$coefficients[2,1], "\n")
cat("exper:", out$coefficients[3,1], "\n")
cat("exper_2:", out$coefficients[4,1], "\n")
cat("south:", out$coefficients[5,1], "\n")
cat("black:", out$coefficients[6,1], "\n\n")

# Estimation of the covariance matrix:

sum = 0; # initialization
for (i in 1:n) {
  sum = sum + Z[i,] %*% t(X[i,]); # always include column vector of ones for both Z and X
}
Q_hat = (1/n) * sum;

sum = 0; 
for (i in 1:n) {
  sum = sum + Z[i,] %*% t(Z[i,]);
}
W_hat = (1/n) * sum;

sum = 0; 
for (i in 1:n) {
  sum = sum + Z[i,] %*% t(Z[i,]) * (e_hat[i]^2); 
}
omega_hat = (1/n) * sum;

V_hat = solve(t(Q_hat) %*% W_hat %*% Q_hat) %*% t(Q_hat) %*% W_hat %*% omega_hat %*% W_hat %*% Q_hat %*% solve(t(Q_hat) %*% W_hat %*% Q_hat);

cat("Asymptotic Standard Error Estimates:\n")
cat("Intercept:", sqrt(V_hat[1,1]), "\n")
cat("educ:", sqrt(V_hat[2,2]), "\n")
cat("exper:", sqrt(V_hat[3,3]), "\n")
cat("exper_2:", sqrt(V_hat[4,4]), "\n")
cat("south:", sqrt(V_hat[5,5]), "\n")
cat("black:", sqrt(V_hat[6,6]), "\n\n")

cat("# 3 (c)\n\n")

near2 = data$V2;
fatheduc = data$V6;
motheduc = data$V7;
Z = cbind(rep(1,n), exper, exper_2, south, black, near2, near4, fatheduc, motheduc);
educ_hat = Z %*% solve(t(Z) %*% Z) %*% t(Z) %*% educ;
olsspec_c = lm(log_wage ~ educ_hat + exper + exper_2 + south + black);
out = summary(olsspec_c);

beta_2SLS = out$coefficients[,1];
X = cbind(rep(1,n), educ, exper, exper_2, south, black);
e_hat = log_wage - X %*% beta_2SLS; # residual from fitting model using 2SLS estimates (NOT the residual from the 2nd stage!)


cat("Coefficient Estimates:\n")
cat("Intercept:", out$coefficients[1,1], "\n")
cat("educ:", out$coefficients[2,1], "\n")
cat("exper:", out$coefficients[3,1], "\n")
cat("exper_2:", out$coefficients[4,1], "\n")
cat("south:", out$coefficients[5,1], "\n")
cat("black:", out$coefficients[6,1], "\n\n")

# Estimation of the covariance matrix:

sum = 0;
for (i in 1:n) {
  sum = sum + Z[i,] %*% t(X[i,]);
}
Q_hat = (1/n) * sum;

sum = 0; 
for (i in 1:n) {
  sum = sum + Z[i,] %*% t(Z[i,]);
}
W_hat = (1/n) * sum;

sum = 0;
for (i in 1:n) {
  sum = sum + Z[i,] %*% t(Z[i,]) * (e_hat[i]^2); 
}
omega_hat = (1/n) * sum;

V_hat = solve(t(Q_hat) %*% W_hat %*% Q_hat) %*% t(Q_hat) %*% W_hat %*% omega_hat %*% W_hat %*% Q_hat %*% solve(t(Q_hat) %*% W_hat %*% Q_hat);

cat("Asymptotic Standard Error Estimates:\n")
cat("Intercept:", sqrt(V_hat[1,1]), "\n")
cat("educ:", sqrt(V_hat[2,2]), "\n")
cat("exper:", sqrt(V_hat[3,3]), "\n")
cat("exper_2:", sqrt(V_hat[4,4]), "\n")
cat("south:", sqrt(V_hat[5,5]), "\n")
cat("black:", sqrt(V_hat[6,6]), "\n\n")

sink();