# Q2

# (a)
Y = log(data[,2]);
X1 = data[,10];
X2 = data[,13];
X = cbind(matrix(1,1173,1), X1, X2);

# (b)
plot(X1, Y, main = "Logged Violent Crime Rate Vs. Per Capita Income", xlab="Per Capita Income", ylab="Logged Violent Crime Rate");
dev.copy(pdf, 'ps1_2b.pdf');
dev.off();
# Scatterplot printed on the following page.

# (c)
X_tilde = cbind(matrix(1,1173,1), X2);
beta_tilde = solve(t(X_tilde) %*% X_tilde) %*% t(X_tilde) %*% Y;
# Interpretation: beta_tilde is the estimated percentage change in the violent crime rate associated with having a shall-carry law, in the sense that it minimizes the sum of squared residuals.
# In this case, it says that states that have a shall-carry law have approximately 44% less violent crime.
# However, nothing can be inferred from this coefficient due to omitted variable bias, in particular due to per capita income and year being omitted from the regression.
# In addition, there is nothing here that establishes a causative relation. 

# (d)
beta_hat = solve(t(X) %*% X) %*% t(X) %*% Y;

# (e)
beta_hat = solve(t(X)%*%X, t(X)%*%Y);

# (f)
lmfit = lm(Y ~ X1 + X2);
print(lmfit);
# As expected, the results are approximately the same as in part (e).