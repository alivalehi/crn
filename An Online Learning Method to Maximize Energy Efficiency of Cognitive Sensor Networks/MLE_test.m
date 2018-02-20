lambda = .1;
n = 3000;
X = exprnd(lambda,1,n);
phat = mle(X,'distribution','exp')