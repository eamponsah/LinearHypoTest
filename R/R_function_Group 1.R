######### linear combination - test ##############


#### t-test
#' Hypothesis testing for one linear combination
#'
#' Used to conduct an hypothesis testing of a single linear combination of parameters in a regression model.
#' @param model fitted model;
#' @param coef.test vector specifying linear combination to test;
#' @param hypo_values hypothesized or null value;
#' @return Returns results of hypothesis testing for a single linear combination
#' @examples
#' lm.fit <- lm(Sepal.Length ~Sepal.Width+Petal.Width, data=iris)
#' t.test_coef(lm.fit,c(0,2,4),4)

#' @export
t.test_coef <- function(model, coef_test, hypo_values){
B <- as.matrix(model$coefficients) #retrieve beta values
Sigma <- vcov(model) #var-cov matrix
A=as.matrix(coef_test) #transformation matrix

df <- model$df.residual
estimate <- t(A) %*% B
st.err <- as.numeric(sqrt(t(A) %*% Sigma %*% A))

test.stat <- as.numeric((estimate - hypo_values)/st.err)
p_value <- pt(q=abs(test.stat), df=df, lower.tail = F, log.p = FALSE)

cat(paste0("     Test for Linear Combination", "\n"))
cat(paste0(model$call, "\n"))
cat(paste0(" ", "\n"))
cat(paste0("Results: ", "\n"))
cat(paste0("t-test with ", df, " degrees of freedom", "\n"))
cat(paste0("Alternative: ", "the linear combination equals the specified value", "\n"))
cat(paste0(" ", "\n"))

result <- cbind(estimate, st.err, test.stat, p_value)
colnames(result) <- c("Estimate", "Std. error", "t value", "P(T>|t|)")

return(result)

}


#### f-test
#' Hypothesis testing for several linear combinations
#'
#' Used to conduct an hypothesis testing of a one or more linear combinations of parameters in a regression model.
#' @param model fitted model;
#' @param coef.test matrix specifying linear combinations to test;
#' @param hypo_values vector specifying hypothesized or null values;
#' @return Returns results of hypothesis testing for multiple linear combinations
#' @examples
#' lm.fit <- lm(Sepal.Length ~Sepal.Width+Petal.Width, data=iris)
#' my_matrix <- matrix(c(0,1,-1,-1,1,0), byrow = TRUE, nrow = 2)
#' my_hyp <- c(-0.6,-3)
#' f.test_coef(lm.fit, my_matrix, my_hyp)
#' @export
f.test_coef <- function(model, # fitted model from user
                         coef_test, # linear combination matrix to be applied to estimate
                         hypo_values # "right hand side" of hypothesis (row vector)
                        ) {
  B <- as.matrix(model$coefficients) # model coefficient, including intercept!
  Sigma <- vcov(model) # variance-covariance matrix for estimates
  A <- coef_test # transformation matrix (obtained from the supplied linear combination matrix)

  df <- model$df.residual # degrees of freedom
  q.num <- nrow(A) # number of linear combinations to be tested
  estimate <- A %*% B # computed estimate for transformation
  var.est <- A %*% Sigma %*% t(A) # variance of linear transformation

  # Compute test statistic and corresponding p-value
  test.stat <- (as.numeric(t(estimate - as.matrix(hypo_values)) %*%
                            solve(var.est) %*% (estimate - as.matrix(hypo_values)))) /
    (q.num * sum(model$residuals^2)/df)
  p_value <- pf(q=abs(test.stat), df1=q.num, df2=df, lower.tail = F, log.p = FALSE)

  cat(paste0("     Test for Linear Combinations", "\n"))
  cat(paste0(model$call, "\n"))
  cat(paste0(" ", "\n"))
  cat(paste0("Results: ", "\n"))
  cat(paste0("F-test with ", q.num, " and ", df, " degrees of freedom", "\n"))
  cat(paste0("Alternative: ", "the linear combination equals the specified vector", "\n"))
  cat(paste0(" ", "\n"))

  result <- cbind(test.stat, p_value)
  colnames(result) <- c("F value", "P(F>f)")

  return(result)
}



# two-sided confidence interval
#' Confidence Interval for one linear combination
#'
#' Used to construct confidence interval interval for a single linear combination of parameters in a regression model.
#' @param model fitted model;
#' @param coef.test vector specifying linear combination to test;
#' @param conf.level confidence level;
#' @return Returns results of confidence interval for a single linear combination
#' @examples
#' lm.fit <- lm(Sepal.Length ~Sepal.Width+Petal.Width, data=iris)
#' coef_conf.int(lm.fit,c(0,2,4),0.99)
coef_conf.int <- function(model, coef_test, conf.level=0.95){
  B <- as.matrix(model$coefficients)
  Sigma <- vcov(model) #varcov matrix
  A=as.matrix(coef_test) #transformation matrix

  st.err <- as.numeric(sqrt(t(A) %*% Sigma %*% A))
  df <- model$df.residual

  tvalue <- qt(p=(1-conf.level)/2, df=df,lower.tail = FALSE, log.p = FALSE)
  margin_err <- tvalue * st.err
  estimate <- (t(A) %*% B)
  lwr_limit <- estimate - margin_err
  upp_limit <- estimate + margin_err
  conf_interval <- c(estimate, lwr_limit, upp_limit)
  names(conf_interval) <- c("Estimate", "lower limit", "upper limit")

  cat(paste0(conf.level*100, "%", " Confidence Interval for linear combination", "\n"))
  cat(paste0(model$call, "\n"))
  cat(paste0(" ", "\n"))
  return(conf_interval)
}

