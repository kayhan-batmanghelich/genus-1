require("fuzzyforest")
require("pROC")

y <- read.csv("../datasets_for_r/response.txt", header=F)
X <- read.csv("../datasets_for_r/XB.csv", sep=",")

# set up 
xrow <- dim(X)[1]
xcol <- dim(X)[2]
test_size <- ceiling(xrow*.1)
train_size <- xrow - test_size
mod_mem <- c(rep("brain", xcol))

# saving the results
auc_save <- c()
feat_save <- list()

for (iter in 1:400) {
    # assign training and testing data
    train_sample_idx <- sample(nrow(X), train_size)
    test_sample_idx <- setdiff(seq(xrow), train_sample_idx)
    X_train <- X[train_sample_idx, ]
    y_train <- y[train_sample_idx, ]
    X_test <- X[test_sample_idx, ]
    y_test <- y[test_sample_idx, ]
   
    # fit the model
    ff_fit <- ff(X_train, as.factor(as.numeric(unlist(y_train))), module_membership=mod_mem)
    # predict on test data
    pred <- predict(ff_fit, new_data=X_test)
    # compute auc of roc
    auc_score <- auc(y_test, as.numeric(pred))
    auc_save[iter] <- auc_score
    feat_save[[iter]] <- ff_fit$feature_list
}


