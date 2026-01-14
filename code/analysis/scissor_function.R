Scissor <- function(bulk_dataset, sc_dataset, phenotype, tag = NULL,
                    alpha = NULL, cutoff = 0.2, family = c("gaussian","binomial","cox"),
                    Save_file = "Scissor_inputs.RData", Load_file = NULL){
    library(Seurat)
    library(Matrix)
    library(preprocessCore)


    if (is.null(Load_file)){
        common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
        sc_exprs <- as.matrix(sc_dataset@assays$RNA@data)
        network  <- as.matrix(sc_dataset@graphs$RNA_snn)
        diag(network) <- 0
        network[which(network != 0)] <- 1
        # 合并数据集并进行分位数归一化
        dataset0 <- cbind(bulk_dataset[common,], sc_exprs[common,])         # Dataset before quantile normalization.
        dataset1 <- normalize.quantiles(dataset0)                           # Dataset after  quantile normalization.
        rownames(dataset1) <- rownames(dataset0)
        colnames(dataset1) <- colnames(dataset0)

        Expression_bulk <- dataset1[,1:ncol(bulk_dataset)]
        Expression_cell <- dataset1[,(ncol(bulk_dataset) + 1):ncol(dataset1)]
        X <- cor(Expression_bulk, Expression_cell)

        quality_check <- quantile(X)
        print("|**************************************************|")
        print("Performing quality-check for the correlations")
        print("The five-number summary of correlations:")
        print(quality_check)
        print("|**************************************************|")
        if (quality_check[3] < 0.01){
            warning("The median correlation between the single-cell and bulk samples is relatively low.")
        }
        rm(common, sc_exprs, dataset0, dataset1)
        gc()
        if (family == "binomial"){
            Y <- as.numeric(phenotype)
            z <- table(Y)
            if (length(z) != length(tag)){
                stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
            }else{
                print(sprintf("Current phenotype contains %d %s and %d %s samples.", z[1], tag[1], z[2], tag[2]))
                print("Perform logistic regression on the given phenotypes:")
            }
        }
        if (family == "gaussian"){
            Y <- as.numeric(phenotype)
            z <- table(Y)
            if (length(z) != length(tag)){
                stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
            }else{
                tmp <- paste(z, tag)
                print(paste0("Current phenotype contains ", paste(tmp[1:(length(z)-1)], collapse = ", "), ", and ", tmp[length(z)], " samples."))
                print("Perform linear regression on the given phenotypes:")
            }
        }
        if (family == "cox"){
            Y <- as.matrix(phenotype)
            if (ncol(Y) != 2){
                stop("The size of survival data is wrong. Please check Scissor inputs and selected regression type.")
            }else{
                print("Perform cox regression on the given clinical outcomes:")
            }
        }
        save(X, Y, network, Expression_bulk, Expression_cell, file = Save_file)
    }else{
        load(Load_file)
    }
    gc()
    if (is.null(alpha)){
        alpha <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    }
    for (i in 1:length(alpha)){
        set.seed(123)
        fit0 <- APML1(X, Y, family = family, penalty = "Net", alpha = alpha[i], Omega = network, nlambda = 100, nfolds = min(10,nrow(X)))
        gc()
        fit1 <- APML1(X, Y, family = family, penalty = "Net", alpha = alpha[i], Omega = network, lambda = fit0$lambda.min)
        gc()
        if (family == "binomial"){
            Coefs <- as.numeric(fit1$Beta[2:(ncol(X)+1)])
        }else{
            Coefs <- as.numeric(fit1$Beta)
        }
        Cell1 <- colnames(X)[which(Coefs > 0)]
        Cell2 <- colnames(X)[which(Coefs < 0)]
        percentage <- (length(Cell1) + length(Cell2)) / ncol(X)
        print(sprintf("alpha = %s", alpha[i]))
        print(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.", length(Cell1), length(Cell2)))
        print(sprintf("The percentage of selected cell is: %s%%", formatC(percentage*100, format = 'f', digits = 3)))

        if (percentage < cutoff){
            break
        }
        cat("\n")
        gc()
    }
    print("|**************************************************|")

    return(list(para = list(alpha = alpha[i], lambda = fit0$lambda.min, family = family),
                Coefs = Coefs,
                Scissor_pos = Cell1,
                Scissor_neg = Cell2))
}