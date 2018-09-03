
############## KM plot
kmPlot <- function(gene, project, sep='median', type='OS') {
    
    prj <- strsplit(project, split = '-', fixed = TRUE)[[1]][2]
    
    rnaExpr <- load(paste(project, paste('rnaExpr',prj, 'rda',sep='.'), sep='/'))
    rnaExpr <- get(rnaExpr)
    
    survDa <- getTCGASurvival(project)
    survDa
    
    samples <- intersect(colnames(rnaExpr), rownames(survDa))
    #expr=rnaExpr[symbol2ensemblFun(gene),samples]
    
    exprDa <- rnaExpr[gene,samples]
    survDa <- survDa[samples,]
    
    if (sep=='1stQu') {
        thresh <- as.numeric(summary(exprDa)[2])
    } else if (sep=='median') {
        thresh <- as.numeric(summary(exprDa)[3])
    } else if (sep=='mean') {
        thresh <- as.numeric(summary(exprDa)[4])
    } else if (sep=='3rdQu') {
        thresh <- as.numeric(summary(exprDa)[5])
    }
    
    exprGroup <- exprDa > thresh
    
    nH <- sum(exprGroup)
    nL <- sum(!exprGroup)
    
    
    if (type == 'OS') {
        survDa <- data.frame(daysToDeath=survDa$OS,vitalStatus=survDa$OS_status, exprGroup)
    } else if (type == 'RFS') {
        survDa <- data.frame(daysToDeath=survDa$RFS,vitalStatus=survDa$RFS_status, exprGroup)
    }
    
    
    sdf <- survdiff(Surv(survDa$daysToDeath, survDa$vitalStatus) ~ exprGroup)
    pValue <- format(pchisq(sdf$chisq, length(sdf$n)-1, 
                            lower.tail = FALSE),digits=3)
    #pValue <- format(1-pchisq(sdf$chisq, df=1),digits=3)
    
    HR = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
    upper95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
    lower95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
    
    HR <- format(HR, digits = 3)
    upper95 <- format(upper95, digits = 3)
    lower95 <- format(lower95, digits = 3)
    
    label1 <- paste('HR = ', HR, '(', lower95, '-', upper95, ')', sep='')
    label2 <- paste('logrank P value = ', pValue, sep='')
    
    fit <- survfit(Surv(daysToDeath, vitalStatus) ~ exprGroup, data=survDa)
    
    lgdXpos <- 1/1.3
    lgdYpos = 0.9
    
    xpos = max(survDa$daysToDeath, na.rm=TRUE)/1.8
    ypos1 = 0.8
    
    ggsurvplot(fit, data=survDa, pval = paste(label1, '\n', label2), pval.coord = c(xpos, ypos1),
               pval.size=5,
               font.main = c(16, 'bold', 'black'), conf.int = FALSE, 
               title = project,
               legend = c(lgdXpos, lgdYpos), 
               #color = c('blue', 'green'),
               palette= c('blue', 'red'),
               legend.labs = c(paste('lowExp (N=',nL,')',sep=''), 
                               paste('highExp  (N=',nH,')',sep='')),  
               legend.title='group',
               xlab = paste(type,'(months)'), ylab = 'Survival probability',
               font.x = c(16), font.y = c(16), ylim=c(0,1), #16
               ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           #panel.border = element_rect(colour='black'),
                                           panel.border = element_blank(),
                                           panel.background = element_blank(),
                                           legend.text = element_text(size=12),#14
                                           legend.title = element_text(size=14),
                                           axis.text = element_text(size=14, color='black')))

}




############## KM analysis
kmTest <- function(exprDa, daysToDeath, vitalStatus, sep='median') {
    
    if(is.vector(exprDa)) {
        DEG <- exprDa
        
        if (sep=='1stQu') {
            thresh <- as.numeric(summary(DEG)[2])
        } else if (sep=='median') {
            thresh <- as.numeric(summary(DEG)[3])
        } else if (sep=='mean') {
            thresh <- as.numeric(summary(DEG)[4])
        } else if (sep=='3rdQu') {
            thresh <- as.numeric(summary(DEG)[5])
        }
        
        exprGroup <- DEG > thresh
        
        sdf <- survdiff(Surv(daysToDeath, vitalStatus) ~ exprGroup)
        pValue <- format(pchisq(sdf$chisq, length(sdf$n)-1, 
                                lower.tail = FALSE),digits=3)
        #p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
        
        HR = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
        upper95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
        lower95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
        
        kmDEGs <- c(HR, lower95, upper95, pValue)
        
    } else {
        kmDEGs <- c()
        for (i in seq_len(nrow(exprDa))) {
            DEG <- unlist(exprDa[i,])
            
            if (sep=='1stQu') {
                thresh <- as.numeric(summary(DEG)[2])
            } else if (sep=='median') {
                thresh <- as.numeric(summary(DEG)[3])
            } else if (sep=='mean') {
                thresh <- as.numeric(summary(DEG)[4])
            } else if (sep=='3rdQu') {
                thresh <- as.numeric(summary(DEG)[5])
            }
            
            exprGroup <- DEG > thresh
            
            sdf <- survdiff(Surv(daysToDeath, vitalStatus) ~ exprGroup)
            pValue <- format(pchisq(sdf$chisq, length(sdf$n)-1, 
                                    lower.tail = FALSE),digits=3)
            #p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
            
            HR = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
            upper95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
            lower95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
            
            kmDEGs <- rbind(kmDEGs, c(HR, lower95, upper95, pValue))
            
        }
        
        rownames(kmDEGs) <- rownames(exprDa)
        colnames(kmDEGs) <- c('HR','lower95','upper95','pValue')
        kmDEGs <- data.frame(symbol=ensembl2symbolFun(rownames(exprDa)), kmDEGs)
        #kmDEGs$FDR <- p.adjust(kmDEGs$pValue, method='fdr')
        
        #o <- order(coxphDEGs$pValue)
        #coxphDEGs <- coxphDEGs[o,]
        
    }
    return (kmDEGs)
}




################ ID conversion
ensembl2symbolFun <- function(ensemblID, info='symbol') {
    geneInfo <- biotype[match(ensemblID, biotype$ensemblID),]
    geneSymbol <- geneInfo$geneSymbol
    
    if (info=='symbol') {
        return (geneSymbol)
    } else if (info=='all') {
        return (geneInfo)
    }
}


symbol2ensemblFun <- function(symbol, info='ensemblID') {
    geneInfo <- biotype[match(symbol, biotype$geneSymbol),]
    ensemblID <- geneInfo$ensemblID
    
    if (info=='ensemblID') {
        return (ensemblID)
    } else if (info=='all') {
        return (geneInfo)
    }
}


################# Get survival data
getTCGASurvival <- function(project) {
    prj <- strsplit(project, split = '-', fixed = TRUE)[[1]][2]
    #db = dbConnect(SQLite(), dbname="TCGA.survival.sqlite")
    db = dbConnect(SQLite(), dbname="TCGA.survival.sqlite")
    dbTable <- paste('survival',prj,sep='.')
    res <- dbReadTable(db, dbTable, row.names=TRUE)
    return (res)
    dbDisconnect(db)
}

