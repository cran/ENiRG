enirg.predict <-
function(enirg.results, qtegv.maps = NULL, qlegv.maps = NULL, load.map = FALSE, 
    method = "normal") {
    if (class(enirg.results) != "enirg") 
        stop("The function predict.enirg needs an object of class 'enirg'!")
    if (enirg.results$nf > 4) 
        stop("The function predict.enirg does not allow nf > 4!")
    egv.maps <- c(qlegv.maps, qtegv.maps)
    HSM <- list()
    if (!is.null(egv.maps)) {
        cat("\n\nCalculating li prediction maps ...\n\n")
        if (method == "normal") {
            for (j in 1:(enirg.results$nf + 1)) {
                execGRASS("r.mapcalculator", outfile = paste("predli_", colnames(enirg.results$co)[j], 
                  sep = ""), formula = "0", flags = "overwrite")
                for (i in 1:length(egv.maps)) {
                  execGRASS("r.mapcalculator", amap = egv.maps[i], outfile = "temp", 
                    formula = paste("A*", enirg.results$co[i, j], sep = ""), flags = c("overwrite"))
                  execGRASS("r.mapcalculator", amap = paste("predli_", colnames(enirg.results$co)[j], 
                    sep = ""), bmap = "temp", outfile = paste("predli_", colnames(enirg.results$co)[j], 
                    sep = ""), formula = "A+B", flags = "overwrite")
                }
            }
            cat("\n\nCleaning GRASS mapset!\n\n")
            execGRASS("g.remove", rast = "temp")
        }
        if (method == "large") {
            for (j in 1:(enirg.results$nf + 1)) {
                cc <- paste("(", enirg.results$co[, j], ")", sep = "")
                pre1 <- paste(egv.maps, cc, sep = " * ")
                calc <- paste("r.mapcalc '", paste("predli_", colnames(enirg.results$co)[j], sep = ""), 
                  " = 0", sep = "")
                for (i in 1:(length(egv.maps) - 1)) calc <- paste(calc, pre1[i], sep = " + ")
                calc <- paste(calc, " + ", pre1[length(egv.maps)], "'", sep = "")
                system(calc)
            }
        }
    }
    cat("\n\nCalculating Mahalanobis distances ...\n\n")
    tmp <- as.matrix(enirg.results$obs.li[, c(3:(enirg.results$nf + 3), 2)])
    mat.presences <- matrix(nrow = sum(tmp[, enirg.results$nf + 2]), ncol = ncol(enirg.results$co))
    colnames(mat.presences) <- colnames(enirg.results$co)
    samples <- 1:nrow(tmp)
    position <- 1
    for (i in samples) {
        variables <- c(tmp[i, 1:(enirg.results$nf + 1)])
        rep <- tmp[i, enirg.results$nf + 2]
        ocupa <- (position):(position + rep - 1)
        for (j in ocupa) {
            mat.presences[j, 1:(enirg.results$nf + 1)] <- variables
        }
        position <- position + rep
    }
    mean.centre <- apply(mat.presences, 2, median)
    cov <- t(as.matrix(mat.presences)) %*% as.matrix(mat.presences)/nrow(mat.presences)
    cova <- solve(cov)
    cat(paste("\n\nCalculating the Habitat Suitability Map (", enirg.results$species, 
        "_hsm...\n\n", sep = ""))
    f3 <- function(i) enirg.results$co[, i]/sqrt(crossprod(enirg.results$co[, i])/length(enirg.results$co[, 
        i]))
    c1 <- matrix(unlist(lapply(1:(enirg.results$nf + 1), f3)), length(enirg.results$egvs))
    if (method == "normal") {
        li.maps <- LETTERS[1:(enirg.results$nf + 1)]
        li_number <- length(li.maps)
        c1 <- paste("(", li.maps, "-", mean.centre, ")", sep = "")
        c <- (paste(c1, "*", cova, sep = ""))
        c3 <- NULL
        for (i in 1:(enirg.results$nf + 1)) {
            c3[i] <- ""
            for (j in 1:enirg.results$nf) c3[i] <- paste(c3[i], c[j + (i - 1) * (enirg.results$nf + 
                1)], "+", sep = "")
            j <- j + 1
            c3[i] <- paste(c3[i], c[j + (i - 1) * (enirg.results$nf + 1)], sep = "")
        }
        calc <- paste("(", c3[1], ")", "*", c1[1], sep = "")
        i <- 1
        if (length(c3) > 2) 
            for (i in 2:(length(c3) - 1)) calc <- paste(calc, "+", "(", c3[i], ")", 
                "*", c1[i], sep = "")
        calc <- paste(calc, "+", "(", c3[length(c3)], ")", "*", c1[i + 1], sep = "")
        if (is.null(egv.maps)) {
            execGRASS("r.mapcalculator", amap = "li_Mar", bmap = "li_Spec1", cmap = "li_Spec2", 
                dmap = "li_Spec3", emap = "li_Spec4", fmap = "li_Spec5", formula = calc, 
                outfile = paste(enirg.results$species, "_pred", sep = ""), flags = "overwrite", 
                legacyExec = TRUE)
        }
        else execGRASS("r.mapcalculator", amap = "predli_Mar", bmap = "predli_Spec1", 
            cmap = "predli_Spec2", dmap = "predli_Spec3", emap = "predli_Spec4", 
            fmap = "predli_Spec5", formula = calc, outfile = paste(enirg.results$species, 
                "_pred", sep = ""), flags = "overwrite", legacyExec = TRUE)
        statistic <- execGRASS("r.univar", map = paste(enirg.results$species, "_pred", 
            sep = ""), flags = "g", intern = TRUE, legacyExec = TRUE)
        map.max <- as.numeric(strsplit(statistic[agrep(statistic, pattern = "max")], 
            "=")[[1]][[2]])
        execGRASS("r.mapcalculator", amap = paste(enirg.results$species, "_pred", 
            sep = ""), formula = paste("1-(A/", map.max, ")", sep = ""), outfile = paste(enirg.results$species, 
            "_hsm", sep = ""), flags = "overwrite", legacyExec = TRUE)
        cat(paste(enirg.results$species, "_hsm", sep = ""), " HSM map was sucessfully created!\n\n", 
            sep = "")
        if (load.map == TRUE) {
            HSM[[paste(enirg.results$species, "_hsm", sep = "")]] <- raster(readRAST6(paste(enirg.results$species, 
                "_hsm", sep = "")))
        }
    }
    if (method == "large") {
        if (is.null(egv.maps)) 
            predli <- paste("li_", colnames(enirg.results$co), sep = "")
        else predli <- paste("predli_", colnames(enirg.results$co), sep = "")
        c1 <- paste("(", predli, "-", mean.centre, ")", sep = "")
        c <- (paste(c1, "*", cova, sep = ""))
        c3 <- NULL
        for (i in 1:(enirg.results$nf + 1)) {
            c3[i] <- ""
            for (j in 1:enirg.results$nf) c3[i] <- paste(c3[i], c[j + (i - 1) * (enirg.results$nf + 
                1)], "+", sep = "")
            j <- j + 1
            c3[i] <- paste(c3[i], c[j + (i - 1) * (enirg.results$nf + 1)], sep = "")
        }
        calc <- paste("r.mapcalc '", enirg.results$species, "_pred=(", c3[1], ")", 
            "*", c1[1], sep = "")
        i <- 1
        if (length(c3) > 2) 
            for (i in 2:(length(c3) - 1)) calc <- paste(calc, "+", "(", c3[i], ")", 
                "*", c1[i], sep = "")
        calc <- paste(calc, "+", "(", c3[length(c3)], ")", "*", c1[i + 1], "'", sep = "")
        system(calc)
        statistic <- system(paste("r.univar -g ", paste(enirg.results$species, "_pred", 
            sep = ""), sep = ""), intern = T)
        map.max <- as.numeric(strsplit(statistic[agrep(statistic, pattern = "max")], 
            "=")[[1]][[2]])
        system(paste("r.mapcalc '", paste(enirg.results$species, "_hsm", sep = ""), 
            " = 1 - (", paste(enirg.results$species, "_pred", sep = ""), "/", map.max, 
            ")'", sep = ""))
        cat(paste(paste(enirg.results$species, "_hsm", sep = ""), " HSM map was sucessfully created!\n\n", 
            sep = ""))
        cat("You can find HSM in your mapset in GRASS. You can see it through QGIS or using 'raster' library instead.")
    }
    execGRASS("g.remove", vect = enirg.results$species)
    presences_map <- SpatialPointsDataFrame(enirg.results$presences, data = data.frame(presences = enirg.results$presences[, 
        3]))
    writeVECT6(presences_map, enirg.results$species, v.in.ogr_flags = c("e", "overwrite", 
        "o", "z"))
    execGRASS("v.db.addcol", map = enirg.results$species, columns = "pred double precision")
    execGRASS("v.what.rast", vector = enirg.results$species, raster = paste(enirg.results$species, 
        "_hsm", sep = ""), column = "pred")
    vect <- readVECT6(enirg.results$species)
    vect <- data.frame(cbind(vect@data$presences, vect@data$pred))
    names(vect) <- c("observed", "predicted")
    HSM[["predictions"]] <- vect
    hsm.report <- strsplit(execGRASS("r.report", map = paste(enirg.results$species, 
        "_hsm", sep = ""), nsteps = 10, units = "c", intern = T)[16:25], split = " ")
    hist.hsm <- c()
    for (i in 1:10) {
        hist.hsm <- c(hist.hsm, as.numeric(gsub("([0-9]+).*$", "\\1", hsm.report[[i]][length(hsm.report[[i]])])))
    }
    intervals <- seq(0, 1, by = 0.1)
    hist.predicted <- hist(vect$predicted, breaks = intervals, plot = FALSE)$counts
    HSM[["validation"]] <- cbind(intervals[-11], intervals[-1], hist.hsm, hist.predicted)
    colnames(HSM[["validation"]])[1:2] <- c("int.inf", "int.sup")
    class(HSM) <- c("hsm")
    return(HSM)
}
