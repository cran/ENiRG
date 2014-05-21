enirg <-
function(presences.table, qtegv.maps, qlegv.maps = NULL, col.w = NULL, scannf = TRUE, 
    nf = 1, method = "normal", load.maps = FALSE, species.name = "species") {
    egv.maps <- c(qtegv.maps, qlegv.maps)
    number_maps <- length(egv.maps)
    execGRASS("g.region", rast = egv.maps[1])
    region.parameters <- gmeta6()
    null_cells <- execGRASS("r.stats", input = egv.maps[1], flags = "c", intern = TRUE)
    null_cells <- agrep(null_cells, pattern = "*", max.distance = list(all = 0), 
        value = T)
    null_cells <- as.numeric(strsplit(null_cells, " ")[[1]][2])
    numberpixel <- as.numeric(region.parameters$cells) - null_cells
    cat(paste("The number of pixels for the calculated area is", numberpixel, "\n\n", 
        sep = " "))
    execGRASS("g.remove", vect = species.name)
    presences_map <- SpatialPointsDataFrame(presences.table, data = data.frame(presences = presences.table[, 
        3]))
    writeVECT6(presences_map, species.name, v.in.ogr_flags = c("e", "overwrite", 
        "o", "z"))
    execGRASS("v.to.rast", parameters = list(input = species.name, output = species.name, 
        use = "attr", column = "presences"), flags = "overwrite")
    cat(paste(species.name, " was sucessfully convert into a raster map and loaded into GRASS ...\n\n", 
        sep = ""))
    presences_map <- species.name
    number_presences <- sum(presences.table[, 3])
    combine_maps <- t(combn(1:number_maps, 2))
    map.cor <- cbind(egv.maps[combine_maps[, 1]], egv.maps[combine_maps[, 2]])
    combine_maps2 <- rbind(t(combn(1:number_maps, 2)), cbind(1:number_maps, 1:number_maps))
    map.cor2 <- cbind(egv.maps[combine_maps2[, 1]], egv.maps[combine_maps2[, 2]])
    if (method == "normal") {
        cat("\n\nCalculating the species covariance matrix (Rs) ...\n\n")
        grass.cor <- NULL
        for (i in 1:nrow(map.cor)) {
            tmp <- execGRASS("r.stats", input = paste(map.cor[i, 1], map.cor[i, 2], 
                presences_map, sep = ","), flags = c("1", "A"), intern = TRUE, legacyExec = TRUE)
            wrapping <- function(tmp) as.numeric(strsplit(tmp, " ")[[1]])
            tmp <- t(apply(matrix(tmp[-1 * agrep(tmp, pattern = "* *")], ncol = 1), 
                1, wrapping))
            tmp <- tmp[which(tmp[, 3] > 0), ]
            np <- sum(tmp[, 3])
            tmp <- sum(apply(tmp, 1, prod))
            grass.cor[i] <- tmp/np
        }
        grass.cor2 <- NULL
        for (i in 1:number_maps) {
            tmp <- execGRASS("r.stats", input = paste(egv.maps[i], egv.maps[i], presences_map, 
                sep = ","), flags = c("1", "A"), intern = TRUE, legacyExec = TRUE)
            wrapping <- function(tmp) as.numeric(strsplit(tmp, " ")[[1]])
            tmp <- t(apply(matrix(tmp[-1 * agrep(tmp, pattern = "* *")], ncol = 1), 
                1, wrapping))
            tmp <- tmp[which(tmp[, 3] > 0), ]
            np <- sum(tmp[, 3])
            tmp <- sum(apply(tmp, 1, prod))
            grass.cor2[i] <- tmp/np
        }
        Rs <- matrix(1, number_maps, number_maps)
        diag(Rs) <- grass.cor2
        Rs[t(combn(1:number_maps, 2))] <- grass.cor
        Rs[t(combn(1:number_maps, 2))[, 2:1]] <- grass.cor
        rownames(Rs) <- colnames(Rs) <- egv.maps
        cat("\n\nCalculating Ze ...\n\n")
        grass.cor = NULL
        for (i in 1:nrow(map.cor)) {
            tmp <- execGRASS("r.stats", input = paste(map.cor[i, 1], map.cor[i, 2], 
                sep = ","), flags = c("1", "A"), intern = TRUE, legacyExec = TRUE)
            wrapping <- function(tmp) as.numeric(strsplit(tmp, " ")[[1]])
            tmp <- t(apply(matrix(tmp[-1 * agrep(tmp, pattern = "* *")], ncol = 1), 
                1, wrapping))
            grass.cor[i] <- sum(apply(tmp, 1, prod))
        }
        Ze <- matrix(1, number_maps, number_maps)
        Ze[t(combn(1:number_maps, 2))] <- grass.cor
        Ze[t(combn(1:number_maps, 2))[, 2:1]] <- grass.cor
        Ze <- Ze/numberpixel
        rownames(Ze) <- colnames(Ze) <- egv.maps
        cat("\n\nCalculating the coordinates of the marginality vector ...\n\n")
        mar <- NULL
        for (i in egv.maps) {
            tmp <- execGRASS("r.stats", input = paste(i, presences_map, sep = ","), 
                flags = c("1"), intern = TRUE, legacyExec = TRUE)
            wrapping <- function(tmp) as.numeric(strsplit(tmp, " ")[[1]])
            tmp <- (t(apply(matrix(tmp[-1 * agrep(tmp, pattern = "* *")], ncol = 1), 
                1, wrapping)))
            tmp <- tmp[which(tmp[, 2] > 0), ]
            np <- sum(tmp[, 2])
            tmp <- apply(tmp, 1, prod)
            mar[i] <- sum(tmp)/np
        }
        names(mar) <- egv.maps
        cat("\n\nCalculating the matrix Ge ...\n\n")
        grass.cor <- NULL
        for (i in 1:nrow(map.cor)) {
            tmp <- execGRASS("r.stats", input = paste(map.cor[i, 1], map.cor[i, 2], 
                sep = ","), flags = "1", intern = TRUE, legacyExec = TRUE)
            wrapping <- function(tmp) as.numeric(strsplit(tmp, " ")[[1]])
            tmp <- t(apply(matrix(tmp[-1 * agrep(tmp, pattern = "* *")], ncol = 1), 
                1, wrapping))
            tmp <- sum(apply(tmp, 1, prod))
            grass.cor[i] <- tmp/numberpixel
        }
        Ge <- matrix(1, number_maps, number_maps)
        Ge[t(combn(1:number_maps, 2))] <- grass.cor
        Ge[t(combn(1:number_maps, 2))[, 2:1]] <- grass.cor
        rownames(Ge) <- colnames(Ge) <- egv.maps
        cat("\n\nCalculating the matrix Se ...\n\n")
        for (i in egv.maps) {
            output.names <- tmp <- execGRASS("r.mapcalculator", amap = i, bmap = presences_map, 
                formula = paste("A*B", "/", number_presences, sep = ""), outfile = paste("pond_", 
                  i, sep = ""), flags = "overwrite", legacyExec = TRUE)
        }
        grass.cor <- NULL
        for (i in 1:nrow(map.cor2)) {
            tmp <- execGRASS("r.stats", input = paste(map.cor2[i, 1], ",pond_", map.cor2[i, 
                2], sep = ""), flags = c("1", "A"), intern = TRUE, legacyExec = TRUE)
            wrapping <- function(tmp) as.numeric(strsplit(tmp, " ")[[1]])
            tmp <- t(apply(matrix(tmp[-1 * agrep(tmp, pattern = "* *")], ncol = 1), 
                1, wrapping))
            grass.cor[i] <- sum(apply(tmp, 1, prod))
        }
        Se <- matrix(1, number_maps, number_maps)
        Se[combine_maps2] <- grass.cor
        Se[combine_maps2[, 2:1]] <- grass.cor
        rownames(Se) <- colnames(Se) <- egv.maps
    }
    if (method == "large") {
        cat("\n\nCalculating the species covariance matrix (Rs) ...\n\n")
        grass.cor.map <- NULL
        for (i in 1:nrow(map.cor)) grass.cor.map[i] <- paste("r.stats -1 -A ", map.cor[i, 
            1], map.cor[i, 2], presences_map, sep = ",")
        grass.cor.p <- NULL
        for (i in 1:nrow(map.cor)) grass.cor.p[i] <- as.numeric(system(paste(grass.cor.map[i], 
            " | awk '$3>0 {sum=sum+($1*$3*$2);np=np+$3} END{print sum/np}'"), intern = T))
        grass.cor.map2 <- NULL
        for (i in 1:number_maps) grass.cor.map2[i] <- paste("r.stats -1 -A ", egv.maps[i], 
            egv.maps[i], presences_map, sep = ",")
        grass.cor.p2 <- NULL
        for (i in 1:number_maps) grass.cor.p2[i] <- as.numeric(system(paste(grass.cor.map2[i], 
            " | awk '$3>0 {sum=sum+($1*$3*$2);np=np+$3} END{print sum/np}'"), intern = T))
        Rs <- matrix(1, number_maps, number_maps)
        diag(Rs) <- grass.cor.p2
        Rs[t(combn(1:number_maps, 2))] <- grass.cor.p
        Rs[t(combn(1:number_maps, 2))[, 2:1]] <- grass.cor.p
        rownames(Rs) <- colnames(Rs) <- egv.maps
        cat("\n\nCalculating Ze ...\n\n")
        grass.cor.map <- NULL
        for (i in 1:nrow(map.cor)) grass.cor.map[i] <- paste("r.stats -1 -A ", map.cor[i, 
            1], map.cor[i, 2], sep = ",")
        grass.cor.p <- NULL
        for (i in 1:nrow(map.cor)) grass.cor.p[i] <- as.numeric(system(paste(grass.cor.map[i], 
            " | awk '{sum=sum+($1*$2)} END{print sum}'"), intern = T))
        Ze <- matrix(1, number_maps, number_maps)
        Ze[t(combn(1:number_maps, 2))] <- grass.cor.p
        Ze[t(combn(1:number_maps, 2))[, 2:1]] <- grass.cor.p
        Ze <- Ze/numberpixel
        rownames(Ze) <- colnames(Ze) <- egv.maps
        cat("\n\nCalculating the coordinates of the marginality vector ...\n\n")
        mar <- NULL
        for (i in egv.maps) mar[i] <- as.numeric(system(paste("r.stats -1 ", i, ",", 
            presences_map, " | awk '{s=s+$1*$2;n=n+$2} END{print s/n}'", sep = ""), 
            intern = T))
        names(mar) <- egv.maps
        cat("\n\nCalculating the matrix Ge ...\n\n")
        grass.cor.map <- NULL
        for (i in 1:nrow(map.cor)) grass.cor.map[i] <- paste("r.stats -1 ", map.cor[i, 
            1], map.cor[i, 2], sep = ",")
        grass.cor.p <- NULL
        for (i in 1:nrow(map.cor)) grass.cor.p[i] <- as.numeric(system(paste(grass.cor.map[i], 
            " | awk '{sum=sum+($1*$2)/", numberpixel, "} END{print sum}'"), intern = T))
        Ge <- matrix(1, number_maps, number_maps)
        Ge[t(combn(1:number_maps, 2))] <- grass.cor.p
        Ge[t(combn(1:number_maps, 2))[, 2:1]] <- grass.cor.p
        rownames(Ge) <- colnames(Ge) <- egv.maps
        cat("\n\nCalculating the matrix Se ...\n\n")
        for (i in egv.maps) system(paste("r.mapcalc 'pond_", i, " = ", i, " *", presences_map, 
            " / ", number_presences, "'", sep = ""))
        grass.cor.map <- NULL
        for (i in 1:nrow(map.cor2)) grass.cor.map[i] <- paste("r.stats -1A ", map.cor2[i, 
            1], ",pond_", map.cor2[i, 2], sep = "")
        grass.cor.p <- NULL
        for (i in 1:nrow(map.cor2)) grass.cor.p[i] <- as.numeric(system(paste(grass.cor.map[i], 
            " | awk '{sum=sum+($1*$2)} END{print sum}'"), intern = T))
        Se <- matrix(1, number_maps, number_maps)
        Se[combine_maps2] <- grass.cor.p
        Se[combine_maps2[, 2:1]] <- grass.cor.p
        rownames(Se) <- colnames(Se) <- egv.maps
    }
    eS <- eigen(Se)
    kee <- (eS$values > 1e-09)
    S12 <- eS$vectors[, kee] %*% diag(eS$values[kee]^(-0.5)) %*% t(eS$vectors[, kee])
    W <- S12 %*% Ge %*% S12
    if (is.null(col.w)) 
        col.w <- rep(1, number_maps)
    me <- mar * sqrt(col.w)
    x <- S12 %*% me
    b <- x/sqrt(sum(x^2))
    H <- (diag(ncol(Ze)) - b %*% t(b)) %*% W %*% (diag(ncol(Ze)) - b %*% t(b))
    s <- eigen(H)$values[-number_maps]
    if (scannf) {
        barplot(s)
        cat("Select the number of specialization axes: ")
        nf <- as.integer(readLines(n = 1))
    }
    if (nf <= 0 | nf > (ncol(Ze) - 1)) 
        nf <- 1
    co <- matrix(nrow = number_maps, ncol = nf + 1)
    tt <- data.frame((S12 %*% eigen(H)$vectors)[, 1:nf])
    ww <- apply(tt, 2, function(x) x/sqrt(col.w))
    norw <- sqrt(diag(t(as.matrix(tt)) %*% as.matrix(tt)))
    co[, 2:(nf + 1)] <- sweep(ww, 2, norw, "/")
    m <- me/sqrt(col.w)
    co[, 1] <- m/sqrt(sum(m^2))
    m <- sum(m^2)
    marginalities <- matrix(co[, 1], ncol = 1)
    rownames(marginalities) <- egv.maps
    colnames(marginalities) <- "Marginality"
    total_marginality <- sqrt(sum((marginalities/sqrt(col.w))^2))/1.96
    rownames(co) <- egv.maps
    colnames(co) <- c("Mar", paste("Spec", 1:nf, sep = ""))
    specializations <- matrix(co[, 2:(nf + 1)], ncol = nf)
    rownames(specializations) <- egv.maps
    colnames(specializations) <- paste("Spec", 1:nf, sep = "")
    total_specialization <- sqrt(sum(abs(s)))/number_maps
    call <- match.call()
    results_ENFA <- list(call = call, nf = nf, cw = col.w, species = species.name, 
        egvs = egv.maps, qt.egvs = qtegv.maps, ql.egvs = qlegv.maps, presences = presences.table, 
        total.marginality = total_marginality, marginalities = marginalities, total.specialization = total_specialization, 
        specializations = specializations, co = co, mar = mar, m = m, s = s)
    class(results_ENFA) <- c("enirg")
    print(results_ENFA)
    cat("\n\nCalculating li maps ...\n\n")
    if (method == "normal") {
        for (j in 1:(nf + 1)) {
            execGRASS("r.mapcalculator", outfile = paste("li_", colnames(co)[j], 
                sep = ""), formula = "0", flags = "overwrite")
            for (i in 1:number_maps) {
                execGRASS("r.mapcalculator", amap = egv.maps[i], outfile = "temp", 
                  formula = paste("A*", results_ENFA$co[i, j], sep = ""), flags = c("overwrite"))
                execGRASS("r.mapcalculator", amap = paste("li_", colnames(co)[j], 
                  sep = ""), bmap = "temp", outfile = paste("li_", colnames(co)[j], 
                  sep = ""), formula = "A+B", flags = "overwrite")
            }
        }
    }
    if (method == "large") {
        for (j in 1:(nf + 1)) {
            cc <- paste("(", co[, j], ")", sep = "")
            pre1 <- paste(egv.maps, cc, sep = " * ")
            calc <- paste("r.mapcalc '", paste("li_", colnames(co)[j], sep = ""), 
                " = 0", sep = "")
            for (i in 1:(number_maps - 1)) calc <- paste(calc, pre1[i], sep = " + ")
            calc <- paste(calc, " + ", pre1[number_maps], "'", sep = "")
            system(calc)
        }
        load.maps <- FALSE
    }
    if (load.maps) {
        for (i in 1:(nf + 1)) {
            results_ENFA[[paste("li_", colnames(co)[i], sep = "")]] <- raster(readRAST6(paste("li_", 
                colnames(co)[i], sep = "")))
        }
        plot(results_ENFA$li_Mar)
        contour(results_ENFA$li_Mar, add = T)
    }
    for (i in colnames(co)) {
        execGRASS("v.db.addcol", map = results_ENFA$species, columns = paste("li_", 
            i, " double precision", sep = ""))
        execGRASS("v.what.rast", vector = results_ENFA$species, raster = paste("li_", 
            i, sep = ""), column = paste("li_", i, sep = ""))
    }
    vectorial.data <- readVECT6(results_ENFA$species)
    results_ENFA[["obs.li"]] <- vectorial.data@data
    cat("\n\nCalculations done!\n\n")
    return(results_ENFA)
}
