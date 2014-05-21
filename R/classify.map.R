classify.map <-
function(enirg.results, suit.classes, output.name = NULL, load.map = FALSE) {
    if (class(enirg.results) != "enirg") 
        stop("The function classify.map needs an object of class 'enirg'!")
    if (class(suit.classes) != "CBI") 
        stop("The function classify.map needs an object of class 'CBI'!")
    if (is.null(output.name)) 
        output.name <- paste(enirg.results$species, "_hsm_class", sep = "")
    breaks <- suit.classes$intervals[2:4]
    print(breaks)
    settings_text <- " 0 thru x1x = 1 unsuitable\n x2x thru x3x = 2 marginal\n x4x thru x5x = 3 suitable\n x6x thru 1 = 4 optimal"
    settings_text <- sub("x1x", breaks[1], settings_text)
    settings_text <- sub("x2x", breaks[1], settings_text)
    settings_text <- sub("x3x", breaks[2], settings_text)
    settings_text <- sub("x4x", breaks[2], settings_text)
    settings_text <- sub("x5x", breaks[3], settings_text)
    settings_text <- sub("x6x", breaks[3], settings_text)
    write(settings_text, file = "classes.txt")
    execGRASS("r.reclass", input = paste(enirg.results$species, "_hsm", sep = ""), 
        output = output.name, rules = "classes.txt", flags = "overwrite")
    cat(paste(output.name, " raster map was sucessfully classify!\n", sep = ""))
    if (load.map) {
        map <- list()
        map[[paste(enirg.results$species, "HSM", sep=" ")]] <- raster(readRAST6(output.name))
        class(map) <- c("hsm", "classified")
        return(map)
    }
}
