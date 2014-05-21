list.maps <-
function(type = c("rast", "vect"), prefix = "*") {
    maps.list <- NULL
    for (i in type) {
        maps.list[[i]] <- execGRASS("g.mlist", type = i, pattern = prefix, intern = TRUE)
    }
    return(maps.list)
}
