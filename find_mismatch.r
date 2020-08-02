tested <- read.csv("tested-instances-list.csv", header=F)
optimal <- read.csv("optimal-values.csv", header=F, sep=";")

found_mismatch <- FALSE

for (i in 1:dim(tested)[1]) {
    found_mismatch <- TRUE
    for (j in 1:dim(optimal)[1]) {
        if (tested[i,1] == optimal[j,2]) {
            found_mismatch <- FALSE
        }
    }
    if (found_mismatch) {
        print(paste("Missing optimal of",tested[i,1]))
    }
}