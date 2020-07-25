#! COMPUTE AGGREGATE MEASURES OF A TEST-BED FROM THE TEST_RESULTS.csv FILE

# Read test results data and set column names
data <- read.csv("test_results_nogap.csv", header=T, sep=";")
colnames(data) <- c("Instance","Seed","Cost","Fractionality","Rounds","LPtime(ms)","ZItime(ms)","LP+ZItime(ms)")

#! Compute success rate = number of rows with zero fractionality / total number of rows

# Get fractionalities
fracs <- data$Fractionality
# Compute number of successes and finally the success rate (%)
num_succ <- sum(fracs == 0)
succ_rate <- num_succ / length(fracs) * 100

#! Compute the shifted geometric means (SGM) of the LP solve and ZI-Round times, and their ratio ZI/LP

# Get LP solve times and ZI-Round times
lptimes <- data$LPtime
zitimes <- data$ZItime
# Shifts for the shifted geometric means (in milliseconds)
lpshift <- 1000
zishift <- 100
# Apply shifts
shifted_lptimes <- lptimes + lpshift
shifted_zitimes <- zitimes + zishift
# Compute geometric mean
geom_lptime <- exp(mean(log(shifted_lptimes)))
geom_zitime <- exp(mean(log(shifted_zitimes)))
# Revert initial shifts
sgm_lptime <- geom_lptime - lpshift
sgm_zitime <- geom_zitime - zishift
# Compute shifted geometric means ratio ZI/LP (%)
sgm_ratio <- sgm_zitime / sgm_lptime * 100

#! Compute the gap column and the average gap w.r.t. the optimal solutions

#Add the gap column to the test results dataframe
append(colnames(data),"Gap(%)")
# Read optimal solutions data and set column names
optimal <- read.csv("optimal-miplib2003.csv", header=F, sep=";")
colnames(optimal) <- c("Flag","Name","Opt")
# Build unnamed data
named <- optimal[, c("Flag","Opt")]
# Set names of the instances as names of the rows
rownames(named) <- optimal[,"Name"]
# Scan rows of test results and optimal solutions to
# add the gap column to the test results dataframe
for (i in 1:dim(data)[1]) {
    inst <- data[i,"Instance"]
    if (fracs[i] == 0) {
        # ZI-Round succeeded: compute actual gap
        flag <- optimal[i,"Flag"]
        if (flag == "opt" || flag == "best") {
            # Optimal or best solution exists
            opt <- optimal[i,"Opt"]
            gap <- (data[i,"Cost"] - opt) / abs(opt) * 100
            data[i,"Gap(%)"] <- min(round(gap, digits=2), 100)
        } else {
            # No optimal or best solution exists
            data[i,"Gap(%)"] <- 0
        }
    } else {
        # ZI-Round failed: gap is 100%
        data[i,"Gap(%)"] <- 100
    }
}
# Compute the average gap
avg_gap <- mean(data[,"Gap(%)"])
# Print the complete test results file
write.table(data, file="test_results.csv", row.names=FALSE, dec=".", sep=";", quote=FALSE)

#! Print aggregate measures to file
names <- c("SuccessRate(%)","SGM-LPtime(ms)","SGM-ZItime(ms)","SGM-ratio(ZI/LP)","AvgGap(%)")
aggr <- c(succ_rate, sgm_lptime, sgm_zitime, sgm_ratio, avg_gap)
df <- t(data.frame(names, aggr))
write.table(df, file="aggregate_measures.csv", row.names=FALSE, col.names=FALSE, dec=".", sep=";", quote=FALSE)


#! Print aggregate measures to video
print(paste("ZI-Round success rate (%):", succ_rate, "%"), quote=FALSE)
print(paste("SGM LP solve time (ms):", sgm_lptime), quote=FALSE)
print(paste("SGM ZI-Round time (ms):", sgm_zitime), quote=FALSE)
print(paste("SGM ratio ZI/LP (%):", sgm_ratio, "%"), quote=FALSE)
print(paste("ZI-Round avg opt gap(%):", avg_gap, "%"), quote=FALSE)