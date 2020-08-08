#! COMPUTE AGGREGATE MEASURES OF A TEST-BED FROM THE TEST_RESULTS.csv FILE

# Test type string
testype <- "fractieworstobj"
# Random seed used (string)
rseed <- "37029"

# Read test results data and set column names
data <- read.csv(paste("test_results_nogap(",testype,")(seed_",rseed,").csv",sep=""), header=T, sep=";")
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
zishift <- 10
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
optimal <- read.csv("optimal-values.csv", header=F, sep=";")
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
        flag <- named[inst,"Flag"]
        if (flag == "opt" || flag == "best") {
            # Optimal or best solution exists
            opt <- named[inst,"Opt"]
            cost <- data[i,"Cost"]
            if (cost < opt) {
                data[i,"Gap(%)"] <- 0
                print(paste("Found cost",cost,"<",opt,"for",inst))
            } else {
                if (opt != 0) {
                    gap <- (cost - opt) / abs(opt) * 100
                    data[i,"Gap(%)"] <- min(round(gap, digits=2), 100)
                } else {
                     # Opt is zero
                    if (cost == 0) {
                        data[i,"Gap(%)"] <- 0
                    } else {
                        data[i,"Gap(%)"] <- 100
                    }
                    gap <- (cost - opt) / abs(opt) * 100
                }
            }
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
write.table(data, file=paste("test_results(",testype,")(seed_",rseed,").csv",sep=""), row.names=FALSE, dec=".", sep=";", quote=FALSE)

#! Print aggregate measures to file
names <- c("SuccessRate(%)","SGM-LPtime(ms)","SGM-ZItime(ms)","SGM-ratio(ZI/LP)","AvgGap(%)")
aggr <- c(succ_rate, sgm_lptime, sgm_zitime, sgm_ratio, avg_gap)
df <- t(data.frame(names, aggr))
write.table(df, file=paste("aggregate_measures(",testype,")(seed_",rseed,").csv",sep=""), row.names=FALSE, col.names=FALSE, dec=".", sep=";", quote=FALSE)


#! Print aggregate measures to video
print(paste("Test type:", testype), quote=FALSE)
print(paste("Seed:", rseed), quote=FALSE)
print(paste("ZI-Round success rate (%):", succ_rate, "%"), quote=FALSE)
print(paste("SGM LP solve time (ms):", sgm_lptime), quote=FALSE)
print(paste("SGM ZI-Round time (ms):", sgm_zitime), quote=FALSE)
print(paste("SGM ratio ZI/LP (%):", sgm_ratio, "%"), quote=FALSE)
print(paste("ZI-Round avg opt gap(%):", avg_gap, "%"), quote=FALSE)