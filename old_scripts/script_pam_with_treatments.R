#take cleaned data from Python script and run through phytotools package to get ek and alpha etc
library("phytotools", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")

alpha_ek_alga <- read.csv("/Users/Angela/Library/Mobile Documents/com~apple~CloudDocs/research_plakobranchus/plako_paper/crj/plako_data_122021_for_final_revised/plako_cleaned.csv", sep = ";")

# add a new column with date and specimen ID as unique key
alpha_ek_alga$uid <- paste(alpha_ek_alga$Date,alpha_ek_alga$ID, sep = "_")

# add a new column with a number to indicate just the plant ID (temperature), to be used to color the lines in the plot
alpha_ek_alga$temp = substr(alpha_ek_alga$ID, 2, 2)

# add a new column with a number to indicate just the fragment ID (treatment), to be used to color the lines in the plot
alpha_ek_alga$treatment = substr(alpha_ek_alga$ID, 1, 2)

# remove all rows where rETR is null or negative
selectedData <- subset(alpha_ek_alga, rETR > 0) 

# load up the unique row IDs into a vector
uniqueIds <- unique(selectedData$uid)

n <- length(uniqueIds)

#add a new column for rETRmax
rETRMaxes = array(NA,c(n,1))
#temperatures = array(dim = n)
for (i in 1:n){
  sub <- subset(selectedData, uid == uniqueIds[i])
  subMaxETR <- max(sub$rETR)
  selectedData$rETRmax[selectedData$uid == uniqueIds[i]] <- subMaxETR
  rETRMaxes[i] = subMaxETR
  #temperatures[i] = max(sub$Temp)
}

# prepare empty matrixes to hold output from fitWebb
alpha <- array(NA,c(n,4))
ek <- array(NA,c(n,4))

# Uncomment these to use quartz, if needed
# quartz(width = 24, height = 12)
# par(mfrow = c(7, 9),  lwd = 1.3)

# Use this to switch from one plot with all the curves overlaid together to one individual plot for each specimen
individual_plots = FALSE

if (!individual_plots)plot
plot(NA,NA, main = "rETR", 
     xlab = "Epar (μmols photons m-2 s-1)", 
     ylab = "rETR (μmols electrons m-2 s-1)", 
     xlim = c(0,1000), 
     ylim = c(0,400))

for (i in 1:n){
  #Get ith data
  Epar <- selectedData$Epar[selectedData$uid==uniqueIds[i]]
  rETR <- selectedData$rETR[selectedData$uid==uniqueIds[i]]
  
  #Call function
  myfit <- fitWebb(Epar,rETR)
  #Store the 4-values outputs into the matrix 
  alpha[i,] <- myfit$alpha
  ek[i,] <- myfit$ek 
  
  # plot the data points
  if (individual_plots)
    plot(Epar,rETR, main = uniqueIds[i], xlab = "Epar (μmols photons m-2 s-1)", ylab = "rETR (μmols electrons m-2 s-1)", xlim = c(0,1000), ylim = c(0,500))
  
  E <- seq(0,1000,by=1)
  with(myfit, {
    P <- alpha[1] * ek[1] * (1 - exp (-E / ek[1]))
    lines(E, P)
  })
}

# attach the unique IDs, the alpha and the ek result sets side-by-side, columns with date and specimen ID
dates = substr(uniqueIds, 1, 10)
specimens = substr(uniqueIds, 12, 16)

#add a new column for ETRmax and calculate from rETRmax
ETRmax = round(rETRMaxes*0.5*0.84, digits = 2)

#add a new column for ETRmax based on rETRmax for Dr. Smith to vizualize (this will add to selectedData)
selectedData$ETRmax <- round(selectedData$rETRmax*0.5*0.84, digits = 2)

#add a new column for treatments assigned from the 16th number in uid
#C 35‰/14μM", "28‰/27μM", "18‰/53μM", "11‰/80μM (kept for ulva/hypnea)
treatment = c("C 35‰", "28‰", "18‰", "11‰")[as.integer(substr(uniqueIds, 16, 16))]

#add a new column for species assigned from the 12th character in uid
species = as.character(substr(uniqueIds, 12, 13))

rlc_days_by_date = array(dim = length(dates))
hash_of_rlc_days_by_date <- hash(rlc_day_assign$Date, rlc_day_assign$RLC.Day)
for (i in 1:length(dates)) {
  rlc_days_by_date[i] = hash_of_rlc_days_by_date[[dates[i]]]
}
lanai_side_by_date = array(dim = length(dates))
hash_of_lanai_side_by_date <- hash(rlc_day_assign$Date, rlc_day_assign$Lanai.side)
for (i in 1:length(dates)) {
  lanai_side_by_date[i] = as.character(hash_of_lanai_side_by_date[[dates[i]]])
}

output <- cbind(
  #cbind(
    #cbind(
      #cbind(
        #cbind(
          #cbind(
            cbind(
              cbind(
                cbind(
                  cbind(
                    cbind(
                      cbind(
                      dates, 
                     specimens), 
                   uniqueIds), 
                 round(alpha, digits = 3)), 
               round(ek, digits = 1)), 
             rETRMaxes), 
           substr(specimens, 1, 3)), 
         ETRmax) 
       #treatment),
     #species)
   #rlc_days_by_date), 
 #lanai_side_by_date),
#temperatures)

# convert to a data frame
output_df <- data.frame(output)
# assign the names to the columns
names(output_df) <- c("Date", "Specimen ID", "uid", "alpha_est", "alpha_stderr", "alpha_t", 
                      "alpha_p", "ek_est", "ek_stderr", "ek_t", "ek_p", "rETRmax", "plant ID", "ETRmax")#, "Treatment", "Species", "RLC Day", "Lanai Side", "Temperature (°C)")

# save to file
write.csv(output_df, "/Users/Angela/Library/Mobile Documents/com~apple~CloudDocs/research_plakobranchus/plako_paper/crj/plako_data_122021_for_final_revised/pam_files/plako_alpha_ek_rounded.csv")








