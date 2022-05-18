#multiple dataset from data removing empty rETR lines and renaming working dataset to ae
#can be used for other algae with simple treatments

library("phytotools", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")

alpha_ek_hypnea <- read.csv("//Users/Angela/Library/Mobile Documents/com~apple~CloudDocs/research_limu/ulva_downloaded_data_from_drive/output/ulva-cleaned.csv", sep = ";")

# add a new column with date and specimen ID as unique key
alpha_ek_hypnea$uid <- paste(alpha_ek_hypnea$Date,alpha_ek_hypnea$ID, sep = "_")


# add a new column with a number to indicate just the plant ID, to be used to color the lines in the plot
#alpha_ek_hypnea$plantID = substr(alpha_ek_hypnea$ID, 3, 3)

# remove all rows where rETR is null or negative
selectedData <- subset(alpha_ek_hypnea, rETR > 0) 

# load up the unique row IDs into a vector
uniqueIds <- unique(selectedData$uid)

n <- length(uniqueIds)

#add a new column for rETRmax
rETRMaxes = array(NA,c(n,1))
for (i in 1:n){
  sub <- subset(selectedData, uid == uniqueIds[i])
  subMaxETR <- max(sub$rETR)
  selectedData$rETRmax[selectedData$uid == uniqueIds[i]] <- subMaxETR
  rETRMaxes[i] = subMaxETR
}

# prepare empty matrixes to hold output from fitWebb
alpha <- array(NA,c(n,4))
ek <- array(NA,c(n,4))

# Uncomment these to use quartz, if needed
# quartz(width = 24, height = 12)
# par(mfrow = c(7, 9),  lwd = 1.3)

# Use this to switch from one plot with all the curves overlaid together to one individual plot for each specimen
individual_plots = TRUE

if (!individual_plots)plot
plot(NA,NA, main = "rETR", 
     xlab = "Epar (μmols photons m-2 s-1)", 
     ylab = "rETR (μmols electrons m-2 s-1)", 
     xlim = c(0,1000), 
     ylim = c(0,600))

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

# attach the unique IDs, the alpha and the ek result sets side-by-side
dates = substr(uniqueIds, 1, 10)
specimens = substr(uniqueIds, 12, 16)
#add a new column for ETRmax and calculate from rETRmax
ETRmax = round(rETRMaxes*0.5*0.84, digits = 2)

output <- cbind(cbind(cbind(dates, cbind(specimens, cbind(cbind(cbind(uniqueIds, round(alpha, digits = 3)), round(ek, digits = 1)))), rETRMaxes), substr(specimens, 3,3), ETRmax)) 

# convert to a data frame
output_df <- data.frame(output)
# assign the names to the columns
names(output_df) <- c("Date", "Specimen ID", "uid", "alpha_est", "alpha_stderr", "alpha_t", 
                      "alpha_p", "ek_est", "ek_stderr", "ek_t", "ek_p", "rETRmax", "plant ID", "ETRmax")

# save to file
write.csv(output_df, "/Users/Angela/Library/Mobile Documents/com~apple~CloudDocs/research_limu/ulva_downloaded_data_from_drive/output/ulva_alpha_ek_rounded.csv")








