time_start <- Sys.time()


resampleN <- 100  #number of resamples for the Monte Carlo error analysis; if you wish to bypass resampling, set to 1
percentile_values <- c(0.025, 0.16, 0.5, 0.84, 0.975) #list of percentiles (of any length) used to evaluate a resampled data set; the mean (0.5) is outputted by default and doesn't need to be included here
Godderis <- TRUE  #set to "TRUE" to run time arrays of fA, fAw/fA, fD, and GEOGset to "FALSE" to run standard time arrays from GEOCARBSULF.
input_distribution <- FALSE #set to "TRUE" to return the means and standard deviations of the input parameter choices that are associated with successful (non-failed) runs; only works when loop_parameters is set to "FALSE"; this section of code was not used in Royer et al (2014) (i.e., parameter set to "FALSE")
loop_parameters <- FALSE #set to "TRUE" to test the effect on calculated CO2 and O2 by sequentially varying one input parameter at a time (figures 3-4 in Royer et al 2014); BEWARE: this will take a long time to run (68X longer than a single resampled run)
iteration_threshold <- 10 #maximum number of times the convergence equation for CO2 will iterate before signaling a failed run; in a test with reasonably-well-constrained input parameters (similar to simulations presented in Royer et al, 2014), the number of iterations never exceeded 7; this variable is here mostly as a failsafe stop-gap.


#reading in the two input files
# input <- read.csv("GEOCARB_input_summaries_10%_reservoir_isotope.csv")
input <- read.csv("GEOCARB_input_summaries_50%_sulfate.csv")
scenario = "smO2_5_berner"

time_arrays <- read.csv("GEOCARB_input_arrays_new.csv")
ageN <- length(time_arrays$age) #number of time-steps

#preparing the output files
age <- matrix (time_arrays$age, nrow=ageN, ncol=1); colnames(age) <- "age (Myrs ago)"
failed_runs <- matrix (nrow=ageN, ncol=1); colnames(failed_runs) <- "failed runs (%)" #percent of resampled runs for a given time-step that failed because at least one of the input fluxes was <0, estimated oxygen was <5% or >50% (or <19% or >23% at time t=0), or estimated CO2 was <150 ppm >50000 ppm (or <200 ppm or >300 ppm at time=0)
O2 <- matrix(nrow=ageN, ncol=1); colnames(O2) <- "O2 (%)" #O2 output
CO2 <- matrix(nrow=ageN, ncol=1); colnames(CO2) <- "CO2 (ppm)"  #CO2 output


percentiles_O2 <- matrix(nrow=ageN, ncol=length(percentile_values))
colnames(percentiles_O2) <- paste(percentile_values, '_', 'O2', sep = '')


percentiles_CO2 <- matrix(nrow=ageN, ncol=length(percentile_values))
colnames(percentiles_CO2) <- paste(percentile_values, '_', 'CO2', sep = '')


if (loop_parameters==TRUE)  {
  y <- (4+2*length(percentile_values)) #number of columns needed in summary output file for each parameter loop
  z <- length(input[,"parameter"])  #number of input parameters
  column_header <- matrix(nrow=1,ncol=y*z)  #setting up the column headers for the GEOCARB_output file
  GEOCARB_output <- matrix(nrow=ageN, ncol=y*z)
} else {
  z <- 1}


# !!! Store org burial 
orgb <- matrix(nrow=ageN, ncol=1)
colnames(orgb) <- 'orgb' 
percentiles_orgb <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_orgb) <- paste(percentile_values, '_', 'orgb', sep = '')


# !!! Store carbonate burial 
carb <- matrix(nrow=ageN, ncol=1)
colnames(carb) <- 'carb' 
percentiles_carb <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_carb) <- paste(percentile_values, '_', 'carb', sep = '')


# !!! Store forg
forg <- matrix(nrow=ageN, ncol=1)
colnames(forg) <- 'forg' 
percentiles_forg <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_forg) <- paste(percentile_values, '_', 'forg', sep = '')


# !!! Store pyrite burial 
pyb <- matrix(nrow=ageN, ncol=1)
colnames(pyb) <- 'pyb' 
percentiles_pyb <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_pyb) <- paste(percentile_values, '_', 'pyb', sep = '')



# !!! Store sulfate burial 
sulfb <- matrix(nrow=ageN, ncol=1)
colnames(sulfb) <- 'sulfb' 
percentiles_sulfb <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_sulfb) <- paste(percentile_values, '_', 'sulfb', sep = '')


# !!! Store fpy
fpy <- matrix(nrow=ageN, ncol=1)
colnames(fpy) <- 'fpy' 
percentiles_fpy <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_fpy) <- paste(percentile_values, '_', 'fpy', sep = '')


# !!! Store young org weathering 
orgwy <- matrix(nrow=ageN, ncol=1)
colnames(orgwy) <- 'orgwy' 
percentiles_orgwy <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_orgwy) <- paste(percentile_values, '_', 'orgwy', sep = '')

# !!! Store old org weathering 
orgwa <- matrix(nrow=ageN, ncol=1)
colnames(orgwa) <- 'orgwa' 
percentiles_orgwa <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_orgwa) <- paste(percentile_values, '_', 'orgwa', sep = '')


# !!! Store total org weathering 
orgw_all <- matrix(nrow=ageN, ncol=1)
colnames(orgw_all) <- 'orgw_all' 
percentiles_orgw_all <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_orgw_all) <- paste(percentile_values, '_', 'orgw_all', sep = '')



# !!! Store young carbonate weathering 
carwy <- matrix(nrow=ageN, ncol=1)
colnames(carwy) <- 'carwy' 
percentiles_carwy <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_carwy) <- paste(percentile_values, '_', 'carwy', sep = '')


# !!! Store old carbonate weathering
carwa <- matrix(nrow=ageN, ncol=1)
colnames(carwa) <- 'carwa' 
percentiles_carwa <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_carwa) <- paste(percentile_values, '_', 'carwa', sep = '')


# !!! Store total carbonate weathering 
carw_all <- matrix(nrow=ageN, ncol=1)
colnames(carw_all) <- 'carw_all' 
percentiles_carw_all <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_carw_all) <- paste(percentile_values, '_', 'carw_all', sep = '')


# !!! Store young pyrite weathering 
pywy <- matrix(nrow=ageN, ncol=1)
colnames(pywy) <- 'pywy' 
percentiles_pywy <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_pywy) <- paste(percentile_values, '_', 'pywy', sep = '')


# !!! Store old pyrite weathering 
pywa <- matrix(nrow=ageN, ncol=1)
colnames(pywa) <- 'pywa' 
percentiles_pywa <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_pywa) <- paste(percentile_values, '_', 'pywa', sep = '')

# !!! Store total pyrite weathering 
pyw_all <- matrix(nrow=ageN, ncol=1)
colnames(pyw_all) <- 'pyw_all'
percentiles_pyw_all <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_pyw_all) <- paste(percentile_values, '_', 'pyw_all', sep = '')


# !!! Store young sulfate weathering 
sulfwy <- matrix(nrow=ageN, ncol=1)
colnames(sulfwy) <- 'sulfwy' 
percentiles_sulfwy <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_sulfwy) <- paste(percentile_values, '_', 'sulfwy', sep = '')


# !!! Store old sulfate weathering 
sulfwa <- matrix(nrow=ageN, ncol=1)
colnames(sulfwa) <- 'sulfwa' 
percentiles_sulfwa <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_sulfwa) <- paste(percentile_values, '_', 'sulfwa', sep = '')


# !!! Store total sulfate weathering 
sulfw_all <- matrix(nrow=ageN, ncol=1)
colnames(sulfw_all) <- 'sulfw_all'
percentiles_sulfw_all <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_sulfw_all) <- paste(percentile_values, '_', 'sulfw_all', sep = '')


# !!! Store degassing of org 
morg <- matrix(nrow=ageN, ncol=1)
colnames(morg) <- 'morg' 
percentiles_morg <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_morg) <- paste(percentile_values, '_', 'morg', sep = '')

# !!! Store degassing of carbonate 
mcar <- matrix(nrow=ageN, ncol=1)
colnames(mcar) <- 'mcar' 
percentiles_mcar <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_mcar) <- paste(percentile_values, '_', 'mcar', sep = '')

# !!! Store degassing of pyrite 
mpy <- matrix(nrow=ageN, ncol=1)
colnames(mpy) <- 'mpy' 
percentiles_mpy <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_mpy) <- paste(percentile_values, '_', 'mpy', sep = '')

# !!! Store degassing of sulfate 
msulf <- matrix(nrow=ageN, ncol=1)
colnames(msulf) <- 'msulf' 
percentiles_msulf <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_msulf) <- paste(percentile_values, '_', 'msulf', sep = '')


# !!! Store mass of young org 
orgy <- matrix(nrow=ageN, ncol=1)
colnames(orgy) <- 'orgy' 
percentiles_orgy <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_orgy) <- paste(percentile_values, '_', 'orgy', sep = '')

# !!! Store mass of old org 
orga <- matrix(nrow=ageN, ncol=1)
colnames(orga) <- 'orga' 
percentiles_orga <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_orga) <- paste(percentile_values, '_', 'orga', sep = '')


# !!! Store mass of young carbonate 
cary <- matrix(nrow=ageN, ncol=1)
colnames(cary) <- 'cary' 
percentiles_cary <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_cary) <- paste(percentile_values, '_', 'cary', sep = '')



# !!! Store mass of old carbonate 
cara <- matrix(nrow=ageN, ncol=1)
colnames(cara) <- 'cara' 
percentiles_cara <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_cara) <- paste(percentile_values, '_', 'cara', sep = '')


# !!! Store mass of young py 
pyy <- matrix(nrow=ageN, ncol=1)
colnames(pyy) <- 'pyy' 
percentiles_pyy <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_pyy) <- paste(percentile_values, '_', 'pyy', sep = '')


# !!! Store mass of old py 
pya <- matrix(nrow=ageN, ncol=1)
colnames(pya) <- 'pya' 
percentiles_pya <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_pya) <- paste(percentile_values, '_', 'pya', sep = '')


# !!! Store mass of young sulfate 
sulfy <- matrix(nrow=ageN, ncol=1)
colnames(sulfy) <- 'sulfy' 
percentiles_sulfy <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_sulfy) <- paste(percentile_values, '_', 'sulfy', sep = '')


# !!! Store mass of old sulfate 
sulfa <- matrix(nrow=ageN, ncol=1)
colnames(sulfa) <- 'sulfa' 
percentiles_sulfa <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_sulfa) <- paste(percentile_values, '_', 'sulfa', sep = '')


# !!! Store total C 
TC <- matrix(nrow=ageN, ncol=1)
colnames(TC) <- 'TC' 
percentiles_TC <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_TC) <- paste(percentile_values, '_', 'TC', sep = '')


# !!! Store total S
TS <- matrix(nrow=ageN, ncol=1)
colnames(TS) <- 'TS' 
percentiles_TS <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_TS) <- paste(percentile_values, '_', 'TS', sep = '')


# !!! Store isotope of young org 
dorgy <- matrix(nrow=ageN, ncol=1)
colnames(dorgy) <- 'dorgy' 
percentiles_dorgy <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_dorgy) <- paste(percentile_values, '_', 'dorgy', sep = '')

# !!! Store isotope of old org 
dorga <- matrix(nrow=ageN, ncol=1)
colnames(dorga) <- 'dorga' 
percentiles_dorga <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_dorga) <- paste(percentile_values, '_', 'dorga', sep = '')


# !!! Store isotope of young carbonate 
dcary <- matrix(nrow=ageN, ncol=1)
colnames(dcary) <- 'dcary' 
percentiles_dcary <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_dcary) <- paste(percentile_values, '_', 'dcary', sep = '')


# !!! Store isotope of old carbonate 
dcara <- matrix(nrow=ageN, ncol=1)
colnames(dcara) <- 'dcara' 
percentiles_dcara <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_dcara) <- paste(percentile_values, '_', 'dcara', sep = '')

# !!! Store isotope of young pyrite 
dpyy <- matrix(nrow=ageN, ncol=1)
colnames(dpyy) <- 'dpyy' 
percentiles_dpyy <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_dpyy) <- paste(percentile_values, '_', 'dpyy', sep = '')

# !!! Store isotope of old pyrite 
dpya <- matrix(nrow=ageN, ncol=1)
colnames(dpya) <- 'dpya' 
percentiles_dpya <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_dpya) <- paste(percentile_values, '_', 'dpya', sep = '')


# !!! Store isotope of young sulfate 
dsulfy <- matrix(nrow=ageN, ncol=1)
colnames(dsulfy) <- 'dsulfy' 
percentiles_dsulfy <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_dsulfy) <- paste(percentile_values, '_', 'dsulfy', sep = '')


# !!! Store isotope of old sulfate 
dsulfa <- matrix(nrow=ageN, ncol=1)
colnames(dsulfa) <- 'dsulfa' 
percentiles_dsulfa <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_dsulfa) <- paste(percentile_values, '_', 'dsulfa', sep = '')



# !!! Store isotope fractionation between carbonate and org 
CAPdC <- matrix(nrow=ageN, ncol=1)
colnames(CAPdC) <- 'CAPdC' 
percentiles_CAPdC <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_CAPdC) <- paste(percentile_values, '_', 'CAPdC', sep = '')

# !!! Store isotope fractionation between sulfate and py 
CAPdS <- matrix(nrow=ageN, ncol=1)
colnames(CAPdS) <- 'CAPdS' 
percentiles_CAPdS <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_CAPdS) <- paste(percentile_values, '_', 'CAPdS', sep = '')


# !!! Store total carbon input
input.C <- matrix(nrow=ageN, ncol=1)
colnames(input.C) <- 'input.C' 
percentiles_input.C <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_input.C) <- paste(percentile_values, '_', 'input.C', sep = '')

# !!! Store total carbon input isotope
input.dC <- matrix(nrow=ageN, ncol=1)
colnames(input.dC) <- 'input.dC' 
percentiles_input.dC <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_input.dC) <- paste(percentile_values, '_', 'input.dC', sep = '')


# !!! Store total sulfur input
input.S <- matrix(nrow=ageN, ncol=1)
colnames(input.S) <- 'input.S' 
percentiles_input.S <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_input.S) <- paste(percentile_values, '_', 'input.S', sep = '')

# !!! Store total sulfur input isotope
input.dS <- matrix(nrow=ageN, ncol=1)
colnames(input.dS) <- 'input.dS' 
percentiles_input.dS <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_input.dS) <- paste(percentile_values, '_', 'input.dS', sep = '')



# !!! Store forced total carbonat carbon isotope
dCcar <- matrix(nrow=ageN, ncol=1)
colnames(dCcar) <- 'dCcar' 
percentiles_dCcar <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_dCcar) <- paste(percentile_values, '_', 'dCcar', sep = '')


# !!! Store total carbonate weathering flux to force O2 to follow the line
T_car_weather <- matrix(nrow=ageN, ncol=1)
colnames(T_car_weather) <- 'T_car_weather' 
percentiles_T_car_weather <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_T_car_weather) <- paste(percentile_values, '_', 'T_car_weather', sep = '')


# !!! Store total org weathering flux to force O2 to follow the line
T_org_weather <- matrix(nrow=ageN, ncol=1)
colnames(T_org_weather) <- 'T_org_weather' 
percentiles_T_org_weather <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_T_org_weather) <- paste(percentile_values, '_', 'T_org_weather', sep = '')


# !!! Store total sulfate weathering flux to force O2 to follow the line
T_sulf_weather <- matrix(nrow=ageN, ncol=1)
colnames(T_sulf_weather) <- 'T_sulf_weather' 
percentiles_T_sulf_weather <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_T_sulf_weather) <- paste(percentile_values, '_', 'T_sulf_weather', sep = '')

# !!! Store total pyrite weathering flux to force O2 to follow the line
T_py_weather <- matrix(nrow=ageN, ncol=1)
colnames(T_py_weather) <- 'T_py_weather' 
percentiles_T_py_weather <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_T_py_weather) <- paste(percentile_values, '_', 'T_py_weather', sep = '')



# !!! Store organic carbon weathering scaling factor
soo <- matrix(nrow=ageN, ncol=1)
colnames(soo) <- 'soo' 
percentiles_soo <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_soo) <- paste(percentile_values, '_', 'soo', sep = '')

# !!! Store carbonate carbon weathering scaling factor
scc <- matrix(nrow=ageN, ncol=1)
colnames(scc) <- 'scc' 
percentiles_scc <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_scc) <- paste(percentile_values, '_', 'scc', sep = '')

# !!! Store organic carbon degassing scaling factor
smoo <- matrix(nrow=ageN, ncol=1)
colnames(smoo) <- 'smoo' 
percentiles_smoo <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_smoo) <- paste(percentile_values, '_', 'smoo', sep = '')

# !!! Store carbonate carbon degassing  scaling factor
smcc <- matrix(nrow=ageN, ncol=1)
colnames(smcc) <- 'smcc' 
percentiles_smcc <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_smcc) <- paste(percentile_values, '_', 'smcc', sep = '')

# !!! Store pyrite sulfur weathering scaling factor
spp <- matrix(nrow=ageN, ncol=1)
colnames(spp) <- 'spp' 
percentiles_spp <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_spp) <- paste(percentile_values, '_', 'spp', sep = '')

# !!! Store sulfate sulfur weathering scaling factor
sss <- matrix(nrow=ageN, ncol=1)
colnames(sss) <- 'sss' 
percentiles_sss <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_sss) <- paste(percentile_values, '_', 'sss', sep = '')


# !!! Store pyrite sulfur degassing scaling factor
smpp <- matrix(nrow=ageN, ncol=1)
colnames(smpp) <- 'smpp' 
percentiles_smpp <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_smpp) <- paste(percentile_values, '_', 'smpp', sep = '')

# !!! Store sulfate sulfur degassing scaling factor
smss <- matrix(nrow=ageN, ncol=1)
colnames(smss) <- 'smss' 
percentiles_smss <- matrix(nrow=ageN, ncol=length(percentile_values));
colnames(percentiles_smss) <- paste(percentile_values, '_', 'smss', sep = '')


######################################################
### code for generating and filling input matrices ###
######################################################

#########
#function for generating resamples for time-dependent arrays, including the potential clipping of distributions
resamples_array <- function(x, name, row_position) {
  x <- matrix(nrow=ageN, ncol=resampleN)  #create empty matrix
  col_num <- which(colnames(time_arrays)==name) #matches a column number in the input arrays file to the array in question
  
  #if resampling=1 or if resampling is turned off for the time-dependent array in question (a full matrix is generated but with resampleN mean values at each time-step)
  if (resampleN==1 || input[row_position, "resample"]==FALSE) {
    x <- matrix(time_arrays[, col_num], nrow=ageN, ncol=resampleN)  }
  
  #resampling following a normal distribution
  if (input[row_position, "resample"]==TRUE & input[row_position, "distribution_type"]=="gaussian" & resampleN>1) {
    for (i in 1:ageN) {
      temp_resample <- NULL
      temp_resample <- rnorm(n=resampleN, mean=time_arrays[i,col_num], sd=time_arrays[i,col_num+1]/2)
      temp_resample <- replace(temp_resample,temp_resample<=input[row_position,"lower_limit"],input[row_position,"lower_limit"]+0.0001)
      temp_resample <- replace(temp_resample,temp_resample>=input[row_position,"upper_limit"],input[row_position,"upper_limit"]-0.0001)   
      x[i,] <- matrix(temp_resample, nrow=1, ncol=resampleN)
    } #end of time-sequence loop (i)      
  }    
  
  #resampling following a lognormal distribution
  if (input[row_position, "resample"]==TRUE & input[row_position, "distribution_type"]=="lognormal" & resampleN>1) {
    for (i in 1:ageN) {
      temp_resample <- NULL
      temp_resample <- rlnorm(n=resampleN, meanlog=log(time_arrays[i,col_num]), sdlog=log(time_arrays[i,col_num+1])/2)
      temp_resample <- replace(temp_resample,temp_resample<=input[row_position,"lower_limit"],input[row_position,"lower_limit"]+0.0001)
      temp_resample <- replace(temp_resample,temp_resample>=input[row_position,"upper_limit"],input[row_position,"upper_limit"]-0.0001)   
      x[i,] <- matrix(temp_resample, nrow=1, ncol=resampleN)
    } #end of time-sequence loop (i) 
  }
  x <- x    
} #end of function "resamples_array"

#########

#function for generating resamples for constant parameters, including the potential clipping of distributions
resamples_constant <- function(x, row_position) {
  
  #if resampling=1 or if resampling is turned off for the constant parameter in question (a full matrix is filled with the mean value)
  if (resampleN==1 || input[row_position, "resample"]==FALSE)  {
    x <- matrix(input[row_position, "mean"], nrow=ageN, ncol=resampleN)  }
  
  #resampling following a normal distribution
  if (input[row_position, "resample"]==TRUE & input[row_position, "distribution_type"]=="gaussian" & resampleN>1) {
    temp_resample <- NULL
    temp_resample <- rnorm(n=resampleN, mean=input[row_position,"mean"], sd=input[row_position,"two_sigma"]/2)
    temp_resample <- replace(temp_resample,temp_resample<=input[row_position,"lower_limit"],input[row_position,"lower_limit"]+0.0001)
    temp_resample <- replace(temp_resample,temp_resample>=input[row_position,"upper_limit"],input[row_position,"upper_limit"]-0.0001)
    x <- matrix(temp_resample, nrow=ageN, ncol=resampleN, byrow=TRUE)
  }
  
  #resampling following a lognormal distribution
  if (input[row_position, "resample"]==TRUE & input[row_position, "distribution_type"]=="lognormal" & resampleN>1) {
    temp_resample <- NULL
    temp_resample <- rlnorm(n=resampleN, meanlog=log(input[row_position,"mean"]), sdlog=log(input[row_position,"two_sigma"])/2)
    temp_resample <- replace(temp_resample,temp_resample<=input[row_position,"lower_limit"],input[row_position,"lower_limit"]+0.0001)
    temp_resample <- replace(temp_resample,temp_resample>=input[row_position,"upper_limit"],input[row_position,"upper_limit"]-0.0001)
    x <- matrix(temp_resample, nrow=ageN, ncol=resampleN, byrow=TRUE)
  }    
  x <- x
} #end of function "resamples_constant"

#y=young; a=old; p=pyrite; s=sulfate; c=carbonate; si=silicates; g=organic matter; b=burial; m=degassing; w=weathering
#masses are in units of 10^18 mol
#fluxes ("F" prefix) are in units of 10^18 mol Myrs-1
#rates ("k" prefix) are in units of Myrs-1
#stable isotopic compositions ("d" prefix) are in per mil units

#generate random distributions for the time-dependent parameters using the function "resamples_array"
Sr <- resamples_array(Sr, "Sr", row_position=which(input[,"parameter"]=="Sr"))  #87Sr/86Sr of shallow-marine carbonate ([87Sr/86Sr - 0.7] x 10^4)
d34S <- resamples_array(d34S, "d34S", row_position=which(input[,"parameter"]=="d34S"))  #d34S of marine sulfate sulfur (per mil) (from Wu et al, 2010); called "DLSOC" in BASIC code
fR <- resamples_array(fR, "fR", row_position=which(input[,"parameter"]=="fR"))  #effect of relief on chemical weathering at time (t) relative to the present-day; calculated from equation (5) in Berner (2006b)
fL <- resamples_array(fL, "fL", row_position=which(input[,"parameter"]=="fL"))  #land area covered by carbonates at time (t) relative to the present-day
if (Godderis==TRUE) {
  fA <- resamples_array(fA, "fA_Godderis", row_position=which(input[,"parameter"]=="fA")) #land area at time (t) relative to the present-day
  fAw_fA <- resamples_array(fAw_fA, "fAw_fA_Godderis", row_position=which(input[,"parameter"]=="fAw_fA"))  #fraction of land area experiencing chemical weathering (runoff > 0)
  fD <- resamples_array(fD, "fD_Godderis", row_position=which(input[,"parameter"]=="fD")) #change in global river runoff at time (t) relative to the present-day in the absence of changes in solar luminosity and CO2 (i.e., mainly due to changes in paleogeography)
  GEOG <- resamples_array(GEOG, "GEOG_Godderis", row_position=which(input[,"parameter"]=="GEOG")) #change in land mean surface temperature for areas experiencing chemical weathering (runoff > 0) at time (t) relative to the present-day in the absence of changes in solar luminosity and CO2 (i.e., mainly due to changes in paleogeography) (K)
} else {
  fA <- resamples_array(fA, "fA", row_position=which(input[,"parameter"]=="fA"))
  fAw_fA <- resamples_array(fAw_fA, "fAw_fA", row_position=which(input[,"parameter"]=="fAw_fA"))
  fD <- resamples_array(fD, "fD", row_position=which(input[,"parameter"]=="fD"))
  GEOG <- resamples_array(GEOG, "GEOG", row_position=which(input[,"parameter"]=="GEOG"))
} #end of if...else loop for G
RT <- resamples_array(RT, "RT", row_position=which(input[,"parameter"]=="RT"))  #coefficient relating continental runoff to temperature change (runoff/runoff(0)=1+RT*(T-T0)), as determined from the GCM simulations of Gois called "Y" in Berner 2004 and "RUN" in GEOCARB III; in the BASIC scripts, RT is assigned a value of 0.045 during times with large ice sheets and a value of 0.025 for all other times; there is little difference in estimated CO2 between these two approaches for RT
fSR <- resamples_array(fSR, "fSR", row_position=which(input[,"parameter"]=="fSR"))  #seafloor creation rate at time (t) relative to the present-day
fC <- resamples_array(fC, "fC", row_position=which(input[,"parameter"]=="fC"))  #effect of carbonate content of subducting oceanic crust on CO2 degassing rate at time (t) relative to the present-day

# !!! add S fractionation
sCAPd34S <- resamples_array(sCAPd34S, "sCAPd34S", row_position=which(input[,"parameter"]=="sCAPd34S"))  #effect of carbonate content of subducting oceanic crust on CO2 degassing rate at time (t) relative to the present-day

# !!! add C fractionation
sCAPd13C <- resamples_array(sCAPd13C, "sCAPd13C", row_position=which(input[,"parameter"]=="sCAPd13C"))  #effect of carbonate content of subducting oceanic crust on CO2 degassing rate at time (t) relative to the present-day


# # !!! update d13C  
# d13C <- resamples_array(d13C, "d13C_berner", row_position=which(input[,"parameter"]=="d13C"))  #d13C of shallow-marine carbonate (per mil); called "DLCOC" in BASIC code

# !!! update d13C  
d13C <- resamples_array(d13C, "d13C_10_smooth", row_position=which(input[,"parameter"]=="d13C"))  #d13C of shallow-marine carbonate (per mil); called "DLCOC" in BASIC code


# !!! add smO2 and need to update the parameter name (like "smo2_5_mills") for different scenarios
smO2 <- resamples_array(smO2, scenario, row_position=which(input[,"parameter"]=="smO2"))  #d13C of shallow-marine carbonate (per mil); called "DLCOC" in BASIC code


#generate random distributions for the constant parameters using the function "resamples_constant"
ACT <- resamples_constant(ACT, row_position=which(input[,"parameter"]=="ACT"))  #activation energy (E) for dissolution of Ca- and Mg-silicates on land, where ACT = E/RTT0 (1/K); called "Z" in Berner (2004)
ACTcarb <- resamples_constant(ACTcarb, row_position=which(input[,"parameter"]=="ACTcarb"))  #activation energy (E) for dissolution of carbonates on land (see p. 53 in Berner 2004)
VNV <- resamples_constant(VNV, row_position=which(input[,"parameter"]=="VNV"))  #rate ratio of chemical weathering in volcanic to non-volcanic silicate rocks (called "Wv/Wnv in Berner 2006b, and "basalt/granite" in Berner 2008)
NV <- resamples_constant(NV, row_position=which(input[,"parameter"]=="NV")) #coefficient relating physical erosion to the mean 87Sr/86Sr of non-volcanic silicate rocks
exp_NV <- resamples_constant(exp_NV, row_position=which(input[,"parameter"]=="exp_NV")) #exponent relating physical erosion to the mean 87Sr/86Sr of non-volcanic silicate rocks (see Berner 2008)
LIFE <- resamples_constant(LIFE, row_position=which(input[,"parameter"]=="LIFE")) #rate ratio of chemical weathering in a minimally-vegetated to present-day (angiosperm dominated) world
GYM <- resamples_constant(GYM, row_position=which(input[,"parameter"]=="GYM"))  #rate ratio of chemical weathering by gymnosperms to angiosperms
FERT <- resamples_constant(FERT, row_position=which(input[,"parameter"]=="FERT")) #exponent reflecting the fraction of vegetation whose growth is stimulated by elevated CO2; FERT is related to enhanced chemical weathering by the Michaelis-Menton expression [2RCO2/(1+RCO2)]^FERT, which is called fBb(CO2) in Berner 2004 (p. 24)
exp_fnBb <- resamples_constant(exp_fnBb, row_position=which(input[,"parameter"]=="exp_fnBb")) #exponent used to describe the effect of climate on silicate or carbonate weathering in the absence of vascular plants at time (t) relative to the present-day (see pp. 25 & 53-54 in Berner 2004 and pp. 67-68 in Berner 1994)
deltaT2X <- resamples_constant(deltaT2X, row_position=which(input[,"parameter"]=="deltaT2X")) #climate sensitivity (K per CO2 doubling)
GLAC <- resamples_constant(GLAC, row_position=which(input[,"parameter"]=="GLAC")) #factor by which deltaT2X changes during times with large continental ice sheets
J <- resamples_constant(J, row_position=which(input[,"parameter"]=="J"))  #coefficient used to calculate CAPd13C (called "alphac" in BASIC code), the stable carbon isotopic fractionation between shallow-marine carbonate and shallow-marine organic matter
n <- resamples_constant(n, row_position=which(input[,"parameter"]=="n"))  #coefficient used to calculate CAPd34S (called "alphas" in BASIC code), the stable sulfur isotopic fractionation between marine sulfate sulfur and marine pyrite sulfur
Ws <- resamples_constant(Ws, row_position=which(input[,"parameter"]=="Ws")) #effect on temperature from the linear increase in solar luminosity over time (K per 570 Myrs)
exp_fD <- resamples_constant(exp_fD, row_position=which(input[,"parameter"]=="exp_fD")) #exponent that scales the dilution of dissolved HCO3- with runoff (fD) (see pp. 29-31 & 34-36 in Berner 2004)
#GEOCARBSULF is run forward in time; the following are present-day (t = 0 Myrs ago) or initial values (t = 570 Myrs ago) of parameters that are subsequently recalculated at each time-step (see Berner 2004, 2006a for details)
Fwpa_0 <- resamples_constant(Fwpa_0, row_position=which(input[,"parameter"]=="Fwpa_0")) #sulfate flux from oxidative weathering of old pyrite at present-day
Fwsa_0 <- resamples_constant(Fwsa_0, row_position=which(input[,"parameter"]=="Fwsa_0")) #sulfate flux from weathering of CaSO4 sulfur at present-day
Fwga_0 <- resamples_constant(Fwga_0, row_position=which(input[,"parameter"]=="Fwga_0")) #carbon flux from weathering of old sedimentary organic matter at present-day
Fwca_0 <- resamples_constant(Fwca_0, row_position=which(input[,"parameter"]=="Fwca_0")) #carbon flux from weathering of old Ca and Mg carbonates at present-day
Fmg_0 <- resamples_constant(Fmg_0, row_position=which(input[,"parameter"]=="Fmg_0"))  #carbon degassing flux from volcanism, metamorphism, and diagenesis of organic matter at present-day
Fmc_0 <- resamples_constant(Fmc_0, row_position=which(input[,"parameter"]=="Fmc_0"))  #carbon degassing flux from volcanism, metamorphism, and diagenesis of carbonates at present-day
Fmp_0 <- resamples_constant(Fmp_0, row_position=which(input[,"parameter"]=="Fmp_0"))  #sulfur degassing flux from volcanism, metamorphism, and diagenesis of pyrite at present-day
Fms_0 <- resamples_constant(Fms_0, row_position=which(input[,"parameter"]=="Fms_0"))  #sulfur degassing flux from volcanism, metamorphism, and diagenesis of CaSO4 sulfur at present-day
Fwsi_0 <- resamples_constant(Fwsi_0, row_position=which(input[,"parameter"]=="Fwsi_0")) #weathering flux for all Ca and Mg silicates at present-day
Xvolc_0 <- resamples_constant(Xvolc_0, row_position=which(input[,"parameter"]=="Xvolc_0"))  #fraction of total Ca and Mg silicate weathering drived from volcanic rocks at present-day
CAPd13C_0 <- resamples_constant(CAPd13C_0, row_position=which(input[,"parameter"]=="CAPd13C_0"))  #stable carbon isotopic fractionation between shallow-marine carbonate and shallow-marine organic matter at present-day
CAPd34S_0 <- resamples_constant(CAPd34S_0, row_position=which(input[,"parameter"]=="CAPd34S_0"))  #stable sulfur isotopic fractionation between shallow-marine CaSO4 sulfur and pyrite sulfur at present-day
oxygen_570 <- resamples_constant(oxygen_570, row_position=which(input[,"parameter"]=="oxygen_570")) #mass of atmospheric O2 at 570 Myrs ago
Gy_570 <- resamples_constant(Gy_570, row_position=which(input[,"parameter"]=="Gy_570")) #mass of young crustal organic carbon at 570 Myrs ago
Cy_570 <- resamples_constant(Cy_570, row_position=which(input[,"parameter"]=="Cy_570")) #mass of young crustal carbonate carbon at 570 Myrs ago
Ca_570 <- resamples_constant(Ca_570, row_position=which(input[,"parameter"]=="Ca_570")) #mass of old crustal carbonate carbon at 570 Myrs ago
Ssy_570 <- resamples_constant(Ssy_570, row_position=which(input[,"parameter"]=="Ssy_570"))  #mass of young CaSO4 sulfur at 570 Myrs ago
Spy_570 <- resamples_constant(Spy_570, row_position=which(input[,"parameter"]=="Spy_570"))  #mass of young pyrite sulfur at 570 Myrs ago
dlsy_570 <- resamples_constant(dlsy_570, row_position=which(input[,"parameter"]=="dlsy_570")) #d34S of young CaSO4 sulfur at 570 Myrs ago
dlcy_570 <- resamples_constant(dlcy_570, row_position=which(input[,"parameter"]=="dlcy_570")) #d13C of young carbonate carbon at 570 Myrs ago
dlpy_570 <- resamples_constant(dlpy_570, row_position=which(input[,"parameter"]=="dlpy_570")) #d34S of young pyrite sulfur at 570 Myrs ago
dlpa_570 <- resamples_constant(dlpa_570, row_position=which(input[,"parameter"]=="dlpa_570")) #d34S of old pyrite sulfur at 570 Myrs ago
dlgy_570 <- resamples_constant(dlgy_570, row_position=which(input[,"parameter"]=="dlgy_570")) #d13C of young organic matter at 570 Myrs ago
dlga_570 <- resamples_constant(dlga_570, row_position=which(input[,"parameter"]=="dlga_570")) #d13C of old organic matter at 570 Myrs ago
Rcy_570 <- resamples_constant(Rcy_570, row_position=which(input[,"parameter"]=="Rcy_570"))  #87Sr/86Sr of young carbonates undergoing weathering at 570 Myrs ago
Rca_570 <- resamples_constant(Rca_570, row_position=which(input[,"parameter"]=="Rca_570"))  #87Sr/86Sr of old carbonates undergoing weathering at 570 Myrs ago
Rv_570 <- resamples_constant(Rv_570, row_position=which(input[,"parameter"]=="Rv_570")) #87Sr/86Sr of sub-aerial and submarine volcanic rocks at 570 Myrs ago
Rg_570 <- resamples_constant(Rg_570, row_position=which(input[,"parameter"]=="Rg_570")) #87Sr/86Sr of non-volcanic silicates (granites) at 570 Myrs ago
Fob <- resamples_constant(Fob, row_position=which(input[,"parameter"]=="Fob"))  #Ca and Mg flux between basalt and seawater
COC <- resamples_constant(COC, row_position=which(input[,"parameter"]=="COC"))  #mass of carbon in ocean
Ga <- resamples_constant(Ga, row_position=which(input[,"parameter"]=="Ga")) #mass of old crustal organic carbon
Ssa <- resamples_constant(Ssa, row_position=which(input[,"parameter"]=="Ssa"))  #mass of old CaSO4 sulfur
Spa <- resamples_constant(Spa, row_position=which(input[,"parameter"]=="Spa"))  #mass of old pyrite sulfur
ST <- resamples_constant(ST, row_position=which(input[,"parameter"]=="ST")) #mass of sulfur in oceans + "interacting rocks" (i.e., carbon in rocks undergoing weathering, burial, etc.)
dlst <- resamples_constant(dlst, row_position=which(input[,"parameter"]=="dlst")) #d34S of ST
CT <- resamples_constant(CT, row_position=which(input[,"parameter"]=="CT")) #mass of carbon in oceans + "interacting rocks" (i.e., carbon in rocks undergoing weathering, burial, etc.)
dlct <- resamples_constant(dlct, row_position=which(input[,"parameter"]=="dlct")) #d13C of CT
kwpy <- resamples_constant(kwpy, row_position=which(input[,"parameter"]=="kwpy")) #rate constant expressing mass dependence for young pyrite sulfur
kwsy <- resamples_constant(kwsy, row_position=which(input[,"parameter"]=="kwsy")) #rate constant expressing mass dependence for young CaSO4 sulfur
kwgy <- resamples_constant(kwgy, row_position=which(input[,"parameter"]=="kwgy")) #rate constant expressing mass dependence for young organic matter weathering
kwcy <- resamples_constant(kwcy, row_position=which(input[,"parameter"]=="kwcy")) #rate constant expressing mass dependence for young carbonate weathering


#generate empty matrices for individual calculations of O2 and CO2
O2_resamples <- matrix(nrow=ageN, ncol=resampleN)
CO2_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for orgb
orgb_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for carb
carb_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for forg
forg_resamples <- matrix(nrow=ageN, ncol=resampleN)


# !!! generate empty matrices for pyb
pyb_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for sulfb
sulfb_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for fpy
fpy_resamples <- matrix(nrow=ageN, ncol=resampleN) 

# !!! generate empty matrices for young org weathering
orgwy_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for old org weathering
orgwa_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for young car weathering
carwy_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for old car weathering
carwa_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for young py weathering
pywy_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for old py weathering
pywa_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for young sulfate weathering
sulfwy_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for old sulfate weathering
sulfwa_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for org degassing
morg_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for car degassing
mcar_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for py degassing
mpy_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for sulfate degassing
msulf_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for mass of young org
orgy_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for mass of old org
orga_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for mass of young car
cary_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for mass of old car
cara_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for mass of young py
pyy_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for mass of old py
pya_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for mass of young sulfate
sulfy_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for mass of old sulfate
sulfa_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for total C
TC_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for total S
TS_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for isotope of young org
dorgy_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for isotope of old org
dorga_resamples <- matrix(nrow=ageN, ncol=resampleN)  

# !!! generate empty matrices for isotope of young car
dcary_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for isotope of old car
dcara_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for isotope of young py
dpyy_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for isotope of old py
dpya_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for isotope of young sulfate
dsulfy_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for isotope of old sulfate
dsulfa_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for isotope fractionation between carbonate and org
CAPdC_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for isotope fractionation between sulfate and py
CAPdS_resamples <- matrix(nrow=ageN, ncol=resampleN) 

# !!! generate empty matrices for total carbon input
input.C_resamples <- matrix(nrow=ageN, ncol=resampleN) 

# !!! generate empty matrices for total carbon input isotope
input.dC_resamples <- matrix(nrow=ageN, ncol=resampleN) 


# !!! generate empty matrices for total sulfur input
input.S_resamples <- matrix(nrow=ageN, ncol=resampleN) 

# !!! generate empty matrices for total sulfur input isotope
input.dS_resamples <- matrix(nrow=ageN, ncol=resampleN) 


# !!! generate empty matrices for carbonate isotope
dCcar_resamples <- matrix(nrow=ageN, ncol=resampleN) 

# !!! generate empty matrices for total organic carbon weathering
T_org_weather_resamples <- matrix(nrow=ageN, ncol=resampleN)


# !!! generate empty matrices for total carbonate carbon weathering
T_car_weather_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for total sulfate sulfur weathering
T_sulf_weather_resamples <- matrix(nrow=ageN, ncol=resampleN) 

# !!! generate empty matrices for total pyrite sulfur weathering
T_py_weather_resamples <- matrix(nrow=ageN, ncol=resampleN) 

# !!! generate empty matrices for organic carbon weathering scaling factor
soo_resamples <- matrix(nrow=ageN, ncol=resampleN) 

# !!! generate empty matrices for carbonate carbon weathering scaling factor
scc_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for organic carbon degassing scaling factor
smoo_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for carbonate carbon degassing scaling factor
smcc_resamples <- matrix(nrow=ageN, ncol=resampleN)

# !!! generate empty matrices for pyrite sulfur weathering scaling factor
spp_resamples <- matrix(nrow=ageN, ncol=resampleN) 

# !!! generate empty matrices for sulfate sulfur weathering scaling factor
sss_resamples <- matrix(nrow=ageN, ncol=resampleN) 

# !!! generate empty matrices for pyrite sulfur degassing scaling factor
smpp_resamples <- matrix(nrow=ageN, ncol=resampleN) 

# !!! generate empty matrices for sulfate sulfur degassing scaling factor
smss_resamples <- matrix(nrow=ageN, ncol=resampleN) 


print("Done filling input arrays")


######################################
### code for GEOCARBSULFvolc model ###
######################################

Dt <- 10  #time-step (millions of years, Myrs)
oxygen_0 <- 38 #mass of atmospheric O2 at the present-day

#start of resampling loop (i)
for (i in 1:resampleN) {
  RCO2 <- 10  #atmospheric CO2 (ratio between CO2 at time t to the Pleistocene mean [taken as 250 ppm]); this initial value is a place-holder (it is solved for explicitly)
  #these variables are recalculated at each time-step
  oxygen <- oxygen_570[1,i]
  Gy <- Gy_570[1,i]
  Cy <- Cy_570[1,i]
  Ca <- Ca_570[1,i]
  Ssy <- Ssy_570[1,i]
  Spy <- Spy_570[1,i]
  dlsy <- dlsy_570[1,i]
  dlcy <- dlcy_570[1,i]
  dlpy <- dlpy_570[1,i]
  dlpa <- dlpa_570[1,i]
  dlgy <- dlgy_570[1,i]
  dlga <- dlga_570[1,i]
  dlsa <- (dlst[1,i]*ST[1,i]-(dlpy*Spy+dlsy*Ssy+dlpa*Spa[1,i]))/Ssa[1,i]  #d34S value of old CaSO4 sulfur
  dlca <- (dlct[1,i]*CT[1,i]-(dlgy*Gy+dlcy*Cy+dlga*Ga[1,i]))/Ca  #d13C value of old crustal carbonate carbon   
  Rcy <- Rcy_570[1,i]
  Rca <- Rca_570[1,i] 
  Rcy <- Rcy_570[1,i]
  Rca <- Rca_570[1,i]
  
  #start the nested time loop (j)
  for (j in 1:ageN) {
    failed_run <- FALSE  #flag for failed runs
    t <- time_arrays[j,"age"] #age of time-step, in Myrs ago
    
    #calculate factors that are influenced by glacial vs. non-glacial state ("GCM")
    if ((t<=330 & t>=260) || t<=40) { #for glacial periods between 260 and 330 Myrs ago and between 35 and 0 Myrs ago; BASIC code calls for a 270-340 Myrs ago interval, but see Fielding et al 2008 (GSA Special Paper 441: 343-354) for justification
      GCM <- GLAC[j,i]*deltaT2X[j,i]/log(2) #in GEOCARBSULF, deltaT2X = GCM*ln(2); called capital gamma in Berner (2004)
    } else {
      GCM <- deltaT2X[j,i]/log(2) }
    
    #calculate factors that are influenced by glacial vs. non-glacial state ("GCM")
    if ((t<=330 & t>=260) || t<=40) { #for glacial periods between 260 and 330 Myrs ago and between 35 and 0 Myrs ago; BASIC code calls for a 270-340 Myrs ago interval, but see Fielding et al 2008 (GSA Special Paper 441: 343-354) for justification
      GCM <- GLAC[j,i]*deltaT2X[j,i]/log(2) #in GEOCARBSULF, deltaT2X = GCM*ln(2); called capital gamma in Berner (2004)
    } else {
      GCM <- deltaT2X[j,i]/log(2) }
    
    #calculate factors related to vegetation type
    if (t<=570 & t>380) { #vegetation = domination by non-vascular plants
      fE <- LIFE[j,i] #effect of plants on weathering rate at time (t) to the present-day
      fBB <- (1+ACTcarb[j,i]*GCM*log(RCO2)-ACTcarb[j,i]*Ws[j,i]*(t/570)+ACTcarb[j,i]*GEOG[j,i])*RCO2^exp_fnBb[j,i]  }  #effect of CO2 on plant-assisted weathering for carbonates at time (t) to the present-day
    if (t<=380 & t>350) { #vegetation = ramp-up to gymnosperm domination
      fE <- (GYM[j,i]-LIFE[j,i])*((380-t)/30)+LIFE[j,i]
      fBB <- ((380-t)/30)*((1+ACTcarb[j,i]*GCM*log(RCO2)-ACTcarb[j,i]*Ws[j,i]*(t/570)+ACTcarb[j,i]*GEOG[j,i])*(2*RCO2/(1+RCO2))^FERT[j,i])+((t-350)/30)*(1+ACTcarb[j,i]*GCM*log(RCO2)-ACTcarb[j,i]*Ws[j,i]*(t/570)+ACTcarb[j,i]*GEOG[j,i])*RCO2^exp_fnBb[j,i] }
    if (t<=350)  {
      fBB <- (1+ACTcarb[j,i]*GCM*log(RCO2)-ACTcarb[j,i]*Ws[j,i]*(t/570)+ACTcarb[j,i]*GEOG[j,i])*(2*RCO2/(1+RCO2))^FERT[j,i] }
    if (t<=350 & t>130)  { #vegetation = gymnosperm domination
      fE <- GYM[j,i]  }
    if (t<=130 & t>80) { #vegetation = ramp-up to angiosperm domination
      fE <- (1-GYM[j,i])*((130-t)/50)+GYM[j,i]  }
    if (t<=80)  { #vegetation = angiosperm domination
      fE <- 1 }
    
    
    #calculate source fluxes (weathering and degassing); see pp. 5660-5661 in Berner (2006a) for a discussion of the "Fwp", "Fws", and "Fwg" fluxes
    
    fwr_now = 1 # This is the scaling factor that could influence organic carbon weathering
    
    
    Fwgy_ini <- fA[j,i]*fR[j,i]*fwr_now*kwgy[j,i]*Gy
    Fwga_ini <- fR[j,i]*fwr_now*Fwga_0[j,i]
    Fwcy_ini <- fA[j,i]*fD[j,i]*fL[j,i]*fE*fBB*kwcy[j,i]*Cy
    Fwca_ini <- fA[j,i]*fD[j,i]*fL[j,i]*fE*fBB*Fwca_0[j,i]
    Fmg_ini <- fSR[j,i]*Fmg_0[j,i]
    Fmc_ini <- fSR[j,i]*fC[j,i]*Fmc_0[j,i]
    
    Fwpy_ini <- fA[j,i]*fR[j,i]*kwpy[j,i]*Spy
    Fwpa_ini <- fR[j,i]*Fwpa_0[j,i]
    Fwsy_ini <- fA[j,i]*fD[j,i]*kwsy[j,i]*Ssy
    Fwsa_ini <- fA[j,i]*fD[j,i]*Fwsa_0[j,i]
    Fmp_ini <- fSR[j,i]*Fmp_0[j,i]
    Fms_ini <- fSR[j,i]*Fms_0[j,i]
    
    
    # CAPd34S <- CAPd34S_0[j,i]*((oxygen/oxygen_0)^n[j,i]) #isotopic fractionation between sulfate sulfur and pyrite sulfur (see Berner 2006a)
    
    CAPd34S <- sCAPd34S[j,i] # CAPd34S comes from records
    
    # CAPd13C <- sCAPd13C[j,i]
    
    CAPd13C <- CAPd13C_0[j,i]+J[j,i]*(oxygen/oxygen_0-1) #isotopic fractionation between carbonate and organic matter (see Berner 2006a)
    
    
    Fbp_ini <- (1/CAPd34S)*((d34S[j,i]-dlsy)*Fwsy_ini+(d34S[j,i]-dlsa)*Fwsa_ini+(d34S[j,i]-dlpy)*Fwpy_ini+(d34S[j,i]-dlpa)*Fwpa_ini+(d34S[j,i]-dlsa)*Fms_ini+(d34S[j,i]-dlpa)*Fmp_ini) #burial flux of pyrite
    
    
    
    s_dCcar <- d13C[j,i]
    
    Fbg_ini <- (1/CAPd13C)*((d13C[j,i]-dlcy)*Fwcy_ini+(d13C[j,i]-dlca)*Fwca_ini+(d13C[j,i]-dlgy)*Fwgy_ini+(d13C[j,i]-dlga)*Fwga_ini+(d13C[j,i]-dlca)*Fmc_ini+(d13C[j,i]-dlga)*Fmg_ini) #burial flux of organic carbon
    
    
    func_min <- function(par) {
      
      
      par_Fwgy <- par[1]*Fwgy_ini #weathering flux of young organic carbon
      par_Fwga <- par[1]*Fwga_ini  #weathering flux of old organic carbon
      
      par_Fwcy <- par[2]*Fwcy_ini
      par_Fwca <- par[2]*Fwca_ini
      
      par_Fmg <- par[3]*Fmg_ini
      par_Fmc <- par[4]*Fmc_ini
      
      par_Fwpy <- par[5]*Fwpy_ini
      par_Fwpa <- par[5]*Fwpa_ini
      
      par_Fwsy <- par[6]*Fwsy_ini
      par_Fwsa <- par[6]*Fwsa_ini
      
      
      par_Fmp <- par[7]*Fmp_ini
      
      
      par_Fbg <- (1/CAPd13C)*((d13C[j,i]-dlcy)*par_Fwcy+(d13C[j,i]-dlca)*par_Fwca+(d13C[j,i]-dlgy)*par_Fwgy+(d13C[j,i]-dlga)*par_Fwga+(d13C[j,i]-dlca)*par_Fmc+(d13C[j,i]-dlga)*par_Fmg) #burial flux of organic carbon
      
      
      const_Fbp <- (1/CAPd34S)*((d34S[j,i]-dlsy)*par_Fwsy+(d34S[j,i]-dlsa)*par_Fwsa+(d34S[j,i]-dlpy)*par_Fwpy+(d34S[j,i]-dlpa)*par_Fwpa+(d34S[j,i]-dlpa)*par_Fmp) 
      
      
      par_Fms <- (8/15*(smO2[j,i]/Dt - par_Fbg + (par_Fwgy+par_Fwga+par_Fmg) + (15/8)*(par_Fwpy+par_Fwpa+par_Fmp)) - const_Fbp)/((1/CAPd34S)*(d34S[j,i]-dlsa))
      
      
      
      abs((par_Fwgy + par_Fwga - Fwgy_ini - Fwga_ini)/(Fwgy_ini + Fwga_ini)) + abs((par_Fwcy + par_Fwca - Fwcy_ini - Fwca_ini)/(Fwcy_ini + Fwca_ini)) + 
        abs((par_Fmg-Fmg_ini)/Fmg_ini) + abs((par_Fmc-Fmc_ini)/(Fmc_ini)) + abs((par_Fwpy + par_Fwpa - Fwpy_ini - Fwpa_ini)/(Fwpy_ini + Fwpa_ini)) +
        abs((par_Fwsy + par_Fwsa - Fwsy_ini - Fwsa_ini)/(Fwsy_ini + Fwsa_ini)) + abs((par_Fmp-Fmp_ini)/Fmp_ini) + abs((par_Fms-Fms_ini)/(Fms_ini))
      
    }
    
    
    result <- optim(par = c(1, 1, 1, 1, 1, 1, 1), func_min, lower=c(0,0,0,0,0,0,0), method="L-BFGS-B")
    
    oo <- result$par[1]
    cc <- result$par[2]
    moo <- result$par[3]
    mcc <- result$par[4]
    pp <- result$par[5]
    ss <- result$par[6]
    mpp <- result$par[7]
    
    Fwgy <- oo*Fwgy_ini #weathering flux of young organic carbon
    Fwga <- oo*Fwga_ini  #weathering flux of old organic carbon
    Fwcy <- cc*Fwcy_ini
    Fwca <- cc*Fwca_ini
    Fmg <- moo*Fmg_ini
    Fmc <- mcc*Fmc_ini
    
    Fwpy <- pp*Fwpy_ini
    Fwpa <- pp*Fwpa_ini
    
    Fwsy <- ss*Fwsy_ini
    Fwsa <- ss*Fwsa_ini
    
    Fmp <- mpp*Fmp_ini
    
    
    Fbg <- (1/CAPd13C)*((d13C[j,i]-dlcy)*Fwcy+(d13C[j,i]-dlca)*Fwca+(d13C[j,i]-dlgy)*Fwgy+(d13C[j,i]-dlga)*Fwga+(d13C[j,i]-dlca)*Fmc+(d13C[j,i]-dlga)*Fmg) #burial flux of organic carbon
    
    
    const_Fbp <- (1/CAPd34S)*((d34S[j,i]-dlsy)*Fwsy+(d34S[j,i]-dlsa)*Fwsa+(d34S[j,i]-dlpy)*Fwpy+(d34S[j,i]-dlpa)*Fwpa+(d34S[j,i]-dlpa)*Fmp) 
    
    Fms_old <- Fms_ini
    
    Fms <- (8/15*(smO2[j,i]/Dt - Fbg + (Fwgy+Fwga+Fmg) + (15/8)*(Fwpy+Fwpa+Fmp)) - const_Fbp)/((1/CAPd34S)*(d34S[j,i]-dlsa))
    
    mss <- Fms/Fms_old
    
    Fbp <- (1/CAPd34S)*((d34S[j,i]-dlsy)*Fwsy+(d34S[j,i]-dlsa)*Fwsa+(d34S[j,i]-dlpy)*Fwpy+(d34S[j,i]-dlpa)*Fwpa+(d34S[j,i]-dlsa)*Fms+(d34S[j,i]-dlpa)*Fmp) #burial flux of pyrite
    
    # end inverse method
    
    Fbs <- Fwpy+Fwpa+Fwsy+Fwsa+Fms+Fmp-Fbp  #burial flux of CaSO4 sulfur
    Fbc <- Fwgy+Fwga+Fwcy+Fwca+Fmc+Fmg-Fbg  #burial flux of carbonate carbon
    
    Fyop <- Fwpa+Fmp  #degassing + weathering flux of pyrite
    Fyos <- Fwsa+Fms  #degassing + weathering flux of CaSO4 sulfur
    Fyog <- Fwga+Fmg  #degassing + weathering flux of organic carbon
    Fyoc <- Fwca+Fmc  #degassing + weathering flux of carbonate carbon
    
    #incorporate the volcanic/non-volcanic components (the "volc" in GEOCARBSULFvolc; see Berner 2006b & 2008)
    Roc <- (Sr[j,i]/10000)+0.7  #87Sr/86Sr of seawater as recorded in carbonates
    Rg <- Rg_570[j,i]-NV[j,i]*(1-fR[j,i]^(1/exp_NV[j,i]))  #87Sr/86Sr of non-volcanic silicates (same as "Rnv" in Berner 2006b)
    A <- ((Rv_570[j,i]-Roc)*fSR[j,i]*Fob[j,i])/(Fbc-Fwcy-Fwca) #note: this is always a negative number because Roc > Rv
    B <- Fwcy/(Fbc-Fwcy-Fwca)
    D <- Fwca/(Fbc-Fwcy-Fwca)
    E <- Fbc/(Fbc-Fwcy-Fwca)
    Xvolc <- (A+B*Rcy+D*Rca-E*Roc+Rg)/(Rg-Rv_570[j,i]) #fraction of total Ca and Mg silicate weathering drived from volcanic rocks
    fvolc <- (VNV[j,i]*Xvolc+1-Xvolc)/(VNV[j,i]*Xvolc_0[j,i]+1-Xvolc_0[j,i])  #volcanic weathering effect at time (t) relative to present-day
    
    #calculate oxygen mass for the time-step
    if (t<570)  {
      oxygen <- oxygen+(Fbg+(15/8)*Fbp)*Dt-(Fwgy+Fwga+Fmg)*Dt-(15/8)*(Fwpy+Fwpa+Fmp)*Dt
    }
    
    #expression that summarizes the chemical weathering of silicates at time (t) relative to the present-day in the absence of climatic effects (other than the effect of CO2 on carbonate weathering); it is calculated by dividing total silicate chemical weathering at time (t)--as determined by mass balance between burial and degassing (numerator)--by the non-climatic processes that affect silicate chemical weathering multiplied by the present-day chemical weathering flux, Fwsi_0 (denominator) (see equation 5.2 in Berner 2004) (called "fB" in BASIC scripts)
    fwsi_no_climate <- (Fbc-Fwcy-Fwca)/(((fAw_fA[j,i]*fA[j,i]*fD[j,i])^exp_fD[j,i])*fE*fR[j,i]*fvolc*Fwsi_0[j,i])
    
    #test for failed runs (negative numbers, NaN's, or oxygen <5% or >50%)
    # !!! change the lower limit is 0. Now delete the limit, saying that O2 is 8 if Oxygen is less than 0. Change <=0 to 0 as weathering can be zero
    if (min(B,D,E,fwsi_no_climate,Fwpy,Fwsy,Fwgy,Fwcy,Fwpa,Fwsa,Fwga,Fwca,Fmp,Fms,Fmg,Fmc,CAPd13C,CAPd34S,Fbp,Fbg,Fbs,Fbc,Spy,Ssy,Gy,Cy,Ca,Xvolc,fvolc)<0 || is.nan(fwsi_no_climate+A+B+D+E+Fwpy+Fwsy+Fwgy+Fwcy+Fwpa+Fwsa+Fwga+Fwca+Fmp+Fms+Fmg+Fmc+CAPd13C+CAPd34S+Fbp+Fbg+Fbs+Fbc+Spy+Ssy+Gy+Cy+Ca+dlpy+dlpa+dlsy+dlsa+dlgy+dlga+dlcy+dlca+Rcy+Rca+Xvolc+fvolc+GCM+oxygen)|| 100*(oxygen/(oxygen+143))<0 || 100*(oxygen/(oxygen+143))>50)  {
      failed_run <- TRUE
      break()
      
    }
    
    #calculate atmospheric CO2 (ppm) through iterative convergence (see FORTRAN scripts of Park & Royer 2011 for details); this is done by calculating the climatic factors that affect silicate weathering rates normalized to the present-day (fwsi_climate, calculated below; called "Fbbs" in BASIC scripts); fwsi_climate combines the expressions "fBt(CO2)" and "fBb(CO2)" in Berner (2004) (see his equations 2.6-7, 2.29, & 5.2); through inversion of fwsi_climate ("W", "V", & "X" below) and comparison to fwsi_no_climate (calculated above), RCO2 can be determined; see pp. 72-76 in Berner (2004) for details; iterative convergence is necessary because RCO2--taken from the previous time-step--is needed to calculate some of the dependent parameters; the initial calculation of the new time-step for RCO2 is therefore not likely to be correct
    RCO2_old <- 2*RCO2  #this is here simply to ensure that the convergence isn't satisfied on the first step
    iteration_count <- 0
    if (t<=570 & t>380 & failed_run==FALSE)  {
      while  (abs(RCO2/RCO2_old-1) > 0.01) {
        iteration_count <- iteration_count+1
        RCO2_old <- RCO2
        fwsi_climate <- (RCO2^(exp_fnBb[j,i]+ACT[j,i]*GCM))*((1+RT[j,i]*GCM*log(RCO2)-RT[j,i]*Ws[j,i]*(t/570)+RT[j,i]*GEOG[j,i])^exp_fD[j,i])*exp(-ACT[j,i]*Ws[j,i]*(t/570))*exp(ACT[j,i]*GEOG[j,i])          
        W <- ((exp_fnBb[j,i]+ACT[j,i]*GCM)*RCO2^(-exp_fnBb[j,i]+ACT[j,i]*GCM))*((1+RT[j,i]*GCM*log(RCO2)-RT[j,i]*Ws[j,i]*(t/570)+RT[j,i]*GEOG[j,i])^exp_fD[j,i])*exp(-ACT[j,i]*Ws[j,i]*(t/570))*exp(ACT[j,i]*GEOG[j,i])
        V <- (RCO2^(exp_fnBb[j,i]+ACT[j,i]*GCM))*exp_fD[j,i]*((1+RT[j,i]*GCM*log(RCO2)-RT[j,i]*Ws[j,i]*(t/570)+RT[j,i]*GEOG[j,i])^-(1-exp_fD[j,i]))*(RT[j,i]*GCM/RCO2)*exp(-ACT[j,i]*Ws[j,i]*t/570)*exp(ACT[j,i]*GEOG[j,i])
        if (is.nan(fwsi_climate+W+V)==TRUE || iteration_count==iteration_threshold)  {
          failed_run <- TRUE
          break()
        }
        if (RCO2 > ((fwsi_climate - fwsi_no_climate)/(W+V))) {
          RCO2 <- RCO2-0.9*((fwsi_climate - fwsi_no_climate)/(W+V))  #damp the iteration to avoid overshoot
        } else {
          RCO2 <- 0.2*RCO2  #damp the iteration (convert the iteration to geometric shrinkage to avoid nonpositive value in overshoot)
        }
      } #end of while loop
    } #end of t>380 loop
    
    if (t<=380 & t>350 & failed_run==FALSE) { #the expressions for this time interval are more complex because the effects of plants on weathering are linearly mixed in; this helps to prevent model failure
      while  (abs(RCO2/RCO2_old-1) > 0.01) {
        iteration_count <- iteration_count+1
        RCO2_old <- RCO2
        fwsi_climate_old <- (RCO2^(exp_fnBb[j,i]+ACT[j,i]*GCM))*((1+RT[j,i]*GCM*log(RCO2)-RT[j,i]*Ws[j,i]*(t/570)+RT[j,i]*GEOG[j,i])^exp_fD[j,i])*exp(-ACT[j,i]*Ws[j,i]*(t/570))*exp(ACT[j,i]*GEOG[j,i])
        W_old <- (exp_fnBb[j,i]+ACT[j,i]*GCM)*RCO2^(-exp_fnBb[j,i]+ACT[j,i]*GCM)*((1+RT[j,i]*GCM*log(RCO2)-RT[j,i]*Ws[j,i]*(t/570)+RT[j,i]*GEOG[j,i])^exp_fD[j,i])*exp(-ACT[j,i]*Ws[j,i]*(t/570))*exp(ACT[j,i]*GEOG[j,i])
        V_old <- RCO2^(exp_fnBb[j,i]+ACT[j,i]*GCM)*exp_fD[j,i]*((1+RT[j,i]*GCM*log(RCO2)-RT[j,i]*Ws[j,i]*(t/570)+RT[j,i]*GEOG[j,i])^-(1-exp_fD[j,i]))*(RT[j,i]*GCM/RCO2)*exp(-ACT[j,i]*Ws[j,i]*(t/570))*exp(ACT[j,i]*GEOG[j,i])
        
        fwsi_climate_new <- ((2^FERT[j,i])*RCO2^(FERT[j,i]+ACT[j,i]*GCM))*((1+RCO2)^(-FERT[j,i]))*((1+RT[j,i]*GCM*log(RCO2)-RT[j,i]*Ws[j,i]*(t/570)+RT[j,i]*GEOG[j,i])^exp_fD[j,i])*exp(-ACT[j,i]*Ws[j,i]*(t/570))*exp(ACT[j,i]*GEOG[j,i])
        W_new <- (2^FERT[j,i])*(FERT[j,i]+ACT[j,i]*GCM)*(RCO2^(FERT[j,i]+ACT[j,i]*GCM-1))*((1+RCO2)^-FERT[j,i])*((1+RT[j,i]*GCM*log(RCO2)-RT[j,i]*Ws[j,i]*(t/570)+RT[j,i]*GEOG[j,i])^exp_fD[j,i])*exp(-ACT[j,i]*Ws[j,i]*(t/570))*exp(ACT[j,i]*GEOG[j,i])
        V_new <- (-FERT[j,i]*(1+RCO2)^-(1+FERT[j,i]))*((2^FERT[j,i])*RCO2^(FERT[j,i]+ACT[j,i]*GCM))*((1+RT[j,i]*GCM*log(RCO2)-RT[j,i]*Ws[j,i]*(t/570)+RT[j,i]*GEOG[j,i])^exp_fD[j,i])*exp(-ACT[j,i]*Ws[j,i]*(t/570))*exp(ACT[j,i]*GEOG[j,i])
        X_new <- exp_fD[j,i]*((1+RT[j,i]*GCM*log(RCO2)-RT[j,i]*Ws[j,i]*(t/570)+RT[j,i]*GEOG[j,i])^-(1-exp_fD[j,i]))*(RT[j,i]*GCM/RCO2)*(2^FERT[j,i]*RCO2^(FERT[j,i]+ACT[j,i]*GCM))*((1+RCO2)^-FERT[j,i])*exp(-ACT[j,i]*Ws[j,i]*(t/570))*exp(ACT[j,i]*GEOG[j,i])
        
        fwsi_climate <- (t-350)/31*fwsi_climate_old+(381-t)/31*fwsi_climate_new
        Fw_v_x <- (t-350)/31*(W_old + V_old)+(381-t)/31*(W_new + V_new + X_new)
        
        if (is.nan(fwsi_climate_old + W_old + V_old + fwsi_climate_new + W_new + V_new + X_new + fwsi_climate + Fw_v_x)==TRUE  || iteration_count==iteration_threshold)  {
          failed_run <- TRUE
          break()
        }
        if (RCO2>((fwsi_climate - fwsi_no_climate)/Fw_v_x))  {
          RCO2 <- RCO2-0.9*((fwsi_climate - fwsi_no_climate)/Fw_v_x)  #damp the iteration to avoid overshoot
        } else {
          RCO2 <- 0.2*RCO2  #damp the iteration (convert the iteration to geometric shrinkage to avoid nonpositive value in overshoot)
        }
      } #end of while loop
    } #end of t<=380 & t>350 loop
    
    if (t<=350 & failed_run==FALSE) {
      while  (abs(RCO2/RCO2_old-1) > 0.01) {
        iteration_count <- iteration_count+1
        RCO2_old <- RCO2
        fwsi_climate <- ((2^FERT[j,i])*RCO2^(FERT[j,i]+ACT[j,i]*GCM))*((1+RCO2)^(-FERT[j,i]))*((1+RT[j,i]*GCM*log(RCO2)-RT[j,i]*Ws[j,i]*(t/570)+RT[j,i]*GEOG[j,i])^exp_fD[j,i])*exp(-ACT[j,i]*Ws[j,i]*(t/570))*exp(ACT[j,i]*GEOG[j,i])
        W <- (2^FERT[j,i])*(FERT[j,i]+ACT[j,i]*GCM)*(RCO2^(FERT[j,i]+ACT[j,i]*GCM-1))*((1+RCO2)^(-FERT[j,i]))*((1+RT[j,i]*GCM*log(RCO2)-RT[j,i]*Ws[j,i]*(t/570)+RT[j,i]*GEOG[j,i])^exp_fD[j,i])*exp(-ACT[j,i]*Ws[j,i]*(t/570))*exp(ACT[j,i]*GEOG[j,i])
        V <- (-FERT[j,i]*(1+RCO2)^-(1+FERT[j,i]))*((2^FERT[j,i])*RCO2^(FERT[j,i]+ACT[j,i]*GCM))*((1+RT[j,i]*GCM*log(RCO2)-RT[j,i]*Ws[j,i]*(t/570)+RT[j,i]*GEOG[j,i])^exp_fD[j,i])*exp(-ACT[j,i]*Ws[j,i]*(t/570))*exp(ACT[j,i]*GEOG[j,i])
        X <- exp_fD[j,i]*((1+RT[j,i]*GCM*log(RCO2)-RT[j,i]*Ws[j,i]*(t/570)+RT[j,i]*GEOG[j,i])^-(1-exp_fD[j,i]))*(RT[j,i]*GCM/RCO2)*(2^FERT[j,i]*RCO2^(FERT[j,i]+ACT[j,i]*GCM))*((1+RCO2)^-FERT[j,i])*exp(-ACT[j,i]*Ws[j,i]*(t/570))*exp(ACT[j,i]*GEOG[j,i])
        if (is.nan(fwsi_climate + W+V+X)==TRUE  || iteration_count==iteration_threshold)  {
          failed_run <- TRUE
          break()
        }
        if (RCO2 > ((fwsi_climate - fwsi_no_climate)/(W+V+X))) {
          RCO2 <- RCO2-0.9*((fwsi_climate - fwsi_no_climate)/(W+V+X))  #damp the iteration to avoid overshoot
        } else {
          RCO2 <- 0.2*RCO2  #damp the iteration (convert the iteration to geometric shrinkage to avoid nonpositive value in overshoot)
        }
      } #end of while loop
    } #end of t<=350 loop
    
    if (failed_run == TRUE) {
      break()
    }
    
    
    if (RCO2<0.6 || RCO2>200) { #test for failed runs when CO2 < 150 ppm or > 50000 ppm
      failed_run <- TRUE
      break()
      
      
    }
    
    #calculate new masses for the next time-step
    Spy <- Spy+(Fbp-Fwpy-Fyop)*Dt  #mass of young pyrite sulfur
    Ssy <- Ssy+(Fbs-Fwsy-Fyos)*Dt  #mass of young CaSO4 sulfur
    Gy <- Gy+(Fbg-Fwgy-Fyog)*Dt  #mass of young crustal organic carbon
    Cy <- Cy+(Fbc-Fwcy-Fyoc)*Dt  #mass of young crustal carbonate carbon
    Ca <- CT[j,i]-Gy-Ga[j,i]-Cy-COC[j,i]  #mass of old crustal carbonate carbon
    
    #calculate new isotopic values for the next time-step
    dlpy <- dlpy+((d34S[j,i]-dlpy-CAPd34S)*Fbp/Spy)*Dt  #d34S of young pyrite sulfur
    dlpa <- dlpa+(Fyop*(dlpy-dlpa)/Spa[j,i])*Dt  #d34S of old pyrite sulfur
    dlsy <- dlsy+((d34S[j,i]-dlsy)*Fbs/Ssy)*Dt  #d34S value of young CaSO4 sulfur
    dlsa <- dlsa+(Fyos*(dlsy-dlsa)/Ssa[j,i])*Dt  #d34S value of old CaSO4 sulfur
    dlgy <- dlgy+((d13C[j,i]-dlgy-CAPd13C)*Fbg/Gy)*Dt  #d13C of young organic matter
    dlga <- dlga+(Fyog*(dlgy-dlga)/Ga[j,i])*Dt  #d13C of old organic matter
    dlcy <- dlcy+((d13C[j,i]-dlcy)*Fbc/Cy)*Dt  #d13C value of young crustal carbonate carbon
    dlca <- dlca+(Fyoc*(dlcy-dlca)/Ca)*Dt  #d13C value of old crustal carbonate carbon
    Rcy <- Rcy+((Roc-Rcy)*Fbc/Cy)*Dt  #87Sr/86Sr of young carbonates undergoing weathering over the time-step
    Rca <- Rca+((Rcy-Rca)*Fyoc/Ca)*Dt  #87Sr/86Sr of old carbonates undergoing weathering over the time-step
    
    
    #fill in the calculated values for O2 (%) and CO2 (ppm) (by default, failed runs will be filled with "NAs")        
    if (failed_run==FALSE) {
      CO2_resamples[j,i] <- RCO2*250  #RCO2 is multiplied by 250 (Pleistocene mean CO2) to convert to ppm units
      O2_resamples[j,i] <- 100*(oxygen/(oxygen+143))  #converts O2 mass to O2 (%)   
      
      orgb_resamples[j,i] <- Fbg  
      
      carb_resamples[j,i] <- Fbc
      
      forg_resamples[j,i] <- Fbg / (Fbg + Fbc)
      
      pyb_resamples[j,i] <- Fbp
      
      sulfb_resamples[j,i] <- Fbs
      
      fpy_resamples[j,i] <- Fbp / (Fbp + Fbs)
      
      orgwy_resamples[j,i] <- Fwgy
      
      orgwa_resamples[j,i] <- Fwga
      
      carwy_resamples[j,i] <- Fwcy   
      
      carwa_resamples[j,i] <- Fwca
      
      pywy_resamples[j,i] <- Fwpy
      
      pywa_resamples[j,i] <- Fwpa 
      
      sulfwy_resamples[j,i] <- Fwsy
      
      sulfwa_resamples[j,i] <- Fwsa
      
      morg_resamples[j,i] <- Fmg
      
      mcar_resamples[j,i] <- Fmc
      
      mpy_resamples[j,i] <- Fmp
      
      msulf_resamples[j,i] <- Fms
      
      orgy_resamples[j,i] <- Gy
      
      orga_resamples[j,i] <- Ga[1,i]
      
      cary_resamples[j,i] <- Cy
      
      cara_resamples[j,i] <- Ca
      
      pyy_resamples[j,i] <- Spy
      
      pya_resamples[j,i] <- Spa[1,i]
      
      sulfy_resamples[j,i] <- Ssy
      
      sulfa_resamples[j,i] <- Ssa[1,i]
      
      TC_resamples[j,i] <- Gy + Ga[1,i] + Cy + Ca
      
      TS_resamples[j,i] <- Spy + Spa[1,i] + Ssy + Ssa[1,i]
      
      dorgy_resamples[j,i] <- dlgy
      
      dorga_resamples[j,i] <- dlga
      
      dcary_resamples[j,i] <- dlcy
      
      dcara_resamples[j,i] <- dlca
      
      dpyy_resamples[j,i] <- dlpy
      
      dpya_resamples[j,i] <- dlpa 
      
      dsulfy_resamples[j,i] <- dlsy 
      
      dsulfa_resamples[j,i] <- dlsa 
      
      CAPdC_resamples[j,i] <- CAPd13C
      
      CAPdS_resamples[j,i] <- CAPd34S
      
      input.C_resamples[j,i] <- Fwgy + Fwga + Fwcy + Fwca + Fmg + Fmc
      
      input.dC_resamples[j,i] <- (Fwgy * dlgy +  (Fwga + Fmg) * dlga +  Fwcy * dlcy + (Fwca + Fmc) * dlca) / input.C_resamples[j,i]
      
      input.S_resamples[j,i] <- Fwpy + Fwpa + Fwsy + Fwsa + Fmp + Fms
      
      input.dS_resamples[j,i] <- (Fwpy * dlpy +  (Fwpa + Fmp) * dlpa +  Fwsy * dlsy + (Fwsa + Fms) * dlsa) / input.S_resamples[j,i]        
      
      dCcar_resamples[j,i] <- s_dCcar
      
      T_org_weather_resamples[j,i] <- Fwgy + Fwga
      
      T_car_weather_resamples[j,i] <- Fwcy + Fwca
      
      T_py_weather_resamples[j,i] <- Fwpy + Fwpa
      
      T_sulf_weather_resamples[j,i] <- Fwsy + Fwsa
      
      soo_resamples[j,i] <- oo  
      
      scc_resamples[j,i] <- cc
      
      smoo_resamples[j,i] <- moo
      
      smcc_resamples[j,i] <- mcc
      
      spp_resamples[j,i] <- pp
      
      sss_resamples[j,i] <- ss
      
      smpp_resamples[j,i] <- mpp
      
      smss_resamples[j,i] <- mss
      
    }  else {
      RCO2 <- 10 #reset RCO2 seed for next run
    }  #end of if..else loop
    
    if (t==0 & (is.nan(oxygen+RCO2) || 100*(oxygen/(oxygen+143))<19 || 100*(oxygen/(oxygen+143))>23 || RCO2<0.8 || RCO2>1.2))  {
      failed_run <- TRUE
      CO2_resamples[,i] <- NA; O2_resamples[,i] <- NA; orgb_resamples[,i] <- NA; carb_resamples[,i] <- NA; forg_resamples[,i] <- NA; pyb_resamples[,i] <- NA; sulfb_resamples[,i] <- NA; fpy_resamples[,i] <- NA; orgwy_resamples[,i] <- NA; orgwa_resamples[,i] <- NA; carwy_resamples[,i] <- NA; carwa_resamples[,i] <- NA; pywy_resamples[,i] <- NA; pywa_resamples[,i] <- NA; sulfwy_resamples[,i] <- NA; sulfwa_resamples[,i] <- NA; morg_resamples[,i] <- NA; mcar_resamples[,i] <- NA; mpy_resamples[,i] <- NA; msulf_resamples[,i] <- NA; orgy_resamples[,i] <- NA; orga_resamples[,i] <- NA; cary_resamples[,i] <- NA; cara_resamples[,i] <- NA; pyy_resamples[,i] <- NA; pya_resamples[,i] <- NA; sulfy_resamples[,i] <- NA; sulfa_resamples[,i] <- NA; TC_resamples[,i] <- NA; TS_resamples[,i] <- NA; dorgy_resamples[,i] <- NA; dorga_resamples[,i] <- NA; dcary_resamples[,i] <- NA; dcara_resamples[,i] <- NA; dpyy_resamples[,i] <- NA; dpya_resamples[,i] <- NA; dsulfy_resamples[,i] <- NA; dsulfa_resamples[,i] <- NA; CAPdC_resamples[,i] <- NA; CAPdS_resamples[,i] <- NA; input.C_resamples[,i] <- NA; input.dC_resamples[,i] <- NA; input.S_resamples[,i] <- NA; input.dS_resamples[,i] <- NA; dCcar_resamples[,i] <- NA; T_org_weather_resamples[,i] <- NA; T_py_weather_resamples[,i] <- NA; T_car_weather_resamples[,i] <- NA; T_sulf_weather_resamples[,i] <- NA; soo_resamples[,i] <- NA; scc_resamples[,i] <- NA; smoo_resamples[,i] <- NA; smcc_resamples[,i] <- NA; spp_resamples[,i] <- NA; sss_resamples[,i] <- NA; smpp_resamples[,i] <- NA; smss_resamples[,i] <- NA
      
      break()
      
    }
    
    
  } #end of age loop (i)
  
  cat(i, "finished", "\n")
} # end resample loop


########################################
### code for outputting summary file ###
########################################

#calculate mean and percentiles for O2 and CO2
for (j in 1:ageN) {
  CO2[j,] <- mean(CO2_resamples[j,], na.rm=TRUE)
  O2[j,] <- mean(O2_resamples[j,], na.rm=TRUE)
  
  orgb[j,] <- mean(orgb_resamples[j,], na.rm=TRUE)
  
  carb[j,] <- mean(carb_resamples[j,], na.rm=TRUE)    
  
  forg[j,] <- mean(forg_resamples[j,], na.rm=TRUE) 
  
  pyb[j,] <- mean(pyb_resamples[j,], na.rm=TRUE)
  
  sulfb[j,] <- mean(sulfb_resamples[j,], na.rm=TRUE)  
  
  fpy[j,] <- mean(fpy_resamples[j,], na.rm=TRUE) 
  
  orgwy[j,] <- mean(orgwy_resamples[j,], na.rm=TRUE)
  
  orgwa[j,] <- mean(orgwa_resamples[j,], na.rm=TRUE) 
  
  carwy[j,] <- mean(carwy_resamples[j,], na.rm=TRUE)
  
  carwa[j,] <- mean(carwa_resamples[j,], na.rm=TRUE)
  
  pywy[j,] <- mean(pywy_resamples[j,], na.rm=TRUE)
  
  pywa[j,] <- mean(pywa_resamples[j,], na.rm=TRUE)
  
  sulfwy[j,] <- mean(sulfwy_resamples[j,], na.rm=TRUE)
  
  sulfwa[j,] <- mean(sulfwa_resamples[j,], na.rm=TRUE)
  
  morg[j,] <- mean(morg_resamples[j,], na.rm=TRUE)
  
  mcar[j,] <- mean(mcar_resamples[j,], na.rm=TRUE)
  
  mpy[j,] <- mean(mpy_resamples[j,], na.rm=TRUE)
  
  msulf[j,] <- mean(msulf_resamples[j,], na.rm=TRUE)
  
  orgy[j,] <- mean(orgy_resamples[j,], na.rm=TRUE)
  
  orga[j,] <- mean(orga_resamples[j,], na.rm=TRUE) 
  
  cary[j,] <- mean(cary_resamples[j,], na.rm=TRUE)
  
  cara[j,] <- mean(cara_resamples[j,], na.rm=TRUE)
  
  pyy[j,] <- mean(pyy_resamples[j,], na.rm=TRUE)
  
  pya[j,] <- mean(pya_resamples[j,], na.rm=TRUE)
  
  sulfy[j,] <- mean(sulfy_resamples[j,], na.rm=TRUE)
  
  sulfa[j,] <- mean(sulfa_resamples[j,], na.rm=TRUE)
  
  TC[j,] <- mean(TC_resamples[j,], na.rm=TRUE)
  
  TS[j,] <- mean(TS_resamples[j,], na.rm=TRUE)
  
  dorgy[j,] <- mean(dorgy_resamples[j,], na.rm=TRUE)
  
  dorga[j,] <- mean(dorga_resamples[j,], na.rm=TRUE)
  
  dcary[j,] <- mean(dcary_resamples[j,], na.rm=TRUE)
  
  dcara[j,] <- mean(dcara_resamples[j,], na.rm=TRUE)
  
  dpyy[j,] <- mean(dpyy_resamples[j,], na.rm=TRUE)
  
  dpya[j,] <- mean(dpya_resamples[j,], na.rm=TRUE)
  
  dsulfy[j,] <- mean(dsulfy_resamples[j,], na.rm=TRUE)
  
  dsulfa[j,] <- mean(dsulfa_resamples[j,], na.rm=TRUE)
  
  CAPdC[j,] <- mean(CAPdC_resamples[j,], na.rm=TRUE)
  
  CAPdS[j,] <- mean(CAPdS_resamples[j,], na.rm=TRUE)
  
  input.C[j,] <- mean(input.C_resamples[j,], na.rm=TRUE)
  
  input.dC[j,] <- mean(input.dC_resamples[j,], na.rm=TRUE)
  
  input.S[j,] <- mean(input.S_resamples[j,], na.rm=TRUE)
  
  input.dS[j,] <- mean(input.dS_resamples[j,], na.rm=TRUE) 
  
  dCcar[j,] <- mean(dCcar_resamples[j,], na.rm=TRUE)
  
  T_org_weather[j,] <- mean(T_org_weather_resamples[j,], na.rm=TRUE)
  
  T_car_weather[j,] <- mean(T_car_weather_resamples[j,], na.rm=TRUE)
  
  T_py_weather[j,] <- mean(T_py_weather_resamples[j,], na.rm=TRUE)
  
  T_sulf_weather[j,] <- mean(T_sulf_weather_resamples[j,], na.rm=TRUE)
  
  soo[j,] <- mean(soo_resamples[j,], na.rm=TRUE)
  
  scc[j,] <- mean(scc_resamples[j,], na.rm=TRUE)
  
  smoo[j,] <- mean(smoo_resamples[j,], na.rm=TRUE)
  
  smcc[j,] <- mean(smcc_resamples[j,], na.rm=TRUE)
  
  spp[j,] <- mean(spp_resamples[j,], na.rm=TRUE)
  
  sss[j,] <- mean(sss_resamples[j,], na.rm=TRUE)
  
  smpp[j,] <- mean(smpp_resamples[j,], na.rm=TRUE)
  
  smss[j,] <- mean(smss_resamples[j,], na.rm=TRUE)
  
  
  failed_runs[j,] <- sum(is.na(CO2_resamples[j,]))/resampleN*100
  
  if (resampleN>1)  {
    for (k in 1:length(percentile_values))  {
      percentiles_CO2[j,k] <- quantile(x=CO2_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      percentiles_O2[j,k] <- quantile(x=O2_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_orgb[j,k] <- quantile(x=orgb_resamples[j,], probs=percentile_values[k], na.rm=TRUE)     
      
      percentiles_carb[j,k] <- quantile(x=carb_resamples[j,], probs=percentile_values[k], na.rm=TRUE) 
      
      percentiles_forg[j,k] <- quantile(x=forg_resamples[j,], probs=percentile_values[k], na.rm=TRUE)  
      
      percentiles_pyb[j,k] <- quantile(x=pyb_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_sulfb[j,k] <- quantile(x=sulfb_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_fpy[j,k] <- quantile(x=fpy_resamples[j,], probs=percentile_values[k], na.rm=TRUE) 
      
      percentiles_orgwy[j,k] <- quantile(x=orgwy_resamples[j,], probs=percentile_values[k], na.rm=TRUE)  
      
      percentiles_orgwa[j,k] <- quantile(x=orgwa_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_carwy[j,k] <- quantile(x=carwy_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_carwa[j,k] <- quantile(x=carwa_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_pywy[j,k] <- quantile(x=pywy_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_pywa[j,k] <- quantile(x=pywa_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_sulfwy[j,k] <- quantile(x=sulfwy_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_sulfwa[j,k] <- quantile(x=sulfwa_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_morg[j,k] <- quantile(x=morg_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_mcar[j,k] <- quantile(x=mcar_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_mpy[j,k] <- quantile(x=mpy_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_msulf[j,k] <- quantile(x=msulf_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_orgy[j,k] <- quantile(x=orgy_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_orga[j,k] <- quantile(x=orga_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_cary[j,k] <- quantile(x=cary_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_cara[j,k] <- quantile(x=cara_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_pyy[j,k] <- quantile(x=pyy_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_pya[j,k] <- quantile(x=pya_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_sulfy[j,k] <- quantile(x=sulfy_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_sulfa[j,k] <- quantile(x=sulfa_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_TC[j,k] <- quantile(x=TC_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_TS[j,k] <- quantile(x=TS_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_dorgy[j,k] <- quantile(x=dorgy_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_dorga[j,k] <- quantile(x=dorga_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_dcary[j,k] <- quantile(x=dcary_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_dcara[j,k] <- quantile(x=dcara_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_dpyy[j,k] <- quantile(x=dpyy_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_dpya[j,k] <- quantile(x=dpya_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_dsulfy[j,k] <- quantile(x=dsulfy_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_dsulfa[j,k] <- quantile(x=dsulfa_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_CAPdC[j,k] <- quantile(x=CAPdC_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_CAPdS[j,k] <- quantile(x=CAPdS_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_input.C[j,k] <- quantile(x=input.C_resamples[j,], probs=percentile_values[k], na.rm=TRUE)     
      
      percentiles_input.dC[j,k] <- quantile(x=input.dC_resamples[j,], probs=percentile_values[k], na.rm=TRUE)  
      
      percentiles_input.S[j,k] <- quantile(x=input.S_resamples[j,], probs=percentile_values[k], na.rm=TRUE)     
      
      percentiles_input.dS[j,k] <- quantile(x=input.dS_resamples[j,], probs=percentile_values[k], na.rm=TRUE) 
      
      percentiles_dCcar[j,k] <- quantile(x=dCcar_resamples[j,], probs=percentile_values[k], na.rm=TRUE) 
      
      
      percentiles_T_org_weather[j,k] <- quantile(x=T_org_weather_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_T_car_weather[j,k] <- quantile(x=T_car_weather_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_T_py_weather[j,k] <- quantile(x=T_py_weather_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_T_sulf_weather[j,k] <- quantile(x=T_sulf_weather_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      
      percentiles_soo[j,k] <- quantile(x=soo_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_scc[j,k] <- quantile(x=scc_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_smoo[j,k] <- quantile(x=smoo_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_smcc[j,k] <- quantile(x=smcc_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_spp[j,k] <- quantile(x=spp_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_sss[j,k] <- quantile(x=sss_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_smpp[j,k] <- quantile(x=smpp_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
      percentiles_smss[j,k] <- quantile(x=smss_resamples[j,], probs=percentile_values[k], na.rm=TRUE)
      
    }
  } #end of if statement
} #end of loop for filling the mean and percentile arrays

#merge results and export summary file to working directory


if (resampleN==1)  {
  # !!! add pyb and sulb
  GEOCARB_output <- cbind(age, failed_runs, CO2, O2, orgb, carb, forg, pyb, sulfb, fpy, orgwy, orgwa, carwy, carwa, pywy, pywa, sulfwy, sulfwa, morg, mcar, mpy, msulf, orgy, orga, cary, cara, pyy, pya, sulfy, sulfa, TC, TS, dorgy, dorga, dcary, dcara, dpyy, dpya, dsulfy, dsulfa, CAPdC, CAPdS, input.C, input.dC, input.S, input.dS, dCcar, T_org_weather, T_car_weather, T_py_weather, T_sulf_weather, soo, scc, smoo, smcc, spp, sss, smpp, smss)
} else {
  # !!! add pyb and sulb and everything else
  GEOCARB_output <- cbind(age, failed_runs, CO2, percentiles_CO2, O2, percentiles_O2, orgb, percentiles_orgb, carb, percentiles_carb, forg, percentiles_forg, pyb, percentiles_pyb, sulfb, percentiles_sulfb, fpy, percentiles_fpy, orgwy, percentiles_orgwy, orgwa, percentiles_orgwa, carwy, percentiles_carwy, carwa, percentiles_carwa, pywy, percentiles_pywy, pywa, percentiles_pywa, sulfwy, percentiles_sulfwy, sulfwa, percentiles_sulfwa, morg, percentiles_morg, mcar, percentiles_mcar, mpy, percentiles_mpy, msulf, percentiles_msulf, orgy, percentiles_orgy, orga, percentiles_orga, cary, percentiles_cary, cara, percentiles_cara, pyy, percentiles_pyy, pya, percentiles_pya, sulfy, percentiles_sulfy, sulfa, percentiles_sulfa, TC, percentiles_TC, TS, percentiles_TS, dorgy, percentiles_dorgy, dorga, percentiles_dorga, dcary, percentiles_dcary, dcara, percentiles_dcara, dpyy, percentiles_dpyy, dpya, percentiles_dpya, dsulfy, percentiles_dsulfy, dsulfa, percentiles_dsulfa, CAPdC, percentiles_CAPdC, CAPdS, percentiles_CAPdS, input.C, percentiles_input.C, input.dC, percentiles_input.dC, input.S, percentiles_input.S, input.dS, percentiles_input.dS, dCcar, percentiles_dCcar, T_org_weather, percentiles_T_org_weather, T_car_weather, percentiles_T_car_weather, T_py_weather, percentiles_T_py_weather, T_sulf_weather, percentiles_T_sulf_weather, soo, percentiles_soo, scc, percentiles_scc, smoo, percentiles_smoo, smcc, percentiles_smcc, spp, percentiles_spp, sss, percentiles_sss, smpp, percentiles_smpp, smss, percentiles_smss)
  
  
} #end of if...else statement



if (Godderis == FALSE) {
  file_out <- paste0('GEOCARB_output_', toString(input[40, 'mean']), '_', 'Berner_', toString(resampleN), '.csv')  
} else {
  file_out <- paste0('GEOCARB_use_ds_records_10_smooth_50%_sulfate_O2_forced_solve_8_scaling_factors_monte_fail_test_', scenario, '.csv')
}

write.csv(GEOCARB_output, file_out)


plot(age, O2, type='p')
lines(age, O2)
grid()

# plot(age, dCcar, type='p')
# lines(age, dCcar)
# grid()


time_end <- Sys.time()

print(time_end - time_start)
