# This file is the Yasso15 model core code 

# Copyright (C) <2017>  <Finnish Meteorological Institute>
#   
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


# Improved and updated version based on
# Yasso07 description Tuomi & Liski 17.3.2008  (Yasso07.pdf)
# and Taru Palosuo’s code in December 2011
#
# This version uses the separate temperature/precipitation dependencies for the N and H compartments
# The parameters were estimated using e.g. the additional global scale Zinke data set
# 
# 
# Last edited 26.4.2018
# - Toni Viskari


#  Instructions  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

# 1) first run the source code for the function with "source(....r)"
# 2) then you can use the Yasso-function by just calling it Yasso15_R_version(...)
# 3) Input for the function:
# 	1. Yasso15Parameters - Yasso parameters as a vector, length 35
#		1-16 matrix A entries: 4*alpha, 12*p
#		17-21 Leaching parameters: w1,...,w5 IGNORED IN THIS FUNCTION
#		22-23 Temperature-dependence parameters for AWE fractions: beta_1, beta_2
#		24-25 Temperature-dependence parameters for N fraction: beta_N1, beta_N2
#		26-27 Temperature-dependence parameters for H fraction: beta_H1, beta_H2
#		28-30 Precipitation-dependence parameters for AWE, N and H fractions: gamma, gamma_N, gamma_H
#		31-32 Humus decomposition parameters: p_H, alpha_H (Note the order!)
#		33-35 Woody parameters: theta_1, theta_2, r
#	2. SimulationTime - time when the result is requested [a]
#       3. MeanTemperature - mean annual temperature [C]
#       4. TemperatureAmplitude - temperature amplitude i.e. (T_max-T_min)/2, [C] 
#       5. Precipitation - annual precipitation [mm]
#       6. InitialCPool - initial C pools of model compartments, length 5, [whatever]
#       7. LitterInput - mean litter input, 5 columns AWENH, must be the same unit as InitialCpool per year
#       8. WoodySize - size of woody litter (for non-woody litter this is 0) [cm]
#       9. Leaching - Parameter value representing the impact of the holes in the litter decomposition bag. Should be zero in Yasso runs that aren't fitted to those specific measurements.
#       10. SS_pred - Logical variable representing if the steady state is directly calculated. Set to TRUE if determining the steady state, FALSE otherwise.

# 4) The function returns the amount of litter as 5-vector (AWENH compartments) at SimulationTime

# NOTE that this function eats only one type of material at the time. So, non-woody and different woody litter
# materials needs to be calculated separately (and finally count together if desired).

# The output of the function is the vector AWENH compartments at the given time since the simulation start


# Basics  BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

# additional R libraries (as needed)
#library(Matrix)  # tai Matrix tms  


# Function definition   FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

Yasso15_R_version <- function(Yasso15Parameters, SimulationTime, MeanTemperature, TemperatureAmplitude, Precipitation, InitialCPool, LitterInput, WoodySize, Leac, SS_pred)
{

# using shorter names for input
theta <- Yasso15Parameters
t <- SimulationTime
climate <- c(MeanTemperature,Precipitation,TemperatureAmplitude)
init <- InitialCPool
b <- LitterInput
d <- WoodySize
leac <- Leac

# temperature annual cycle approximation
te1 <- climate[1]+4*climate[3]/pi*(1/sqrt(2)-1)
te2 <- climate[1]-4*climate[3]/(sqrt(2)*pi)
te3 <- climate[1]+4*climate[3]/pi*(1-1/sqrt(2))
te4 <- climate[1]+4*climate[3]/(sqrt(2)*pi)

# Average temperature dependence
te <- c(te1,te2,te3,te4)
tem <- mean(exp(theta[22]*te+theta[23]*te^2)) 
temN <- mean(exp(theta[24]*te+theta[25]*te^2)) 
temH <- mean(exp(theta[26]*te+theta[27]*te^2)) 

# Precipitation dependence
tem <- tem*(1.0-exp(theta[28]*climate[2]/1000.0)) # precipitation as m
temN <- temN*(1.0-exp(theta[29]*climate[2]/1000.0))
temH <- temH*(1.0-exp(theta[30]*climate[2]/1000.0))

# Size class dependence -- no effect if d == 0.0
size_dep <- min(1.0,(1+theta[33]*d+theta[34]*d^2)^(-abs(theta[35])))

# check rare case where no decomposition happens for some of the compartments
# (basically, if no rain)
if(tem <= 1e-16) 
{
	xt <- init + b*t;
	return(xt);
}

# Calculating matrix A (will work ok despite the sign of alphas)

alpha <- abs(c(theta[1], theta[2], theta[3], theta[4], theta[32]))  # Vector of decomposition rates

# Creating the matrix A_p
row1 <- c(-1, theta[5], theta[6], theta[7], 0)
row2 <- c(theta[8], -1, theta[9], theta[10], 0)
row3 <- c(theta[11], theta[12], -1, theta[13], 0)
row4 <- c(theta[14], theta[15], theta[16], -1, 0)
row5 <- c(theta[31], theta[31], theta[31], theta[31], -1)
A_p <- matrix(c(row1, row2, row3, row4, row5), 5, 5, byrow=TRUE) # A_p is now computed

# computing the diagonal coefficient matrix k
k <- diag(c(tem*alpha[1:3]*size_dep,temN*alpha[4]*size_dep,temH*alpha[5])) # no size effect in humus

A <- A_p%*%k # coefficient matrix A is now computed

# Leaching (no leaching for humus)
for (ii in 1:4){
	A[ii,ii] = A[ii,ii] + leac*Precipitation/1000.
}

if(SS_pred==TRUE){
	# Solve DE directly in steady state conditions
	xt <- as.array(solve(-A,b))
} else {
	# Solve the differential equation x'(t) = A(theta)*x(t) + b, x(0) = init, using the analytical formula
	z1 <- A %*% init + b;
	mexpAt <- expm(A*t)
	z2 <- mexpAt %*% z1 - b
	xt <- as.array(solve(A,z2))
}
return(xt)
    
} # end of Yasso15 function





