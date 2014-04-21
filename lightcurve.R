# This program models a light curve for a HD52721 (GU CMa) system
# on circles and ellipses

# model parameters:
eps <- 1e-9

# P -orbital period (days)
P <- 1.610155

# apmag -apparent magnitude integral (m)
apmag <- 6.6

# par -parralax (mas)
par <- 2.0

# sep -angular separation of components (10 * mas)
sep <- 6.5

# r -radius of the orbit
r <- 3

# rad1 -radius of star A
rad1 <- 1.1

# rad2 -radius of star B
rad2 <- 1

# lux1 -luminosity coef of star A
lux1 <- 1

# lux2 -luminosity coef of star B
lux2 <- 0.8

# incl -inclination angle
incl <- 69.2 * pi / 180

# initial phase
phase <- 9.2 * pi / 180

## function: vdist
# This function will return a visible distance between stars
# consider an inclination angle.

# Input parameters:
# incl -inclination angle
# pa -positional angle

vdist <- function(incl, pa)
{
	return(r * sqrt( cos(pa)**2 + cos(incl)**2 * sin(pa)**2 ) )
}

## function: v_circl_square
# This function will return a visible square projection (spheric star)
# Input parameters:
# t -time (phase of a period T)
# incl -inclination angle

v_circle_square <- function(t, incl)
{
	# square terms
	s <- 0
	
	s <- lux1 * pi * rad1**2 + lux2 * pi * rad2**2

	# current positional angle of a star B
	phi <- 2 * pi / P * t + phase

	# visible distance
	vd <- vdist(incl, phi)
  #print("vd: ", vd)
  
	#upper edge
	upper_edge <- rad1 + rad2

	# find min of radii
	min <- rad1
	max <- rad2
	is_rad1_max <- FALSE

	if ( rad2 < rad1 )
	{
		min <- rad2
		max <- rad1 
		is_rad1_max <- TRUE
	}		

	# define coef intensity
	light_down <- 0
	light_up <- 0

	if ( phi - pi <= eps )
	{
		light_down <- lux2
		light_up <- lux1
	}
	else
	{
		light_down <- lux1
		light_up <- lux2
	}

	# stars are contiguous to each other
	if ( abs( upper_edge - vd ) <= eps )
	{	
		return(s)
	}

	# stars are far away from each other
	if ( vd > upper_edge )
	{
		return(s)
	}
		
	# full eclipse (case 5)
	if ( is_rad1_max && vd + rad2 - rad1 <= eps )
	{
		if ( phi - pi <= eps )
		   return( light_up * pi * rad1**2 )
		return( light_up * pi * rad2**2 +
              light_down * pi * rad1**2 -
              light_down * pi * rad2**2)
	} 

	# stars are overlapping (inludes subcases )
	
	arg1 <- ( vd**2 + rad1**2 - rad2**2 ) / 2 / vd / rad1
	arg2 <- ( vd**2 + rad2**2 - rad1**2 ) / 2 / vd / rad2

	psi1 <- acos( arg1 )
	psi2 <- acos( arg2 )

	if ( is_rad1_max )
	{
		segm1 <- rad1**2 * ( 2 * psi1 - sin(2 * psi1) ) / 2
		segm2 <- rad2**2 * ( 2 * psi2 - sin(2 * psi2) ) / 2
		semisquare <- pi * rad2**2 / 2 
		obtuse_arc <- rad2**2 * ( 2 * psi2 ) / 2
		triangle <- rad2**2 * sin( -2 * psi2) / 2
		
		if ( psi2 - pi / 2 < eps )
		{
		   return( s - light_down * (segm1 + segm2) )
		}		   

		if ( abs( psi2 - pi / 2 ) <= eps )
		{
		   return( s - light_down * ( semisquare + segm1 ) )
		}

		return( s - light_down * ( segm1 + obtuse_arc + triangle ) )
	}
	
  return(0)
}

# Funcion: model()
# This function calls a square model
#
# Input parameters:
# t -time (phase of a period T)
#
# Returns:
# modelled square at the "t" moment
model <- function(t)
{
  s <- lux1 * pi * rad1**2 + lux2 * pi * rad2**2
  
  val <- -0.37 + apmag + 2.5*log10(s/v_circle_square(t, incl))
  return(val)
}

# Function: readSmoothData()
#
# This function reads smoothed data.
#
# Returns: a data table of smoothed data
#
readSmoothData <- function()
{
  obs <-  read.table("obs/pmm_all_1610155_15")
  
  s <- lux1 * pi * rad1**2 + lux2 * pi * rad2**2
  
  res <- c()
  for (q in 1:nrow(obs))
    res <- rbind(res, c(obs[q,1], model(obs[q,1]),#-0.38 +
    #                    apmag + 2.5*log10(s/v_circle_square(obs[q,1], incl)),
                        obs[q,2]
                         ))
  write.table(res,"model.dat", col.names=F, row.names=F)
  return(res)
  #return(obs)
}
res <- readSmoothData()

# This function reads raw data
# and returns a data fold
readRawData <- function()
{
  w <- read.table("obs/all_10_13")
  sp <- c()
  for (t in 1:nrow(w))
  {
    sp <- rbind(sp, c(
      t,                    # number
      w[t, 1] %% P,         # folded JD - phase
      w[t, 1],              # JD
      model(w[t, 1] %% P),  # C data
      w[t, 2]               # O data
    ))
  }
  write.table(sp,"rawdata", row.names=F, col.names=F)
  return(sp)
  #return(sp[order(sp[,1]),])
}
rawd <- readRawData()

draw <- function()
{
  
  # open new device 
 # dev.new()
  
  # prepare filename for pdf
  filename <- paste("plot",
                   "incl", incl * 180 / pi,
                   "radii",rad1,rad2,
                   "sep",r,
                   "lux",lux1,lux2,
                   "P",P,
                   sep="_")
  
  # set plot parameters for 4 graphs on a page
  par(mfrow=c(3,2))
  
  # plot Smooth vs Model
  plot(res[,1], res[,2], col=rgb(0,0,0), pch = 18, cex=0.5,
       main="GU CMa light curve",
       ylab="apparent magnitude (m)",
       xlab="period (days)",
       ylim=c(7,6))
  points(res[,1], res[,3], col=rgb(0,0,1), pch =18, cex=0.5)
  #text(0.35,7.5,col=rgb(0,0,0),"model")
  #text(0.35,7.7, col=rgb(0,0,1),"data")
  legend("bottomleft", c("model", "data"), pch=c(18, 18), cex=0.5,
         col=c(rgb(0,0,0), rgb(0,0,1)) )
  text(1.1,6.70, "parameters", adj=0)
  text(1.1,6.75, paste("inclination", incl * 180 / pi), adj=0)
  text(1.1,6.80, paste("radii", rad1, ";", rad2), adj=0)
  text(1.1,6.85, paste("separation", r), adj=0)
  text(1.1,6.90, paste("lux factors", lux1, ";", lux2), adj=0)
  text(1.1,6.95, paste("period", P), adj=0)
  text(1.1,7.0, paste("initial phase", phase * 180 / pi), adj=0)
  
  # plot Residuals Smooth vs Model
  plot(res[,1],abs(res[,2]-res[,3]), pch=".", type="l",
       main="Residuals (smoothed data)",
       ylab="app. magn. (m)",
       xlab="period (days)")
  
  # 
  #w <- read.table("obs/all_10_13")
  #sp <- c()
  #for (t in 1:nrow(w))
  #{
  #  sp <- rbind(sp, c(
  #    t,
  #    w[t, 1] %% P,         # folded JD - phase
  #    w[t, 1],              # JD
  #    model(w[t, 1] %% P),  # C data
  #    w[t, 2],               # O data
  #  ))
  #}
  
  # plot O - C  (first 832 Nodes)
  N <- 832 
  plot(rawd[1:N,3], rawd[1:N,5] - rawd[1:N,4], pch=18,
       cex=0.5,
#       type="l",
       main="Residuals (2010 raw data)",
       ylab="apparent magnitude (m)",
       xlab="JD from 2000",)
  
  
  # fit
  vect <- rawd[1:N,5] - rawd[1:N,4]
  t <- rawd[1:N,3]
  lm.s <- lm(vect ~ sin(t) + cos(t))
  #plot(lm.s)
 # write.table(lm.s, "myfit")
  
  
  # Plot scaled residuals
#   indexes.scaled <- 1:100
#   plot(rawd[indexes.scaled,3],
#        rawd[indexes.scaled,5] - rawd[indexes.scaled,4],
#        pch=18,
#        #type="l",
#        main="Residuals",
#        ylab="apparent magnitude (m)",
#        xlab="period (days)",
#        )
  
  # Plot Lomb-Scargle periodogram
  indexes <- 1:832
  GM.dat <- rawd[indexes,5] - rawd[indexes,4]
  GM.time <- rawd[indexes,3]
  GM.ts <- ts(GM.dat, GM.time)
  
  GM.lsfreq <- lsp(GM.dat,
                        time = GM.time,
                        type="frequency",
                        #from=0.1,
                        #to=50,
  )
  text(20,200, paste("peaks at:",
                     round(GM.lsfreq$peak.at[1],4),
                     round(GM.lsfreq$peak.at[2],4)))
  
  
  # Plot 2013 raw data
  indexes <- 833:nrow(rawd)
  GM2013.dat <- rawd[indexes,5] - rawd[indexes,4]
  GM2013.time <- rawd[indexes,3]
  GM2013.ts <- ts(GM.dat, GM.time)
  
  
  plot(GM2013.time, GM2013.dat, pch=18,
       cex=0.5,
       #       type="l",
       main="Residuals (2013 raw data)",
       ylab="apparent magnitude (m)",
       xlab="JD from 2000",)
  GM2013.lsfreq <- lsp(GM2013.dat,
                   time = GM2013.time,
                   type="frequency",
                   #from=0.1,
                   #to=50,
  )
  text(9,120, paste("peaks at:",
                     round(GM2013.lsfreq$peak.at[1],4),
                     round(GM2013.lsfreq$peak.at[2],4)))
  
  
  # save this plot to pdf
  dev.copy2pdf(
    file = paste0("plots/",filename,".pdf"),
    width=24, height=24)
}

draw()
regression <- function()
{
  # Create time series
  indexes <- 1:832
  GM.dat <- rawd[indexes,5] - rawd[indexes,4]
  GM.time <- rawd[indexes,3]
  GM.ts <- ts(GM.dat, GM.time)
  

  sin.t <- sin(2*pi*GM.time)
  cos.t <- cos(2*pi*GM.time)
  lm.GM <- lm(GM.dat ~ sin.t + cos.t)
  plot(lm.GM)
  
  #plot.ts(GM.ts)
  #acf(GM.ts)
}
#regression()

require("lomb")
lombscargle <- function()
{
  # Create time series
  indexes <- 1:832
  GM.dat <- rawd[indexes,5] - rawd[indexes,4]
  GM.time <- rawd[indexes,3]
  GM.ts <- ts(GM.dat, GM.time)
  
  GM.lombscargle <- lsp(GM.dat,
                        time = GM.time,
                        type="frequency",
                        #from=0.1,
                        #to=50,
                        )
  
  return(GM.lombscargle)
}
#lslist <- lombscargle()