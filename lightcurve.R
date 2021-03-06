# This program models a light curve for a HD52721 (GU CMa) system
# with assumption of star's ellipticity

# Orbit parameters
# ORBIT.INCL inclination of the orbit (radians)
# ORBIT.R radius of the component star (conventional units)
# ORBIT.ECC orbit eccentricity
# ORBIT.PHI0 initial angle of the component star at the orbit (radians)
# ORBIT.T0 initial phase of period 
# ORBIT.P period (days)

require("lomb")
require("zoo")

ORBIT.INCL <- 67 * pi / 180
ORBIT.R <- 3
ORBIT.ECC <- 0.0000
#ORBIT.ECC <- 0.5
ORBIT.PHI0 <- 9.2 * pi / 180
ORBIT.T0 <- 0.83
ORBIT.P <- 1.610157

# Component parameters
A.big.axis <- 1.03
A.small.axis <- 1* A.big.axis
A.lux <- 1

B.big.axis <- 0.97
B.small.axis <- 0.93 * B.big.axis
B.lux <- 0.9

# Photometric parameters
app.mag <- 6.6

# Common parameters
ZERO.POINT <- -0.372
indexes.2010 <- 1:832

# model parameters:
eps <- 1e-9

## RAW DATA INDEXES
i.2010 <- 1:832
i.2013 <- 833:2465

i.85.200 <- 85:200
i.24.305 <- 24:305

# P -orbital period (days)
#P <- 1.610155

# apmag -apparent magnitude integral (m)
#apmag <- 6.6

# par -parralax (mas)
#par <- 2.0

# sep -angular separation of components (10 * mas)
#sep <- 6.5

# r -radius of the orbit
#r <- 3

# rad1 -radius of star A
#rad1 <- 1.2

# rad2 -radius of star B
#rad2 <- 1

# lux1 -luminosity coef of star A
#lux1 <- 1

# lux2 -luminosity coef of star B
#lux2 <- 0.85

# incl -inclination angle
#incl <- 69.2 * pi / 180

# initial phase
#phase <- 9.2 * pi / 180

## function: vdist1
# This function will return a visible distance between stars
# consider an inclination angle.

# Input parameters:
# incl -inclination angle
# pa -positional angle

vdist <- function(incl, pa)
{
  a <- ORBIT.R
  b <- sqrt(a^2*(1-ORBIT.ECC^2))
	return(sqrt( a**2 * cos(pa)**2 + b ** 2 * cos(incl)**2 * sin(pa)**2 ) )
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
	
  rad1 <- A.big.axis
  rad2 <- B.big.axis
  lux1 <- A.lux
  lux2 <- B.lux
  
  P <- ORBIT.P
  
  phase <- ORBIT.PHI0
  
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

	if ( rad2 <= rad1 )
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


# Function: v_ellipse_square
# Input:
# t -time (phase of a period T)
# inlc -inlination angle (radians)
# Value: joint square value

v_ellipse_square <- function(t, incl)
{
  # current positional angle of a star B
  phi <- 2.0 * pi / ORBIT.P * t   + ORBIT.PHI0
  
  # axises change its angle
  if ( incl <= pi / 2 )
  {
    # phi.plane - angle of star rotation in plane of orbit
    # from cosines theorem for paralactic triangles
    phi.plane <- acos( cos(phi) * cos( ( (pi/2 - incl)* sin(phi) )))
    
    # calculating axises of A component
    A.a.true <- A.big.axis #* cos(phi.plane)
    A.b.true <- A.small.axis
    
    # calculating axises of B component
    B.a.true <- B.big.axis #* cos(phi.plane)
    B.b.true <- B.small.axis
    
    # phi.new -  of angle of radius-vector with max coordinate
    A.phi.new <- atan(A.b.true / A.a.true * tan(phi.plane))
    B.phi.new <- atan(B.b.true / B.a.true * tan(phi.plane))
    
    # calculating visible axises
    A.a <- sqrt( (A.a.true*cos(A.phi.new))^2 + (A.b.true*sin(A.phi.new))^2  )
    A.b <- A.b.true
    
    B.a <- sqrt( (B.a.true*cos(B.phi.new))^2 + (B.b.true*sin(B.phi.new))^2  )
    B.b <- B.b.true
  }
  
  
  s <- A.lux * pi * A.a * A.b +
    B.lux * pi * B.a * B.b
  
  # visible distance
  vd <- vdist(ORBIT.INCL, phi.plane)

  
  
  #upper edge
  #upper_edge <- A.big.axis + B.big.axis
  upper_edge <- A.a + B.a
  
  # find min of radii
  #min <- A.big.axis
  min <- A.a
  #max <- B.big.axis
  max <- B.a
  is_rad1_max <- FALSE
  
  if ( B.a <= A.a)
  {
    min <- B.a
    max <- A.a
    is_rad1_max <- TRUE
  }		
  
  # define coef intensity
 
  if ( phi - pi <= eps )
  {
    light_down <- B.lux
    light_up <- A.lux
  }
  else
  {
    light_down <- A.lux
    light_up <- B.lux
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
  if ( is_rad1_max && vd + B.a - A.a <= eps)
  {
     if ( phi - pi <= eps )
      return( light_up * pi * A.a * A.b )
     return( light_up * pi * B.a * B.b +
               light_down * pi * A.a * A.b -
               light_down * pi * B.a * B.b)
  } 
  
  # stars are overlapping (inludes subcases )
  
  # Calculate angle for circles first
  #arg1 <- ( vd**2 + rad1**2 - rad2**2 ) / 2 / vd / rad1
  arg1 <- ( vd**2 + A.a**2 - B.a**2 ) / 2 / vd / A.a
  #arg2 <- ( vd**2 + rad2**2 - rad1**2 ) / 2 / vd / rad2
  arg2 <- ( vd**2 + B.a**2 - A.a**2 ) / 2 / vd / B.a
  
  psi1 <- acos( arg1 )
  psi2 <- acos( arg2 )
  
  zeta1 <- atan( A.b / A.a * tan(psi1))
  zeta2 <- atan( B.b / B.a * tan(psi2))
  
  #za <- A.a*cos(zeta1)
  #zb <- B.a*cos(zeta2)
  #za <- sqrt((A.a*cos(zeta1))^2 + (A.b*sin(zeta1))^2) * cos(zeta1)
  #zb <- sqrt((B.a*cos(zeta2))^2 + (B.b*sin(zeta2))^2) * cos(zeta2)
  
  if ( is_rad1_max )
  {
    #segm1 <- rad1**2 * ( 2 * psi1 - sin(2 * psi1) ) / 2
    segm1 <- A.a**2 * (2 * zeta1 - sin(2 * zeta1)) / 2
#     segm1 <- za**2 * (2 * zeta1 - sin(2 * zeta1)) / 2
    
    #segm2 <- rad2**2 * ( 2 * psi2 - sin(2 * psi2) ) / 2
    segm2 <- B.a**2 * (2 * zeta2 -sin(2 * zeta2)) / 2
#     segm2 <- zb**2 * (2 * zeta2 -sin(2 * zeta2)) / 2
    
    #semisquare <- pi * rad2**2 / 2 
    semisquare <- pi * B.a**2 / 2
#     semisquare <- pi * zb**2 / 2
    
    #obtuse_arc <- rad2**2 * ( 2 * psi2 ) / 2
    obtuse_arc <- B.a**2 * (2 * zeta2) / 2
#     obtuse_arc <- zb**2 * (2 * zeta2) / 2
    
    #triangle <- rad2**2 * sin( -2 * psi2) / 2
    triangle <- B.a**2 * sin( -2 * zeta2 ) / 2
#     triangle <- zb**2 * sin( -2 * zeta2 ) / 2
    
    #if ( psi2 - pi / 2 < eps )
    if (abs(zeta2 - pi /2) <= eps)
    {
      return( s - light_down * ( semisquare + segm1 ) )
    }
    if ( zeta2 - pi / 2 < eps)
    {
      return( s - light_down * (segm1 + segm2) )
    }		   
    
    #if ( abs( psi2 - pi / 2 ) <= eps )
    
    
    return( s - light_down * ( segm1 + obtuse_arc + triangle ) )
  }
  
  return(0)
}

# Funcion: ellmodel()
# This function calls a square model
#
# Input parameters:
# t -time (phase of a period T)
#
# Returns:
# modelled square at the "t" moment
ellmodel <- function(t)
{
  #s <- lux1 * pi * rad1**2 + lux2 * pi * rad2**2
  s <- A.lux * pi * A.big.axis * A.small.axis +
      B.lux * pi * B.big.axis * B.small.axis
  #val <- -0.37 + apmag + 2.5*log10(s/v_circle_square(t, incl))
  val <- ZERO.POINT +  app.mag + 2.5*log10(s/v_ellipse_square(t, ORBIT.INCL))
  #val <- 12.95198 +  app.mag + 2.5*log10(s/v_ellipse_square(t, ORBIT.INCL))
  return(val)
}


# Function: cmodel()
# Circle modeling function
#
# Input paramteres:
# t -time (phase of a period T)
cmodel <- function(t)
{
  rad1 <- A.big.axis
  rad2 <- B.big.axis
  
  lux1 <- A.lux
  lux2 <- B.lux
  
  s <- lux1 * pi * rad1**2 + lux2 * pi * rad2**2
#   s <- A.lux * pi * A.big.axis * A.small.axis +
#     B.lux * pi * B.big.axis * B.small.axis
  #val <- -0.37 + apmag + 2.5*log10(s/v_circle_square(t, incl))
  val <- ZERO.POINT +  app.mag + 2.5*log10(s/v_circle_square(t, ORBIT.INCL))
  #val <- 12.95198 +  app.mag + 2.5*log10(s/v_ellipse_square(t, ORBIT.INCL))
  return(val)
}

# SMOOTHED DATA
# column 1: time
# column 2: model data
# column 3: smooth data
smoothd <- c()

# Function: readSmoothData()
# This function reads smoothed data.
# Returns: a data table of smoothed data
readSmoothData <- function()
{
  if (is.null(smoothd))
  {
    obs <-  read.table("obs/pmm_all_1610155_15")
    
    res <- c()
    for (q in 1:nrow(obs))
      res <- rbind(res, c(obs[q,1],
                          ellmodel(obs[q,1]),
                          obs[q,2]
      ))
    write.table(res,"smoothd.data", col.names=F, row.names=T)
    return(res)
  }
  return(NULL)
}
smoothd  <- readSmoothData()

# RAW DATA
# col 1: index number
# col 2: folded JD - phase
# col 3: JD
# col 4: C data
# col 5: O data
rawd <-c()

# This function reads raw data
# and returns a data fold
readRawData <- function()
{
  w <- read.table("obs/all_10_13")
  #w <- w[order(w[1,]),]
  
  sp <- c()
  for (t in 1:nrow(w))
  {
    sp <- rbind(sp, c(
      t,                          # number
      w[t, 1] %% ORBIT.P,         # folded JD - phase
      w[t, 1],                    # JD
      ellmodel(w[t, 1] %% ORBIT.P),  # C data
      w[t, 2]                     # O data
    ))
  }
  write.table(sp,"rawd.data", row.names=F, col.names=F)
  #return(sp)
  return(sp[order(sp[,3]),])
}
rawd <- readRawData()




# PLOT FUNCTIONS

# Function: plot model
plot.model <- function(save=F)
{
  par(mfrow=c(1,1))
  # plot Smooth vs Model
  plot(smoothd[,1], smoothd[,2], col=rgb(0,0,0), pch = 18, cex=0.5,
       main="GU CMa light curve",
       ylab="apparent magnitude (m)",
       xlab="period (days)",
       ylim=c(6.5,6.18))
  points(smoothd[,1], smoothd[,3], col=rgb(0,0,1), pch =18, cex=0.5)
  #text(0.35,7.5,col=rgb(0,0,0),"model")
  #text(0.35,7.7, col=rgb(0,0,1),"data")
  #legend("bottomleft", c("model", "data"), pch=c(18, 18), cex=0.5,
  #       col=c(rgb(0,0,0), rgb(0,0,1)) )
  #text(1.1,6.70, "parameters", adj=0)
  #text(1.1,6.75, paste("inclination", ORBIT.INCL * 180 / pi), adj=0)
  #text(1.1,6.80, paste("radii", A.big.axis, ";", B.big.axis), adj=0)
  #text(1.1,6.85, paste("separation", ORBIT.R), adj=0)
  #text(1.1,6.90, paste("lux factors", A.lux, ";", B.lux), adj=0)
  #text(1.1,6.95, paste("period", ORBIT.P), adj=0)
  #text(1.1,7.0, paste("initial phase", ORBIT.PHI0 * 180 / pi), adj=0)
  
  # save to pdf
  if (save)
  {
    filename <- paste("plot.model",
                      "incl", ORBIT.INCL * 180 / pi,
                      "major",A.big.axis, B.big.axis,
                      "ellipticity", round(A.small.axis / A.big.axis, 4),
                      round(B.small.axis / B.big.axis, 4), 
                      "lux",A.lux,B.lux,
                      "phase",round(ORBIT.PHI0, 4),
                      sep="_")
    
    # save this plot to pdf
    dev.copy2pdf(
      file = paste0("plots/",filename,".pdf"),
      width=16, height=8)
    
    # save to eps
    dev.copy2eps(
      file=paste0("plots/",filename,".eps"),
      width=16, height=8)
  }
}
#plot.model(save=T)


#Function: period.fold
# Paramters:
# indexe -subset of indexes
# period to fold
period.fold <- function(indexes, p, save=F)
{
  res <- c()
  # folding
  for (t in indexes)
  {
    res <- rbind(res,c(
      rawd[t,3],              #JD
      rawd[t,3] %% p,         # folded JD
      rawd[t,5]-rawd[t,4]     # O-C
      ))
  }
  
  #return(sp[order(sp[,1]),])
  res <- res[order(res[,2]),]
  #write.table(res)
  
  
  #fres <- filter(res, rep(1,10)/10)
  
  par(mfrow=c(1,1))
  xticks <- seq(0,max(res[,2]),by=0.01)
  yticks <- seq(-0.1,max(res[,3]),by=0.01)
  axis(1, at=xticks)
  axis(2, at=yticks)
  
  
  
  plot(res[,2], res[,3],
       main=paste("Period fold P=",p),
       xlab="time (days)",
       ylab="app. mag. (m)",
       #type="l",
       pch=18,
       cex=0.6,
       panel.first=grid(length(xticks), length(yticks), equilogs=T)
       )
 
  if (save)
  {
    filename <- paste("period.fold",
                      "indexes", indexes,
                      "p",p,
                      sep="_")
    
    dev.copy2pdf(
      file = paste0("plots/",filename,".pdf"),
      width=16, height=8)
    
    # save to eps
    dev.copy2eps(
      file=paste0("plots/",filename,".eps"),
      width=16, height=8)
  }
  return(invisible(res))
}



# Function: plot.curve
# This functions plots smoothed data foled 1.61 period
plot.curve <- function(save=T)
{
  par(mfrow=c(1,1))
  plot(smoothd[,1], smoothd[,3], col=rgb(0,0,1), pch = 18, cex=0.5,
       main="GU CMa light curve",
       ylab="apparent magnitude (m)",
       xlab="period (days)",
       ylim=c(6.55,6.15))
  
  if (save)
  {
    filename <- paste("curve_folded_1.6",
                      "incl", ORBIT.INCL * 180 / pi,
                      "major",A.big.axis, B.big.axis,
                      "ellipticity", round(A.small.axis / A.big.axis, 4),
                      round(B.small.axis / B.big.axis, 4), 
                      "lux",A.lux,B.lux,
                      "phase",round(ORBIT.PHI0, 4),
                      sep="_")
    
    dev.copy2pdf(
      file = paste0("plots/",filename,".pdf"),
      width=16, height=8)
    
    # save to eps
    dev.copy2eps(
      file=paste0("plots/",filename,".eps"),
      width=16, height=8)
  }
}


draw <- function()
{
  
  # open new device 
 # dev.new()
  
  # prepare filename for pdf
  filename <- paste("plot",
                   "incl", ORBIT.INCL * 180 / pi,
                   "major",A.big.axis, B.big.axis,
                   "ellipticity", round(A.small.axis / A.big.axis, 4),
                                  round(B.small.axis / B.big.axis, 4), 
                   "lux",A.lux,B.lux,
                   "phase",round(ORBIT.PHI0, 4),
                   sep="_")
  
  # set plot parameters for 4 graphs on a page
par(mfrow=c(3,2))
#   par(mfrow=c(2,1))
  
  ### plot.model function body
  # plot Smooth vs Model
  plot(smoothd[,1], smoothd[,2], col=rgb(0,0,0), pch = 2, cex=0.5,
       main="GU CMa light curve",
       ylab="apparent magnitude (m)",
       xlab="time (days)",
       ylim=c(6.55,6.15))
  points(smoothd[,1], smoothd[,3], col=rgb(0,0,1), pch =1, cex=0.5)
  #text(0.35,7.5,col=rgb(0,0,0),"model")
  #text(0.35,7.7, col=rgb(0,0,1),"data")
  #legend("bottomleft", c("model", "data"), pch=c(2, 1), cex=0.5,
   #      col=c(rgb(0,0,0), rgb(0,0,1)) )
  #text(1.1,6.70, "parameters", adj=0)
  #text(1.1,6.75, paste("inclination", ORBIT.INCL * 180 / pi), adj=0)
  #text(1.1,6.80, paste("radii", A.big.axis, ";", B.big.axis), adj=0)
  #text(1.1,6.85, paste("separation", ORBIT.R), adj=0)
  #text(1.1,6.90, paste("lux factors", A.lux, ";", B.lux), adj=0)
  #text(1.1,6.95, paste("period", ORBIT.P), adj=0)
  #text(1.1,7.0, paste("initial phase", ORBIT.PHI0 * 180 / pi), adj=0)

 
  # plot Residuals Smooth vs Model
  plot(smoothd[,1],smoothd[,3]-smoothd[,2], pch=".", type="l",
       main="Residuals (smoothed data)",
       ylab="apparent magnitude (m)",
       xlab="time (days)")
  
  
 # dev.copy2pdf(file = "lc.pdf",width=16, height=16)
 # stopifnot(1==2)

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
       xlab="JD - JD2000",)
  
  
  # fit
#   vect <- rawd[1:N,5] - rawd[1:N,4]
#   t <- rawd[1:N,3]
#   lm.s <- lm(vect ~ sin(t) + cos(t))
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
  
  # Plot Lomb-Scargle periodogram 2010 data
  indexes <- 1:832
  GM.dat <- rawd[indexes,5] - rawd[indexes,4]
  GM.time <- rawd[indexes,3]
  GM.ts <- ts(GM.dat, GM.time)
  
  GM.lsfreq <- lsp(GM.dat,
                        time = GM.time,
                        type="frequency",
                        #from=0.1,
                        to=15,
  )
  text(20,120, paste("peaks at:",
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
       xlab="JD - JD2000",)
  
  
  GM2013.lsfreq <- lsp(GM2013.dat,
                   time = GM2013.time,
                   type="frequency",
                   #from=0.1,
                   #to=50,
  )
  text(8,60, paste("peaks at:",
                     round(GM2013.lsfreq$peak.at[1],4),
                     round(GM2013.lsfreq$peak.at[2],4)))
  
  
  # save this plot to pdf
  dev.copy2pdf(
    file = paste0("plots/",filename,".pdf"),
    width=24, height=24)
  
  # save to eps
  dev.copy2eps(
    file=paste0("plots/",filename,".eps"),
    width=24, height=24)
}

draw()
# regression <- function()
# {
#   # Create time series
#   indexes <- 1:832
#   GM.dat <- rawd[indexes,5] - rawd[indexes,4]
#   GM.time <- rawd[indexes,3]
#   GM.ts <- ts(GM.dat, GM.time)
#   
# 
#   sin.t <- sin(2*pi*GM.time)
#   cos.t <- cos(2*pi*GM.time)
#   lm.GM <- lm(GM.dat ~ sin.t + cos.t)
#   plot(lm.GM)
#   
#   #plot.ts(GM.ts)
#   #acf(GM.ts)
# }
#regression()


# lombscargle <- function()
# {
#   # Create time series
#   indexes <- 1:832
#   GM.dat <- rawd[indexes,5] - rawd[indexes,4]
#   GM.time <- rawd[indexes,3]
#   GM.ts <- ts(GM.dat, GM.time)
#   
#   GM.lombscargle <- lsp(GM.dat,
#                         time = GM.time,
#                         type="frequency",
#                         #from=0.1,
#                         #to=50,
#                         )
#   
#   return(GM.lombscargle)
# }
#lslist <- lombscargle()


# Function: lcmodel(phases)
#
# This function returns a light curve model vector 
# for a given argument vector
# Input: phase of period vector
# Result: light curve vector according to the model
lcmodel <- function(phases, plot=FALSE)
{ 
  v <- c()
  for (q in phases) 
    v <- rbind(v, cmodel(q))
  
  if (plot)
  {
    par(mfrow=c(1,1))
    plot(phases, v,
         main="Model",
         xlab="phases of a period",
         ylab="app. mag. (m)",
         ylim=c(max(v)+0.15,min(v)-0.15),
         #pch=1,
         type="l"
         )
    return(invisible(v))
  }
  return(v)
}
# starmodel <- lcmodel(res[,1], plot=T)



# # Function: minimize
# # This function will fit model to raw data and will find residuals
lcresid<- function(z)
{
  return(sum((res[,3]-z+lcmodel(res[,1]))^2)/nrow(res)^2)
}


solve.mag <- function()
{
  # mag constant
  mc <- 15
  
  # given m(a+b)
  m.ab <- 6.93
  
  # fluxes factor
  k <- (B.lux*B.big.axis^2)/(A.lux*A.big.axis^2)
  
  # flux a
  f.a <- 1/(1+k)*10^((mc- m.ab)/2.5)
  
  m.a <- -2.5*log10(f.a) + mc
  m.b <- m.a + 2.5*log10(1/k)
  return(c(m.a, m.b))
}



analyze.2010 <- function()
{
  ts.2010 <- zoo(rawd[i.2010,5]-rawd[i.2010,4], rawd[i.2010,3])
  
  min.interval <- min(diff(rawd[i.2010,3]))
  even.2010 <- merge(ts.2010,
                     zoo(, seq(start(ts.2010), end(ts.2010), by=min.interval)),
                     all=TRUE)
  even.2010[is.na(even.2010)] <- 0
  even.2010[even.2010!=0] <- 1
  
  # analyze uneven data
  par(mfrow=c(2,1))
  ts.lsp <- lsp(coredata(ts.2010), time=index(ts.2010), type="frequency",
                from=0.0, to=15, ofac=1)
  
  # analyze white noise
  even.lsp <- lsp(coredata(even.2010), time=index(even.2010), type="frequency",
                 from=0.0, to=15, ofac=1)
  dev.copy2pdf(
    file = "plots/lspg_2010.pdf",
    width=16, height=8)
  
  
}

analyze.2013 <- function()
{
  ts.2013 <- zoo(rawd[i.2013,5]-rawd[i.2013,4], rawd[i.2013,3])
  
  min.interval.2013 <- min(diff(rawd[i.2013,3]))
  even.2013 <- merge(ts.2013,
                     zoo(, seq(start(ts.2013), end(ts.2013), by=min.interval)),
                     all=TRUE)
  even.2013[is.na(even.2013)] <- 0
  even.2013[even.2013!=0] <- 1
  
  # analyze uneven data
  par(mfrow=c(2,1))
  ts.lsp.2013 <- lsp(coredata(ts.2013), time=index(ts.2013), type="frequency",
                from=0.0, to=15, ofac=1)
  
  # analyze white noise
  even.lsp.2013 <- lsp(coredata(even.2013), time=index(even.2013), type="frequency",
                  from=0.0, to=15, ofac=1)
  dev.copy2pdf(
    file = "plots/lspg_2013.pdf",
    width=16, height=8)
  
  
}

analyze.all <- function()
{
  ts.all <- zoo(rawd[,5]-rawd[,4], rawd[,3])
  
  min.interval <- min(diff(rawd[,3]))
  even.all <- merge(ts.all,
                     zoo(, seq(start(ts.all), end(ts.all), by=min.interval)),
                     all=TRUE)
  even.all[is.na(even.all)] <- 0
  even.all[even.all!=0] <- 1
  
  # analyze uneven data
  par(mfrow=c(2,1))
  ts.lsp.all <- lsp(coredata(ts.all), time=index(ts.all), type="frequency",
                     from=0.0, to=15, ofac=1)
  
  # analyze white noise
  even.lsp.all <- lsp(coredata(even.all), time=index(even.all), type="frequency",
                       from=0.0, to=15, ofac=1)
  dev.copy2pdf(
    file = "plots/lspg_all.pdf",
    width=16, height=8)
  
  
}