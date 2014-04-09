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
lux2 <- 0.6

# incl -inclination angle
incl <- 70 * pi / 180


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
	phi <- 2 * pi / P * t

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


run <- function()
{
  obs <-  read.table("obs/pmm_all_1610155_15")
  
  s <- lux1 * pi * rad1**2 + lux2 * pi * rad2**2
  
  res <- c()
  for (q in 1:nrow(obs))
    res <- rbind(res, c(obs[q,1], -0.38 +
                        apmag + 2.5*log10(s/v_circle_square(obs[q,1], incl)),
                        obs[q,2]
                         ))
  write.table(res,"model.dat", col.names=F, row.names=F)
  return(res)
}

res <- run()

draw <- function()
{
  
  # plot for file
  dev.new()
  filename <- paste("plot",
                   "incl", incl * 180 / pi,
                   "radii",rad1,rad2,
                   "sep",r,
                   "lux",lux1,lux2,
                   "P",P,
                   
                   sep="_")
  par(mfrow=c(1,2))
  
  # plot O,C
  plot(res[,1], res[,2], col=rgb(0,0,0), pch = 18,
       main=" рива€ блеска GU CMa",
       ylab="видима€ зв. величина",
       xlab="период в дол€х (дни)",
       ylim=c(7,6))
  points(res[,1], res[,3], col=rgb(0,0,1), pch =18)
  #text(0.35,7.5,col=rgb(0,0,0),"model")
  #text(0.35,7.7, col=rgb(0,0,1),"data")
  legend("bottomleft", c("model", "data"), pch=c(18, 18),
         col=c(rgb(0,0,0), rgb(0,0,1)) )
  text(1.1,6.70, "параметры модели", adj=0)
  text(1.1,6.75, paste("наклон", incl * 180 / pi), adj=0)
  text(1.1,6.80, paste("радиусы", rad1, ";", rad2), adj=0)
  text(1.1,6.85, paste("separation", r), adj=0)
  text(1.1,6.90, paste("светимости", lux1, ";", lux2), adj=0)
  text(1.1,6.95, paste("период", P), adj=0)
  
  # plot O-C
  plot(res[,1],abs(res[,2]-res[,3]), pch=18,
       main="—пектр нев€зок",
       ylab="величина нев€зки",
       xlab="период в дол€х (дни)")
  #pdf(paste0("plots/",filename), width = 4, height = 4)
  dev.copy2pdf(
    file = paste0("plots/",filename,".pdf"),
    width=16, height=8)
  #dev.off()
}

draw()