site.data <- data.frame(
	wind.speed = 1.0,
	wind.bearing = 270.0,
	stability.class = 6,
	mixing.height = 1000.0,
	sigma.theta = 15.0,
	temperature = 10.0,
	surface.roughness = 10.0)

link.data <- data.frame(
	x1 = 0, y1 = -5000, x2 = 0, y2 = 5000,
	type = 0,
	vehicles.per.hour = 7500,	# in	 vehicles per hour
	emission.factor = 30.0,		# in g/mile
	height = 0.0,
	width = 30.0)
row.names(link.data) <- c('HIGHWAY 22')

receptor.data <- data.frame(
	x = 30.0,
	y = 0.0,
	z = 1.8)
row.names(receptor.data) <- c('RESTSTOP')
	
receptors <- new("Caline4Receptors", receptor.data)
links <- new("Caline4Links", link.data)
site <- new("Caline4Site", site.data)
pollutant <- new("Caline4Pollutant", "CO")

model <- new("Caline4Model", receptors, links, site, pollutant)	