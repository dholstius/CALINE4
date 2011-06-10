setClass('Caline4Pollutant',
	representation(
		PTYP = 'integer',
		MOWT = 'numeric',
		VS = 'numeric',
		VD = 'numeric'
	)
)

setMethod('initialize', 'Caline4Pollutant', 
	function(.Object, name, molecular.weight, settling.velocity=0.0, deposition.velocity=0.0) {
		.Object@PTYP <- as.integer(switch(name, CO=1, NO2=2, INERT=3, PM=4, PARTICULATES=4))
		if(missing(molecular.weight)) {
			.Object@MOWT <- as.single(switch(name, CO=28., NO2=46.))
		} else {
			.Object@MOWT <- as.single(molecular.weight)
		}
		.Object@VS <- as.single(settling.velocity)
		.Object@VD <- as.single(deposition.velocity)
		.Object
	}
)

