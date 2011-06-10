setClass('Caline4Site',
	representation(
		U = 'numeric',
		BRG = 'numeric',
		CLAS = 'numeric',
		MIXH = 'numeric',
		SIGTH = 'numeric',
		TEMP = 'numeric',
		Z0 = 'numeric'
	)
)

setMethod('initialize', 'Caline4Site', function(.Object, site.data) {
	.Object@U <- as.single(with(site.data, wind.speed))
	.Object@BRG <- as.single(with(site.data, wind.bearing))
	.Object@CLAS <- as.single(with(site.data, stability.class))
	.Object@SIGTH <- as.single(with(site.data, sigma.theta))
	.Object@TEMP <- as.single(with(site.data, temperature))
	.Object@Z0 <- as.single(with(site.data, surface.roughness))
	.Object
})
