setClass('Caline4Links',
	representation(
		NL = 'integer',
		XL1 = 'numeric',
		XL2 = 'numeric',
		YL1 = 'numeric',
		YL2 = 'numeric',
		WL = 'numeric',
		HL = 'numeric',
		TYP = 'integer',
		VPHL = 'numeric',
		EFL = 'numeric'
	)
)

setMethod('initialize', 'Caline4Links', function(.Object, link.data) {
	.Object@NL <- as.integer(nrow(link.data))
	.Object@XL1 <- as.single(with(link.data, x1))
	.Object@YL1 <- as.single(with(link.data, y1))
	.Object@XL2 <- as.single(with(link.data, x2))
	.Object@YL2 <- as.single(with(link.data, y2))
	.Object@WL <- as.single(with(link.data, width))
	.Object@HL <- as.single(with(link.data, height))
	.Object@TYP <- as.integer(with(link.data, type))
	.Object@VPHL <- as.single(with(link.data, vehicles.per.hour))
	.Object@EFL <- as.single(with(link.data, emission.factor))
	.Object
})
