setClass('Caline4Receptors',
	representation(
		NR = 'integer',
		XR = 'numeric',
		YR = 'numeric',
		ZR = 'numeric'
	)
)

setMethod('initialize', 'Caline4Receptors', function(.Object, receptor.data) {
	.Object@NR <- as.integer(nrow(receptor.data))
	.Object@XR <- as.single(with(receptor.data, x))
	.Object@YR <- as.single(with(receptor.data, y))
	.Object@ZR <- as.single(with(receptor.data, z))
	.Object
})