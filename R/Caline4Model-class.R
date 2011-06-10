setClass('Caline4Model',
	representation(
		receptors = 'Caline4Receptors', 
		links = 'Caline4Links',
		site = 'Caline4Site',
		pollutant = 'Caline4Pollutant'
	)
)

setMethod('initialize', 'Caline4Model', 
	function(.Object, receptors, links, site, pollutant) {
		if(!inherits(receptors, 'Caline4Receptors')) receptors <- new('Caline4Receptors', receptors)
		if(!inherits(links, 'Caline4Links')) links <- new('Caline4Links', links)
		if(!inherits(site, 'Caline4Site')) site <- new('Caline4Site', site)
		if(!inherits(pollutant, 'Caline4Pollutant')) pollutant <- new('Caline4Pollutant', pollutant)
		.Object@receptors <- receptors
		.Object@links <- links
		.Object@site <- site
		.Object@pollutant <- pollutant
		.Object
	}
)

predict.Caline4Model <- function(mod) {
	CALINE4.Fortran(
		mod@receptors@NR, mod@receptors@XR, mod@receptors@YR, mod@receptors@ZR,
		mod@links@NL, mod@links@XL1, mod@links@YL1, mod@links@XL2, mod@links@YL2, mod@links@WL, mod@links@HL, mod@links@TYP, mod@links@VPHL, mod@links@EFL,
		mod@site@U, mod@site@BRG, mod@site@CLAS, mod@site@MIXH, mod@site@SIGTH, mod@site@TEMP, mod@site@Z0, 
		mod@pollutant@PTYP, mod@pollutant@MOWT, mod@pollutant@VS, mod@pollutant@VD
	)
}