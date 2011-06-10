CALINE4.Fortran <- function(
	NR, XR, YR, ZR,
	NL, XL1, YL1, XL2, YL2, WL, HL, TYP, VPHL, EFL,
	U, BRG, CLAS, MIXH, SIGTH, TEMP, Z0, 
	PTYP, MOWT, VS, VD
) {	
	array.shape <- c(NR)
	returned_data <- .Fortran(
		"CALINE4", 
		NR, XR, YR, ZR,
		NL, XL1, YL1, XL2, YL2, WL, HL, TYP, VPHL, EFL,
		U, BRG, CLAS, MIXH, SIGTH, TEMP, Z0, 
		PTYP, MOWT, VS, VD,
		C = as.single(array(0.0, dim=array.shape)),
		PACKAGE = "CALINE4"
	)
	predicted <- returned_data$C
	dim(predicted) <- array.shape
	return(predicted)	
}