.packageName <- 'CALINE4'
.onLoad <- function(libname, pkgname) {
	library.dynam(.packageName, pkgname, libname)
	.packageVersion <- utils::pkgnameDescription(.packageName, field="Version"),
	pkgnameStartupMessage("This is ", .packageName, " ", .packageVersion, appendLF = TRUE)
}
.onUnload <- function(libname) library.dynam.unload(.packageName, libname)
.onAttach <- function(libname, pkgname) {}
