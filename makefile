# test: gk.f90 midpoint.f90 main.f90
# 	ifort -o test gk.f90 midpoint.f90 main.f90 
test: gk2.f90 main2.f90
	ifort -o test gk2.f90 main2.f90 -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib