MEX=/Applications/MATLAB_R2015b.app/bin/mex
#MEX_LINUX=/opt/MATLAB/R2013b/bin/mex
mac:  swe_jac_mex.c swe_fun_mex.c
	$(MEX) -O swe_jac_mex.c -output swe_jac_mex; $(MEX) -O swe_fun_mex.c -output swe_fun_mex

#linux:	swe_jac_mex.c swe_fun_mex.c
#	$(MEX_LINUX) -O swe_jac_mex.c -output swe_jac_mex;	$(MEX_LINUX) -O swe_fun_mex.c -output swe_fun_mex

clean:
	rm *.mex*




