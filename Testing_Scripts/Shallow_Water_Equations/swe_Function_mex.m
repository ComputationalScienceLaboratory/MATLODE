function [ yp ] = swe_Function_mex( t, y )

    [NX, NY, dx, dy] = swe_Constant_Parameters();
    yp = swe_fun_mex(y,NX+2,NY+2,dx,dy);    
    %yp = sweCartesian_fun_mex(t,y);    
end

