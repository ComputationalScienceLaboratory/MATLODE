function PrintRSTATUS( RSTATUS )

    Ntexit_formatSpec = '    Ntexit: %6.4f \n';
    Nhacc_formatSpec  = '     Nhacc: %6.4f \n';
    Nhnew_formatSpec  = '     Nhnew: %6.4f \n';
    Nhexit_formatSpec = '    Nhexit: %6.4f \n';
    Etime_formatSpec  = '     Etime: %6.4f \n';
    
    fprintf( '\n\nRSTATUS = \n\n' );
    fprintf( Ntexit_formatSpec, RSTATUS.Ntexit );
    fprintf( Nhacc_formatSpec, RSTATUS.Nhacc );
    fprintf( Nhnew_formatSpec, RSTATUS.Nhnew );
    fprintf( Etime_formatSpec, RSTATUS.Etime );
    
    if ( ~isempty(RSTATUS.Nhexit) )
        fprintf( Nhexit_formatSpec, RSTATUS.Nhexit );
    end

return;

