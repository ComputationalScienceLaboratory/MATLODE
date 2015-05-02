function PrintISTATUS( ISTATUS )

    Nfun_formatSpec = '    Nfun: %6.0f \n';
    Njac_formatSpec = '    Njac: %6.0f \n';
    Nstp_formatSpec = '    Nstp: %6.0f \n';
    Nacc_formatSpec = '    Nacc: %6.0f \n';
    Nrej_formatSpec = '    Nrej: %6.0f \n';
    Ndec_formatSpec = '    Ndec: %6.0f \n';
    Nsol_formatSpec = '    Nsol: %6.0f \n';
    Nsng_formatSpec = '    Nsng: %6.0f \n';
%    Nchk_formatSpec = '    Nchk: %6.0f \n';
    
    fprintf( '\n\nISTATUS = \n\n' );
    fprintf( Nfun_formatSpec, ISTATUS.Nfun );
    fprintf( Njac_formatSpec, ISTATUS.Njac );
    fprintf( Nstp_formatSpec, ISTATUS.Nstp );
    fprintf( Nacc_formatSpec, ISTATUS.Nacc );
    fprintf( Nrej_formatSpec, ISTATUS.Nrej );
    fprintf( Ndec_formatSpec, ISTATUS.Ndec );
    fprintf( Nsol_formatSpec, ISTATUS.Nsol );
    fprintf( Nsng_formatSpec, ISTATUS.Nsng );
%    fprintf( Nchk_formatSpec, ISTATUS.Nchk );

return;

