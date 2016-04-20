function [ H, L, U, p, ISTATUS ] = ROK_PrepareMatrix( M, H, Direction, gam, Harn, ISTATUS, OPTIONS )

    Nconsecutive = 0;
    ISING = true;

    while ( ISING )
        % Construct Ghimj = 1/(H*ham) - Jac0
        gh = (Direction*H*gam);

        % Compute LU decomposition
%        [ ISING, e ] = lss_decomp( NVAR, ghinv, fjac );
        [L, U, p, singular] = lu_decomposition(speye(M) - gh * Harn);
        ISTATUS.Ndec = ISTATUS.Ndec + 1;

        if ( ~singular ) % If successful, done.
            ISING = false;
        else % If unsuccessful halve the step size; If 5 consecutive fails then return.
            ISTATUS.Nsng = ISTATUS.Nsng + 1;
            Nconsecutive = Nconsecutive + 1;
            ISING = true;
            str = [ 'Warning: LU Decomposition returned ISING = ', num2str(ISING) ];
            disp(str);
            if ( Nconsecutive <= 5 ) % Less than 5 consecutive failed decompositions
                H = H*0.5;
                H = max( OPTIONS.Hmin, min( H, OPTIONS.Hmax ) );
            else
                % More than 5 consecutive failed decompositions
                error('Matrix is repeatedly singular.')
                return;
            end % Nconsecutive
        end % ising
    end % while singular

return;

