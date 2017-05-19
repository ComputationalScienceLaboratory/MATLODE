function TEST_MultipleRuns


    v_analysis = -1;
    v_error_mode = -2;
    v_integrator = [ -1 -3 -4 -5 -7 -8 -9 -10 -11 -12 ]
    v_adj_mode = -4;
    
    for i=1:length(v_analysis)
        for j=1:length(v_error_mode)
            for k=1:length(v_integrator)
                for l=1:length(v_adj_mode)
                    u_analysis = v_analysis(i);
                    u_error_mode = v_error_mode(j);
                    u_integrator = v_integrator(k);
                    u_adj_mode = v_adj_mode(l);
                    TEST_vanDerPol;
                end
            end
        end
    end

end

