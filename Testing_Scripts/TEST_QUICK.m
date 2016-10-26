options = MatlOde_OPTIONS('Jacobian',@vanDerPol_Jacobian,'AbsTol', ones(2,1)*1d-5, 'RelTol', ones(2,1)*1d-5 );
tspan = [ 0 20 ];
y0 = [2;-0.66];
load TEST_QUICK_vdp % lambda_vdp_erk_adj lambda_vdp_rk_adj lambda_vdp_ros_adj lambda_vdp_sdirk_adj

options_ref = odeset('Jacobian',@vanDerPol_Jacobian);

[ Tref, Yref ] = ode15s(@vanDerPol_Function,tspan,y0,options_ref);

% FWD
[ Terk, Yerk ] = MatlOde_ERK_FWD_DRIVER_Integrator(@vanDerPol_Function,tspan,y0,options);
[ Trk, Yrk ] = MatlOde_RK_FWD_DRIVER_Integrator(@vanDerPol_Function,tspan,y0,options);
[ Tros, Yros ] = MatlOde_ROS_FWD_DRIVER_Integrator(@vanDerPol_Function,tspan,y0,options); %...still relies on LS_SOLVER
[ Tsdirk, Ysdirk ] = MatlOde_SDIRK_FWD_DRIVER_Integrator(@vanDerPol_Function,tspan,y0,options);

% TLM

options = MatlOde_OPTIONS( options, 'Lambda', eye(2), 'AbsTol_ADJ', ones(2,1)*1d-4, 'RelTol_ADJ', ones(2,1)*1d-4 );
% ADJ
[ tErkAdj, yErkAdj, sensErkAdj ] = MatlOde_ERK_ADJ_DRIVER_Integrator(@vanDerPol_Function,tspan,y0,options);

error_erk = norm(Yerk(:,end)-Yref(end,:)')/norm(Yref(end,:)');
error_rk = norm(Yrk(:,end)-Yref(end,:)')/norm(Yref(end,:)');
error_ros = norm(Yros-Yref(end,:))/norm(Yref(end,:)');
error_sdirk = norm(Ysdirk-Yref(end,:))/norm(Yref(end,:)');

error_erk_adj_lambda = norm(sensErkAdj-lambda_vdp_erk_adj)/norm(lambda_vdp_erk_adj);

output_erk   = [ 'FWD ERK   relative error: ' num2str(error_erk) ];
output_rk    = [ 'FWD RK    relative error: ' num2str(error_rk) ];
output_ros   = [ 'FWD ROS   relative error: ' num2str(error_ros) ];
output_sdirk = [ 'FWD SDIRK relative error: ' num2str(error_sdirk) ];

output_erk_adj_lambda = [ 'ADJ ERK Lambda relative error: ' num2str(error_erk_adj_lambda) ];

disp(output_erk);
disp(output_rk);
disp(output_ros);
disp(output_sdirk);

disp(output_erk_adj_lambda);