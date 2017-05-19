options = struct('evalCode',false,'showCode',false,'maxHeight',400,'maxWidth',400);

publish('ERK_FWD_Integrator.m',options);
publish('RK_FWD_Integrator.m',options);
publish('ROS_FWD_Integrator.m',options);
publish('SDIRK_FWD_Integrator.m',options);