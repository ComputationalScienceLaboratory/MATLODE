%% Publish_UserSupplied_FWD
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
options = struct('evalCode',false,'showCode',false,'maxHeight',400,'maxWidth',400);

publish('OPTIONS_UserSupplied_ERK_FWD.m',options);
publish('OPTIONS_UserSupplied_EXP_FWD.m',options);
publish('OPTIONS_UserSupplied_EXPK_FWD.m',options);
publish('OPTIONS_UserSupplied_RK_FWD.m',options);
publish('OPTIONS_UserSupplied_ROK_FWD.m',options);
publish('OPTIONS_UserSupplied_ROS_FWD.m',options);
publish('OPTIONS_UserSupplied_SDIRK_FWD.m',options);