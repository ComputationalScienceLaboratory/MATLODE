%% Publish_Library
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
options = struct('showCode',false,'evalCode',false,'maxHeight',500,'maxWidth',500);

publish('errorNorm.m',options);
publish('fatOde_ErrorNorm_TLM.m',options);
publish('fatOde_ErrorScale.m',options);
publish('fatOde_ROS_ErrorNorm.m',options);
publish('fatOde_ROS_PrepareMatrix.m',options);
publish('PrintISTATUS.m',options);
publish('PrintRSTATUS.m',options);
publish('saveVariable.m',options);
publish('Library.m',options);