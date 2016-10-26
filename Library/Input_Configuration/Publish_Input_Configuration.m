%% Publish_Input_Configuration
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
options = struct('showCode',false,'evalCode',false,'maxHeight',500,'maxWidth',500);

publish('OPTIONS_Compare.m',options);
publish('OPTIONS_GeneralConfiguration.m',options);
publish('OPTIONS_Merge.m',options);

options = struct('showCode',true,'evalCode',false,'maxHeight',500,'maxWidth',500);
publish('OPTIONS_Configuration.m',options);