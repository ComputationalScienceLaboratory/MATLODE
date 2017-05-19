%% Publish_SingleStepMethod
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
options = struct('showCode',false,'evalCode',false,'maxHeight',500,'maxWidth',500);

publish('ArnoldiAdapt.m',options);
publish('erow4AdaptMatFree.m',options);
publish('exp4AdaptMatFree.m',options);
publish('exp4SingleStep.m',options);