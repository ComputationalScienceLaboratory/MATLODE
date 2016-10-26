%% Publish_Memory
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
options = struct('showCode',false,'evalCode',false,'maxHeight',500,'maxWidth',500);

publish('matlOde_allocateMemory.m',options);
publish('matlOde_appendMemory.m',options);
publish('matlOde_deallocateMemory.m',options);