%% RSTATUS_Struct
%
% <html>
%   <div>
%       <img style="float: right" src="../../../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
%% Syntax
%
%% Input Parameters
%
%% Output Parameters
%
%% Description
%
%% Reference
% [1] Tony D'Augustine, Adrian Sandu. MATLODE: A MATLAB ODE Solver and
%     Sensitivity Analysis Toolbox. Submitted to ACM TOMS.
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function RSTATUS = RSTATUS_Struct(varargin)

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    fprintf( ' Ntexit: [ Time corresponding to the computed Y upon return ]\n' );
    fprintf( '  Nhacc: [ Last accepted step before exit ]\n' );
    fprintf( '  Nhnew: [ Last predicted step (not yet taken) ]\n' );
    fprintf( ' Nhexit: [ Step size at exit ]\n' );
    fprintf('\n');
    return;
end

if ( nargin == 1 && strcmp('default',varargin) )
    RSTATUS.Ntexit = 0;
    RSTATUS.Nhacc = 0;
    RSTATUS.Nhnew = 0;
    RSTATUS.Nhexit = 0;
    RSTATUS.Etime = 0;
    return;
end

% Names = [
%     'Ntexit ' % RSTATUS(1) [ RK_ADJ_MATLAB_Integrator ]
%     'Nhacc  ' % RSTATUS(2) [ RK_ADJ_MATLAB_Integrator ]
%     'Nhnew  ' % RSTATUS(3) [ RK_ADJ_MATLAB_Integrator ]
%     'Nhexit '
%     ];
% m = size(Names,1);
% names = lower(Names);
% 
% % Combine all leading RSTATUS structures o1, o2, ... in odeset(o1,o2,...).
% RSTATUS = [];
% for j = 1:m
%   RSTATUS.(deblank(Names(j,:))) = [];
% end
% i = 1;
% while i <= nargin
%   arg = varargin{i};
%   if ischar(arg)                         % arg is an option name
%     break;
%   end
%   if ~isempty(arg)                      % [] is a valid RSTATUS argument
%     if ~isa(arg,'struct')
%       error(message('MATLAB:odeset:NoPropNameOrStruct', i));
%     end
%     for j = 1:m
%       if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
%         val = arg.(deblank(Names(j,:)));
%       else
%         val = [];
%       end
%       if ~isempty(val)
%         RSTATUS.(deblank(Names(j,:))) = val;
%       end
%     end
%   end
%   i = i + 1;
% end
% 
% % A finite state machine to parse name-value pairs.
% if rem(nargin-i+1,2) ~= 0
%     if ( nargin == 1 && strcmp('default',varargin) )
%         RSTATUS.Ntexit = 0;
%         RSTATUS.Nhacc = 0;
%         RSTATUS.Nhnew = 0;
%         RSTATUS.Nhexit = 0;
%         RSTATUS.Etime = 0;
%         return;
%     else
%         error(message('MATLAB:odeset:ArgNameValueMismatch'));
%     end
% end
% 
% expectval = 0;                          % start expecting a name, not a value
% while i <= nargin
%   arg = varargin{i};
%     
%   if ~expectval
%     if ~ischar(arg)
%       error(message('MATLAB:odeset:NoPropName', i));
%     end
%     
%     lowArg = lower(arg);
%     j = strmatch(lowArg,names);
%     if isempty(j)                       % if no matches
%       error(message('MATLAB:odeset:InvalidPropName', arg));
%     elseif length(j) > 1                % if more than one match
%       % Check for any exact matches (in case any names are subsets of others)
%       k = strmatch(lowArg,names,'exact');
%       if length(k) == 1
%         j = k;
%       else
%             matches = deblank(Names(j(1),:));
%         for k = j(2:length(j))'
%                 matches = [matches ', ' deblank(Names(k,:))]; %#ok<AGROW>
%         end
%             error(message('MATLAB:odeset:AmbiguousPropName',arg,matches));
%       end
%     end
%     expectval = 1;                      % we expect a value next
%     
%   else
%     RSTATUS.(deblank(Names(j,:))) = arg;
%     expectval = 0;
%       
%   end
%   i = i + 1;
% end
% 
% if expectval
%   error(message('MATLAB:odeset:NoValueForProp', arg));
% end


%% Major Modification History
% <html>
% <table border=1>
%   <tr>
%       <td><b>Date</b></td>
%       <td>Developer</td>
%       <td>Email</td>
%       <td>Action</td>
%   </tr>
%   <tr>
%       <td>1/1/2014</td>
%       <td>Tony D'Augustine</td>
%       <td>adaug13@vt.edu</td>
%       <td>Release MATLODE_v2.0.00</td>
%   </tr>
% </table>
% </html>
% 
%%
% <html>
%   <div>
%       <img style="float: right" src="../../../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html>