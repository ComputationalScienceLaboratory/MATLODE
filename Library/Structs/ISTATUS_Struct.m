%% ISTATUS_Struct
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
function options = ISTATUS_Struct(varargin)

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
  fprintf( ' Nfun: [ No. of function calls ] \n' );
  fprintf( ' Njac: [ No. of Jacobian calls ]\n' );
  fprintf( ' Nstp: [ No. of steps ]\n' );
  fprintf( ' Nacc: [ No. of accepted steps ]\n' );
  fprintf( ' Nrej: [ No. of rejected steps (except at very beginning) ]\n' );
  fprintf( ' Ndec: [ No. of LU decompositions ]\n' );
  fprintf( ' Nsol: [ No. of forward/backward substitutions ]\n' );
  fprintf( ' Nsng: [ No. of singular matrix decompositions]\n' );
  fprintf( ' Nchk: [ No. of allocated memory chunks ]\n' );
  fprintf( ' Nkdim: [ No. of Krylov vectors]\n' );
  fprintf('\n');
  return;
end

if ( nargin == 1 && strcmp('default',varargin) )
    options.Nfun = 0;
    options.Njac = 0;
    options.Nstp = 0;
    options.Nacc = 0;
    options.Nrej = 0;
    options.Ndec = 0;
    options.Nsol = 0;
    options.Nsng = 0;
    options.Nchk = 1;
    options.Nkdim = 0;
    return;
end

% Names = [
%     'Nfun ' % ISTATUS(1)
%     'Njac ' % ISTATUS(2)
%     'Nstp ' % ISTATUS(3)
%     'Nacc ' % ISTATUS(4)
%     'Nrej ' % ISTATUS(5)
%     'Ndec ' % ISTATUS(6)
%     'Nsol ' % ISTATUS(7)
%     'Nsng ' % ISTATUS(8)
%     'Nchk '
%     ];
% m = size(Names,1);
% names = lower(Names);
% 
% % Combine all leading options structures o1, o2, ... in odeset(o1,o2,...).
% options = [];
% for j = 1:m
%   options.(deblank(Names(j,:))) = [];
% end
% i = 1;
% while i <= nargin
%   arg = varargin{i};
%   if ischar(arg)                         % arg is an option name
%     break;
%   end
%   if ~isempty(arg)                      % [] is a valid options argument
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
%         options.(deblank(Names(j,:))) = val;
%       end
%     end
%   end
%   i = i + 1;
% end
% 
% % A finite state machine to parse name-value pairs or default option.
% if rem(nargin-i+1,2) ~= 0
%     if ( nargin == 1 && strcmp('default',varargin) )
%         options.Nfun = 0;
%         options.Njac = 0;
%         options.Nstp = 0;
%         options.Nacc = 0;
%         options.Nrej = 0;
%         options.Ndec = 0;
%         options.Nsol = 0;
%         options.Nsng = 0;
%         options.Nchk = 1;
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
%     options.(deblank(Names(j,:))) = arg;
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
