%% saveVariable
%
% <html>
%   <div>
%       <img style="float: right" src="../../MATLODE_LOGO.png" height="150px"></img>
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
function saveVariable( filePath, varargin )

    if ( nargin == 0 || mod(nargin-1,2) ~= 1 ) 
        disp( 'Error (saveVariable): Invalid input.' );
        return;
    end

    enableSave = varargin{nargin-1};
    if ( enableSave == true )
        i=1;
        while i <= ((nargin-1)/2)
            Data = varargin{2*i};
            eval( [ varargin{2*i-1} '= Data;' ]);
            save( [filePath, varargin{2*i-1}], varargin{2*i-1} );
            i = i + 1;
        end
    end
    
return;

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
%       <img style="float: right" src="../../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html>

