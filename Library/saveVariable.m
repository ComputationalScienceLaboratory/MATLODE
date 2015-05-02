function saveVariable( filePath, varargin )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: saveVariable.m
%
% Original Author: Tony D'Augustine
%
% File Create Date: 1 August 2012
%
% Input Arguments:
%   Name        Type
%   varargin    list
%
% Output Arguments:
%   No output
%
% Modification History:
%   Date        Developer       Email       Action
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% saveVariable:
%   Saves a single or list of variables into the current opened folder with
%   the ability to change the variable name. Analogous to saveas except now
%   one can save a list of variables.
%
% saveVariable: INPUT ARGUMENTS
%   varargin (list): variable-length input argument list. Please see
%                    SYNTAX section.
%   
% saveVariable: OUTPUT ARGUMENTS
%   No output
%
% saveVariable: CALLED FUNCTIONS
%   eval(...):
%   save(...):
%
% saveVariable: SYNTAX
%   saveVarialbe( 'varName1', var1, 'varName2', var2, ..., enableSave ); 
%
% saveVariable: EXAMPLE
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

