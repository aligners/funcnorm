function displayLogItem(str, logFile, incTime)  
% FUNCTION displayLogItem(str)
%	Displays str on screen.  Also appends the current time to the end of the line.
%
% FUNCTION displayLogItem(str, logFile [, incTime=1])
%	Displays str on screen, as well as to the specified log file
%	incTime is an optional argument that will append the current time to the end of the line.
%	If logFile=='', then the function will only display str on screen.
%
%
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

    if nargin == 2
        if length(logFile)
            writeLogFile = 1;
        else
            writeLogFile = 0;
        end
        
        incTime = 1;
    elseif nargin == 1
        writeLogFile = 0;
        incTime = 1;
    else
        writeLogFile = 0;
    end

    if incTime
        str = [str, '  **(', datestr(now), ')**'];
    end

    if writeLogFile
        fp = fopen(logFile, 'a');
        fprintf(fp, '%s\n', str);
        fclose(fp);
    end

    % Also write to the screen
    disp(str);
    pause(1e-6);
