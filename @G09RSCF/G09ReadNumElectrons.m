function numElectrons = G09ReadNumElectrons()
logFile = fopen('temp.log');
currLine = '';
while(ischar(currLine))
    if(~isempty(regexp(currLine, ' alpha electrons ', 'ONCE')))
        numCell = regexp(currLine, '[0-9]+', 'match');
        numElectrons = [str2double(numCell{1}), str2double(numCell{2})];
        break;
    end
    currLine = fgetl(logFile);
end
fclose(logFile);
end
