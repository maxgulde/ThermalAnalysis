% reads in csv file
function [data] = ReadCSV(path,fstring,headerlines)

% check if file exists
if (exist(path,'file') == 0)
    fprintf(['Error: Did not find file <' path '>.']);
else
    fid = fopen(path);
    data = textscan(fid,fstring, 'HeaderLines', headerlines);
    fclose(fid);    
end
end