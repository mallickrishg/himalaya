function [ox,d,Cd,flag] = create_inputdataset(filename,datatype)

datainput = readmatrix(filename);
ox = datainput(:,1).*1e3; % locations of data - convert to [m]

d = datainput(:,2);
Cd = diag(datainput(:,3).^2);

if strcmp(datatype, 'horizontal') == 1
    flag = zeros(length(d),1);
elseif strcmp(datatype,'vertical') == 1
    flag = ones(length(d),1);
else
    error('unknown data type: provide as horizontal OR vertical')
end

end