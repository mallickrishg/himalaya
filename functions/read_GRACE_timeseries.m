function [tvec,dm,longr,latgr] = read_GRACE_timeseries(list,datafolder)
% read_GRACE_timeseries to read GRACE grid files
% INPUTS - 
% list - directory list from dir()
% datafolder - name of datafolder
% OUTPUTS -
% tvec - cell-array of (julian day) time vectors for each GRACE point
% dm - cell-array of mass change vectors (GigaTon)
% longr,latgr - (double-array) lon/lat of GRACE data
% Rishav Mallick EOS,2018

for i = 1:length(list)
    fin = [datafolder '/' list(i).name];
    fid = fopen(fin);
    A = textscan(fid,'%f %f','HeaderLines',46,'Delimiter',',');
    fclose(fid);
    fid = fopen(fin);
    B = textscan(fid,'%s %s',1,'HeaderLines',6,'Delimiter',',');
    fclose(fid);
    fid = fopen(fin);
    C = textscan(fid,'%s %s',1,'HeaderLines',5,'Delimiter',',');
    fclose(fid);
    % define variables to be used
    t = A{1};
    % mass change in Gt
    dm{i} = A{2};
    % date vector in julian days
    tvec{i} = round(ConvertSerialYearToDate(t));
    latgr(i) = round(str2double(cell2mat(C{2})),1);
    longr(i) = round(str2double(cell2mat(B{2})),1);
end