function [tsroi] = read_TRMM_timeseries(list,datafolder,longr,latgr)
% read_TRMM_timeseries will read surfacePrecipCOMB for level3-B31 TRMM
% files stored in a directory
% INPUT -
% list - structure created by using dir() containing names of all the files
% datafolder - char-array containing all the data
% longr,latgr - (double array) locations of points where you want toe
% xtract the TRMM timeseries
% Rishav Mallick EOS,2018

%hardcoded definition of gridpoints (see TRMM_README file)
mlon=-179.75:0.5:180;
mlat=-39.75:0.5:40;
[TRlon,TRlat] = meshgrid(mlon,mlat);

tic
for j = 1:length(longr)
    tic
    for i = 1:length(list)
        % read combined surface precipitation from hygrometer and radar sensor
        din = hdfread([datafolder '/' list(i).name],'surfacePrecipCOMB');
        instruct = xml2struct([datafolder '/' list(i).name '.xml']);
        % data needs to be transposed to be in the correct shape
        D = din';
        % get start and end date for each data point
        tst = instruct.S4PAGranuleMetaDataFile.RangeDateTime.RangeBeginningDate.Text;
        tend = instruct.S4PAGranuleMetaDataFile.RangeDateTime.RangeEndingDate.Text;
        
        [r,~] = distance(latgr(j),longr(j),TRlat(:),TRlon(:));
        rind = find(r==min(r),1,'first');
        
        %define a new structure that contains lat,lon, dates and precip values
        %for the points of interest
        tsroi(j).data(i,1) = D(rind);        
        tsroi(j).tvec(i,1) = round(mean(datenum([tst;tend]))); % middle of each month
    end
    tsroi(j).lonTRMM = TRlon(rind);
    tsroi(j).latTRMM = TRlat(rind);
    tsroi(j).lonGR = longr(j);
    tsroi(j).latGR = latgr(j);
    disp(['Finished processing station ' num2str(j)])
    toc
end
disp('Finished extracting TRMM for all stations')
toc
end