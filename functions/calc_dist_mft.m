function [loc_flt,dmft] = calc_dist_mft(infil,gpslon,gpslat)
% Calculate distance of gps points from the MFT
% The MFT is broken up into a series of line segments and i find the
% minimum distance from each line segment
% OUTPUTS  
% loc_flt - the location on the MFT which is closest to the gps station
% dmft - distance from loc_flt to gps station in km
% Rishav Mallick, EOS, 2018
mftin = readtable(infil);

[loc_flt,~,~] = unicycle.utils.distance2curve(mftin{:,1:2},[gpslon,gpslat],'linear');
[alen,az] = distance(gpslat,gpslon,loc_flt(:,2),loc_flt(:,1));
dsgn = az.*0;

% sign convention for station-mft distance
IF = az>=270 | az<=90;
dsgn(IF) = -1;
dsgn(~IF) = 1;
dmft = dsgn.*deg2km(alen);

end