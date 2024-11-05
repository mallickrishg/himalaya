function rneu = GPS_readts(infil)
% function to read a GPS file downlaoded from Geoff Blewitt's website and
% extract data in the following format
% rneu = [jday decyr N E U sign sige sigz] all data in m
% Rishav Mallick, 2018, EOS

[~,yymmmdd,decyr,~,~,~,~,~,E,~,N,~,U,~,se,sn,su,~,~,~] = textread(infil,'%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);

t = datenum(yymmmdd,'yymmmdd');
rneu = [t decyr N E U sn se su];

end