function par = get_xml_par(fnameXml,keyWord,momKeyWord)

%GET_RSC_PAR return parameter value of given keyword from .rsc file
%
%   [par,indx] = get_rsc_par(fnameRsc,keyWord) reads in .rsc file, searhes 
%   for specified keyword and returns the parameter value of the key word.
%   The *.rsc file is in standard format as in the ROI_PAC. If the second 
%   input parameter is omitted, all the keywords in the file are listed.
%
%INPUT:
%   fnameRsc    rsc filename
%   keyWord     keyword to find
%
%OUTPUT:
%   par         parameter value
%   indx        index number (i-th row) of the key word found

% Sang-Ho Yun, 2004.10.27 Original
% Sang-Ho Yun, 2012.02.02 Added string return

if (nargin == 0)
  help get_xml_par;
end

fid = fopen(fnameXml);
tline = fgetl(fid);
if exist('momKeyWord','var')
  while ischar(tline)
    ind = strfind(lower(tline),lower(momKeyWord));
    if ~isempty(ind)
      break;
    end
    tline = fgetl(fid);
  end
end
while ischar(tline)
  ind = strfind(lower(tline),lower(keyWord));
  if ~isempty(ind)
    tline = fgetl(fid);
    ind1 = strfind(tline,'<value>');
    ind2 = strfind(tline,'</value>');
    indBeg = ind1+7;
    indEnd = ind2-1;
    s = tline(indBeg:indEnd);
    par = str2num(s);
    fclose(fid);
    return;
  end
  tline = fgetl(fid);
end
par = [];
fclose(fid);
