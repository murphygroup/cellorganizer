% 2011-11-15 tebuck: channel_as_protein sets the protein image to
% that channel's image and blanks the channel's image. As a result,
% preventing division by zero with the if statements.


% $$$ tmp_idx = find((~prot_re)&(~nuc_re)&(~tub_re)&(~er_re));

tmp_idx = []; 
if ~(optimize && protein_channel_blank)
  if isempty(tmp_idx)
    tmp_idx = ~prot_re;
  else
    tmp_idx = tmp_idx & ~prot_re;
  end
end
if ~(optimize && nuclear_channel_blank)
  if isempty(tmp_idx)
    tmp_idx = ~nuc_re;
  else
    tmp_idx = tmp_idx & ~nuc_re;
  end
end
if ~(optimize && tubulin_channel_blank)
  if isempty(tmp_idx)
    tmp_idx = ~tub_re;
  else
    tmp_idx = tmp_idx & ~tub_re;
  end
end
if ~(optimize && er_channel_blank)
  if isempty(tmp_idx)
    tmp_idx = ~er_re;
  else
    tmp_idx = tmp_idx & ~er_re;
  end
end
tmp_idx = find(tmp_idx);


% $$$ prot_re(tmp_idx) = [];
% $$$ nuc_re(tmp_idx) = [];
% $$$ tub_re(tmp_idx) = [];
% $$$ er_re(tmp_idx) = [];


if ~(optimize && protein_channel_blank)
prot_re(tmp_idx) = [];
end
if ~(optimize && nuclear_channel_blank)
nuc_re(tmp_idx) = [];
end
if ~(optimize && tubulin_channel_blank)
tub_re(tmp_idx) = [];
end
if ~(optimize && er_channel_blank)
er_re(tmp_idx) = [];
end


mi = double(min(prot_re));    ma = double(max(prot_re));
if isempty(mi)
  mi = 0;
  ma = 0;
end
if ma - mi == 0
  ma = mi + 1;
end
prot_re8 = round( (double(prot_re)-mi) * 7/(ma-mi));
prot_re16 = round( (double(prot_re)-mi) * 15/(ma-mi));
prot_re32 = round( (double(prot_re)-mi) * 31/(ma-mi));
prot_re64 = round( (double(prot_re)-mi) * 63/(ma-mi));
prot_re128 = round( (double(prot_re)-mi) * 127/(ma-mi));

mi = double(min(nuc_re));    ma = double(max(nuc_re));
if isempty(mi)
  mi = 0;
  ma = 0;
end
if ma - mi == 0
  ma = mi + 1;
end
nuc_re8 = round( (double(nuc_re)-mi) * 7/(ma-mi));
nuc_re16 = round( (double(nuc_re)-mi) * 15/(ma-mi));
nuc_re32 = round( (double(nuc_re)-mi) * 31/(ma-mi));
nuc_re64 = round( (double(nuc_re)-mi) * 63/(ma-mi));
nuc_re128 = round( (double(nuc_re)-mi) * 127/(ma-mi));

mi = double(min(tub_re));    ma = double(max(tub_re));
if isempty(mi)
  mi = 0;
  ma = 0;
end
if ma - mi == 0
  ma = mi + 1;
end
tub_re8 = round( (double(tub_re)-mi) * 7/(ma-mi));
tub_re16 = round( (double(tub_re)-mi) * 15/(ma-mi));
tub_re32 = round( (double(tub_re)-mi) * 31/(ma-mi));
tub_re64 = round( (double(tub_re)-mi) * 63/(ma-mi));
tub_re128 = round( (double(tub_re)-mi) * 127/(ma-mi));

mi = double(min(er_re));    ma = double(max(er_re));
if isempty(mi)
  mi = 0;
  ma = 0;
end
if ma - mi == 0
  ma = mi + 1;
end
er_re8 = round( (double(er_re)-mi) * 7/(ma-mi));
er_re16 = round( (double(er_re)-mi) * 15/(ma-mi));
er_re32 = round( (double(er_re)-mi) * 31/(ma-mi));
er_re64 = round( (double(er_re)-mi) * 63/(ma-mi));
er_re128 = round( (double(er_re)-mi) * 127/(ma-mi));


if optimize && (protein_channel_blank || nuclear_channel_blank)
  feats1 = nan(1, 5);
else
feats1 = [mutualinformation( [prot_re8 nuc_re8]) ...
          mutualinformation( [prot_re16 nuc_re16]) ...
          mutualinformation( [prot_re32 nuc_re32]) ...
          mutualinformation( [prot_re64 nuc_re64]) ...
          mutualinformation( [prot_re128 nuc_re128])];
end

%keyboard

if optimize && (protein_channel_blank || tubulin_channel_blank)
  feats2 = nan(1, 5);
else
feats2 = [mutualinformation( [prot_re8 tub_re8]) ...
          mutualinformation( [prot_re16 tub_re16]) ...
          mutualinformation( [prot_re32 tub_re32]) ...
          mutualinformation( [prot_re64 tub_re64]) ...
          mutualinformation( [prot_re128 tub_re128])];
end

if optimize && (protein_channel_blank || er_channel_blank)
  feats3 = nan(1, 5);
else
feats3 = [mutualinformation( [prot_re8 er_re8]) ...
          mutualinformation( [prot_re16 er_re16]) ...
          mutualinformation( [prot_re32 er_re32]) ...
          mutualinformation( [prot_re64 er_re64]) ...
          mutualinformation( [prot_re128 er_re128])];
end

if optimize && (protein_channel_blank || nuclear_channel_blank)
  feats4 = nan(1, 5);
else
feats4 = [corr2( prot_re8, nuc_re8), ...
          corr2( prot_re16, nuc_re16), ...
          corr2( prot_re32, nuc_re32), ...
          corr2( prot_re64, nuc_re64), ...
          corr2( prot_re128, nuc_re128)];
end

if optimize && (protein_channel_blank || tubulin_channel_blank)
  feats5 = nan(1, 5);
else
feats5 = [corr2( prot_re8, tub_re8), ...
          corr2( prot_re16, tub_re16), ...
          corr2( prot_re32, tub_re32), ...
          corr2( prot_re64, tub_re64), ...
          corr2( prot_re128, tub_re128)];
end

if optimize && (protein_channel_blank || er_channel_blank)
  feats6 = nan(1, 5);
else
feats6 = [corr2( prot_re8, er_re8), ...
          corr2( prot_re16, er_re16), ...
          corr2( prot_re32, er_re32), ...
          corr2( prot_re64, er_re64), ...
          corr2( prot_re128, er_re128)];
end

feats = [feats1 feats2 feats3 feats4 feats5 feats6];

%{
names = {'nuc:mutual_information_8','nuc:correlation_8', ...
         'nuc:mutual_information_16','nuc:correlation_16', ...
         'nuc:mutual_information_32','nuc:correlation_32', ...
         'nuc:mutual_information_64','nuc:correlation_64', ...
         'nuc:mutual_information_128','nuc:correlation_128', ...
         'tub:mutual_information_8','tub:correlation_8', ...
         'tub:mutual_information_16','tub:correlation_16', ...
         'tub:mutual_information_32','tub:correlation_32', ...
         'tub:mutual_information_64','tub:correlation_64', ...
         'tub:mutual_information_128','tub:correlation_128', ...
         'er:mutual_information_8','er:correlation_8', ...
         'er:mutual_information_16','er:correlation_16', ...
         'er:mutual_information_32','er:correlation_32', ...
         'er:mutual_information_64','er:correlation_64', ...
         'er:mutual_information_128','er:correlation_128'};
%}
names = {'nuc:mutual_information_8','nuc:mutual_information_16','nuc:mutual_information_32','nuc:mutual_information_64','nuc:mutual_information_128',...
         'tub:mutual_information_8','tub:mutual_information_16','tub:mutual_information_32','tub:mutual_information_64','tub:mutual_information_128',...
         'er:mutual_information_8','er:mutual_information_16','er:mutual_information_32','er:mutual_information_64','er:mutual_information_128',...
         'nuc:correlation_8','nuc:correlation_16','nuc:correlation_32','nuc:correlation_64','nuc:correlation_128',...
         'tub:correlation_8','tub:correlation_16','tub:correlation_32','tub:correlation_64','tub:correlation_128',...
         'er:correlation_8','er:correlation_16','er:correlation_32','er:correlation_64','er:correlation_128'};


slfnames = repmat({''},[1 length(names)]);
