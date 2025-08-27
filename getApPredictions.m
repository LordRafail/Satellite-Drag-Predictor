function [ap_out, usedBinUTC, source_out] = getApPredictions(targetDates)
%GETAPPREDICTIONS Exact Ap with strict precedence and nearest-of-two-bin rule.
%
%   [ap, usedBinUTC, source] = getApPredictions(targetDates)
%
%   For times between two 3-hour bins (e.g., 22:27), we compare the two
%   candidate bin starts (floor and ceil) and take the closer one (ties -> floor).
%
%   STRICT PRECEDENCE (per query, no averaging):
%     1) GFZ historical 3-hour (prefer definitive D=1 on duplicates)
%        https://www-app3.gfz-potsdam.de/kp_index/Kp_ap_since_1932.txt
%     2) ESA 72-hour forecast (downloaded JSON file, 3-hour bins)
%     3) SWPC 27-day outlook (daily Ap) → flat across that UTC day
%        https://services.swpc.noaa.gov/text/27-day-outlook.txt
%     4) NASA long-range monthly Ap (flat within/nearest month, persist forward)
%        https://www.nasa.gov/wp-content/uploads/2025/08/aug2025ap-smt.txt
%
%   OUTPUTS
%     ap_out     : Ap used for each query time
%     usedBinUTC : the exact 3-hour bin timestamp actually used
%     source_out : string array: "GFZ", "ESA72", "SWPC27", or "NASA-monthly"
%
%   Notes
%     • Everything is in UTC.
%     • Network fetches are time-bounded to avoid hanging.

  if ~isa(targetDates,'datetime'), error('targetDates must be datetime'); end
  qTimes = toUTC(targetDates(:));

  % ---- Load sources ----
  [gT, gV] = loadGFZ_3h();        % GFZ observed, D-aware
  [eT, eV] = loadESAJson72h();    % ESA JSON 72h forecast (local file)
  [dD, dA] = loadNOAA27day();     % SWPC 27-day (daily Ap)
  [nMonDates, nMonAp] = loadNASAapMonthly(); % NASA monthly Ap

  % ---- Build maps ----
  gMap = make3hMap(gT, gV);
  eMap = make3hMap(eT, eV);
  dMap = makeDayMap(dD, dA);

  % ---- Per-query ----
  n  = numel(qTimes);
  ap = nan(n,1);
  used = NaT(n,1); used.TimeZone = 'UTC';
  src  = strings(n,1);

  nowUTC = datetime('now','TimeZone','UTC');
  warnSingular = false;

  for k = 1:n
      qt = qTimes(k);
      f  = floor3h(qt);
      c  = f + hours(3);

      % 1) GFZ
      [v, tSel] = pickNearestOfTwo(qt, f, c, gMap);
      if ~isnan(v), ap(k)=v; used(k)=tSel; src(k)="GFZ";  continue; end

      % 2) ESA JSON 72h
      [v, tSel] = pickNearestOfTwo(qt, f, c, eMap);
      if ~isnan(v), ap(k)=v; used(k)=tSel; src(k)="ESA72"; continue; end

      % 3) SWPC 27-day (flat per day)
      dk = dayKey(qt);
      if isKey(dMap, dk)
          ap(k)   = dMap(dk);
          used(k) = f; src(k)="SWPC27";
          if hours(qt - nowUTC) > 72, warnSingular = true; end
          continue
      end

      % 4) NASA monthly (persist forward)
      if ~isempty(nMonDates)
          qm = dateshift(qt,'start','month');
          idx = find(year(nMonDates)==year(qm) & month(nMonDates)==month(qm),1);
          if isempty(idx)
              idx = find(nMonDates <= qm,1,'last');
              if isempty(idx)
                  idx = find(nMonDates >= qm,1,'first');
              end
          end
          if ~isempty(idx) && isfinite(nMonAp(idx))
              ap(k)   = nMonAp(idx);
              used(k) = f; src(k)="NASA-monthly";
              if hours(qt - nowUTC) > 72, warnSingular = true; end
              continue
          end
      end
  end

  if warnSingular
      warning('Using singular daily/monthly Ap beyond 72h — no diurnal structure.');
  end

  ap_out     = ap;
  usedBinUTC = used;
  source_out = src;
end

%% ====================== Helpers ======================
function tUTC = toUTC(t)
  if isempty(t)
      tUTC = NaT(size(t)); tUTC.TimeZone = 'UTC'; return
  end
  if ~isdatetime(t), t = datetime(t); end
  t.TimeZone = 'UTC';
  tUTC = t;
end

function t3 = floor3h(t)
  t = toUTC(t);
  t3 = dateshift(t,'start','hour') - hours(mod(hour(t),3));
end

function dk = dayKey(t)
  dk = int32(year(t)*10000 + month(t)*100 + day(t));
end

function M = make3hMap(times, vals)
  M = containers.Map('KeyType','int64','ValueType','double');
  if isempty(times) || isempty(vals), return, end
  times = toUTC(times);
  [times, ia] = unique(times,'last'); vals = vals(ia);
  for i = 1:numel(times)
      if ~isnat(times(i)) && isfinite(vals(i))
          M(int64(posixtime(times(i)))) = vals(i);
      end
  end
end

function M = makeDayMap(dayDates, dayVals)
  M = containers.Map('KeyType','int32','ValueType','double');
  if isempty(dayDates) || isempty(dayVals), return, end
  dayDates = dateshift(toUTC(dayDates),'start','day');
  [dayDates, ia] = unique(dayDates,'last'); dayVals = dayVals(ia);
  for i=1:numel(dayDates)
      M(int32(dayKey(dayDates(i)))) = dayVals(i);
  end
end

function [v, tSel] = pickNearestOfTwo(qt, f, c, map3h)
  v   = nan; tSel = NaT; tSel.TimeZone = 'UTC';
  if isempty(map3h), return, end
  kf = int64(posixtime(f));
  kc = int64(posixtime(c));
  hasF = isKey(map3h, kf);
  hasC = isKey(map3h, kc);
  if ~hasF && ~hasC, return
  elseif hasF && ~hasC, v=map3h(kf); tSel=f;
  elseif ~hasF && hasC, v=map3h(kc); tSel=c;
  else
      df = seconds(abs(qt - f)); dc = seconds(abs(c - qt));
      if dc < df, v=map3h(kc); tSel=c;
      else, v=map3h(kf); tSel=f; end
  end
end

%% ====================== Loaders ======================

function [times, vals] = loadGFZ_3h()
  url = 'https://www-app3.gfz-potsdam.de/kp_index/Kp_ap_since_1932.txt';
  fn  = 'Kp_ap_since_1932.txt';
  needFetch = ~isfile(fn);
  if ~needFetch
      d = dir(fn);
      needFetch = isempty(d) || (now - d.datenum) > (12/24);
  end
  if needFetch
      try, websave(fn, url, weboptions('Timeout',20)); catch, end
  end
  times = NaT(0,1); times.TimeZone='UTC'; vals = [];
  fid = fopen(fn,'r');
  if fid<0, return, end
  C = textscan(fid,'%f %f %f %f %f %f %f %s %f %f',...
               'CommentStyle','#','MultipleDelimsAsOne',true);
  fclose(fid);
  if any(cellfun(@isempty,C)), return, end
  Y=C{1}; Mo=C{2}; D=C{3}; HH=C{4}; ap=C{9}; Def=C{10};
  t0=datetime(Y,Mo,D,HH,0,0,'TimeZone','UTC');
  good=isfinite(ap)&ap~=-1&isfinite(Def);
  T=table(t0(good),ap(good),Def(good),...
          'VariableNames',{'t0','ap','Def'});
  T=sortrows(T,{'t0','Def'});
  [~,ia]=unique(T.t0,'last');
  times=T.t0(ia); vals=T.ap(ia);
  [times,s]=sort(times); vals=vals(s);
end

function [times, vals] = loadESAJson72h(jsonFile)
  if nargin<1, jsonFile='apforecast72hr_download.json'; end
  times=NaT(0,1); times.TimeZone='UTC'; vals=[];
  if ~isfile(jsonFile), return, end
  raw=fileread(jsonFile); J=jsondecode(raw);
  if ~isfield(J,'fields'), return, end
  F=J.fields; n=numel(F);
  T=NaT(n,1); T.TimeZone='UTC'; V=nan(n,1);
  for i=1:n
      try
          yr=F(i).year; doy=F(i).doy; int=F(i).int; ap=double(F(i).ap);
          d0=datetime(yr,1,1,0,0,0,'TimeZone','UTC')+days(doy-1);
          binH=(int-1)*3;
          T(i)=d0+hours(binH); V(i)=ap;
      catch, end
  end
  good=isfinite(V)&~isnat(T);
  T=T(good); V=V(good);
  [T,ia]=unique(T,'last'); V=V(ia);
  [times,s]=sort(T); vals=V(s);
end

function [dailyDates, dailyA] = loadNOAA27day()
  url='https://services.swpc.noaa.gov/text/27-day-outlook.txt';
  dailyDates=NaT(0,1); dailyDates.TimeZone='UTC'; dailyA=[];
  try
      txt=webread(url,weboptions('Timeout',20,'ContentType','text'));
      if isstring(txt), txt=char(txt); end
  catch, return, end
  rowPat='(?<y>\d{4})\s+(?<mon>[A-Za-z]{3})\s+(?<d>\d{1,2})\s+(?<f>\d+)\s+(?<a>\d+)\s+(?<kp>\d+)';
  rows=regexp(txt,rowPat,'names');
  if isempty(rows), return, end
  n=numel(rows); dd=NaT(n,1); dd.TimeZone='UTC'; aa=nan(n,1);
  for i=1:n
      y=str2double(rows(i).y);
      m=month(datetime(['01-' rows(i).mon '-2000'],'InputFormat','dd-MMM-yyyy'));
      d=str2double(rows(i).d);
      dd(i)=datetime(y,m,d,0,0,0,'TimeZone','UTC');
      aa(i)=str2double(rows(i).a);
  end
  [dailyDates,s]=sort(dd); dailyA=aa(s);
end

function [monDates, monAp] = loadNASAapMonthly()
  url='https://www.nasa.gov/wp-content/uploads/2025/08/aug2025ap-smt.txt';
  monDates=NaT(0,1); monDates.TimeZone='UTC'; monAp=[];
  try
      txt=webread(url,weboptions('Timeout',20,'ContentType','text'));
      if isstring(txt), txt=char(txt); end
  catch, return, end
  m=regexp(txt,'(?<yf>\d{4}\.\d+)\s+(?<mean>\d+\.\d+)','names');
  if isempty(m), return, end
  n=numel(m); yf=zeros(n,1); mn=zeros(n,1);
  for i=1:n, yf(i)=str2double(m(i).yf); mn(i)=str2double(m(i).mean); end
  Y=floor(yf); mon=floor((yf-Y)*12)+1; mon=max(1,min(12,mon));
  monDates=datetime(Y,mon,1,0,0,0,'TimeZone','UTC'); monAp=mn;
end
