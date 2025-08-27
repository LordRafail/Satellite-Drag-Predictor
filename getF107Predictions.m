function [f107_avg, f107_daily] = getF107Predictions(targetDates)
% F10.7 (previous-day) and centered 81-day avg using:
%   SWPC 45-day forecast (preferred when available) -> CelesTrak observations
% No monthly JSON. No linear interpolation.
%
% IMPORTANT: targetDates are treated as UTC midnights.

  if ~isa(targetDates,'datetime'), error('targetDates must be datetime'); end
  targetDates = dateshift(targetDates(:),'start','day'); 
  targetDates.TimeZone = 'UTC';
  n = numel(targetDates);

  f107_avg   = nan(n,1);
  f107_daily = nan(n,1);

  % ===== Observations (CelesTrak) =====
  histDates = NaT(0,1); histDates.TimeZone='UTC';
  histObs = []; histAvg = [];
  try
      csvPath = 'SW-Last5Years.csv';
      if ~isfile(csvPath)
          websave(csvPath,'https://celestrak.com/SpaceData/SW-Last5Years.csv');
      end
      opts = detectImportOptions(csvPath,'VariableNamingRule','preserve');
      T    = readtable(csvPath,opts);

      histDates = datetime(string(T{:,1}),'InputFormat','yyyy-MM-dd','TimeZone','UTC');

      obsCol = find(contains(T.Properties.VariableNames,'F10') & ...
                    contains(T.Properties.VariableNames,'OBS') & ...
                   ~contains(T.Properties.VariableNames,'CENTER81'),1,'first');
      avgCol = find(contains(T.Properties.VariableNames,'OBS_CENTER81'),1,'first');

      if ~isempty(obsCol), histObs = T{:,obsCol}; end
      if ~isempty(avgCol), histAvg = T{:,avgCol}; end
  catch
      warning('Could not load CelesTrak observations.');
  end

  % ===== SWPC 45-day F10.7 (daily steps) =====
  [d45, f45] = loadSWPC45DayF107();

  % ===== Outputs =====
  for k = 1:n
      d0 = targetDates(k);

      % daily = previous UTC day, preferring forecast for exact-day match
      prev = d0 - days(1);
      f107_daily(k) = pickF107(prev, histDates, histObs, d45, f45, true);

      % 81-day centered mean over ±40d, preferring forecast where available
      win = (d0-days(40)):(d0+days(40));
      F = nan(size(win));
      for j = 1:numel(win)
          F(j) = pickF107(win(j), histDates, histObs, d45, f45, true);
      end
      f107_avg(k) = mean(F,'omitnan');
  end
end

% ---- choose value with precedence: forecast exact-day -> obs exact-day -> forecast step/hold ----
function v = pickF107(qDate, histDates, histObs, d45, f45, preferForecastExact)
  % 1) If forecast contains THIS exact day, use it (matches website exactly)
  if preferForecastExact && ~isempty(d45)
      idxF = find(d45 == qDate, 1, 'first');
      if ~isempty(idxF), v = f45(idxF); return, end
  end

  % 2) Exact observed day
  if ~isempty(histDates)
      idxO = find(histDates == qDate, 1, 'first');
      if ~isempty(idxO), v = histObs(idxO); return, end
  end

  % 3) Forecast step/hold (last known day <= qDate)
  if ~isempty(d45)
      i = find(d45 <= qDate, 1, 'last');
      if ~isempty(i), v = f45(i); return, end
  end

  v = NaN;
end

% ---- loader for the F10.7 block from SWPC 45-day file ----
function [dates45, f45] = loadSWPC45DayF107()
  dates45 = NaT(0,1); dates45.TimeZone='UTC'; f45 = [];
  try
      raw = webread('https://services.swpc.noaa.gov/text/45-day-ap-forecast.txt', ...
                    weboptions('Timeout',30,'ContentType','text'));
      if isstring(raw), raw = char(raw); end
  catch
      return
  end

  iStart = regexpi(raw,'45-?DAY\s+F10\.?7\s*CM\s*FLUX\s*FORECAST','once');
  if isempty(iStart), return, end
  tail = raw(iStart:end);
  iStop = regexpi(tail,'(FORECASTER:|NNNN|#)','once');
  if ~isempty(iStop), tail = tail(1:iStop-1); end

  toks = regexp(tail,'(\d{2}[A-Za-z]{3}\d{2})\s+(-?\d+(?:\.\d+)?)','tokens');
  if isempty(toks), return, end

  m = numel(toks);
  d = NaT(m,1); d.TimeZone='UTC'; v = zeros(m,1);
  for i=1:m
      d(i) = datetime(toks{i}{1},'InputFormat','ddMMMyy', ...
                      'Locale','en_US','TimeZone','UTC');
      v(i) = str2double(toks{i}{2});
  end

  [d, s] = sort(d); v = v(s);
  [dates45, ia] = unique(d,'last'); f45 = v(ia);
end
