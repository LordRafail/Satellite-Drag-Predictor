function T = getF107DailyUncertainty(method, first7)
%GETF107DAILYUNCERTAINTY  SWPC 45-day F10.7 + per-day % uncertainty
%   T = getF107DailyUncertainty()
%   T = getF107DailyUncertainty(method)               % 'pchip'|'linear'|'previous'|'nearest'
%   T = getF107DailyUncertainty(method, first7)       % 1x7 vector, days 1..7 (%)
%
% Output:
%   T with columns: Date (UTC), DayAhead (1..N), F107, Uncertainty_pct
%
% Notes:
%   • Days 1..7 use your given %: [4.2 6.1 7.8 9.3 10.5 11.7 12.7]
%   • Days 8..45 are extended via interp1(..., method, 'extrap').
%   • Robust to one-line SWPT format like:
%     "45-DAY F10.7 CM FLUX FORECAST 12Aug25 150 13Aug25 145 ... FORECASTER: ...".
%
% Source format confirmed here (single-line with AP + F10.7 blocks). :contentReference[oaicite:0]{index=0}

    if nargin < 1 || isempty(method), method = 'pchip'; end
    if nargin < 2 || isempty(first7), first7 = [4.2 6.1 7.8 9.3 10.5 11.7 12.7]; end
    validateattributes(first7, {'numeric'}, {'real','finite','vector','numel',7,'nonnegative'});

    % Pre-allocate output so it's always assigned
    T = table(datetime.empty(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
              'VariableNames', {'Date','DayAhead','F107','Uncertainty_pct'});

    % 1) Fetch file
    url = 'https://services.swpc.noaa.gov/text/45-day-ap-forecast.txt';
    raw = webread(url, weboptions('Timeout',30));     % returns char on most MATLAB versions
    if isstring(raw), raw = char(raw); end

    % 2) Slice the F10.7 block
    iStart = regexpi(raw, '45-?DAY\s+F10\.?7\s*CM\s*FLUX\s*FORECAST', 'once');
    if isempty(iStart)
        warning('Could not find the "45-DAY F10.7 CM FLUX FORECAST" header. Returning empty table.');
        return
    end
    tail = raw(iStart:end);
    iStop = regexpi(tail, '(FORECASTER:|NNNN|#)', 'once');
    if ~isempty(iStop), tail = tail(1:iStop-1); end

    % 3) Parse pairs like "12Aug25 150"
    toks = regexp(tail, '(\d{2}[A-Za-z]{3}\d{2})\s+(-?\d+(?:\.\d+)?)', 'tokens');
    if isempty(toks)
        warning('Found the header but no date/value pairs. Returning empty table.');
        return
    end

    n = numel(toks);
    if n > 45, toks = toks(1:45); n = 45; end

    % 4) Convert to arrays
    dates = NaT(n,1); dates.TimeZone = 'UTC';
    vals  = zeros(n,1);
    for k = 1:n
        dates(k) = datetime(toks{k}{1}, 'InputFormat','ddMMMyy', 'TimeZone','UTC');
        vals(k)  = str2double(toks{k}{2});
    end

    % 5) Build uncertainty vector: given 1..7, extend to n
    d  = (1:n).';
    yK = first7(:);
    try
        Unc = interp1(1:7, yK, d, method, 'extrap');
    catch
        Unc = interp1(1:7, yK, d, 'pchip', 'extrap');
    end
    Unc = max(Unc, 0);

    % 6) Assemble output
    T = table(dates, d, vals, Unc, 'VariableNames', {'Date','DayAhead','F107','Uncertainty_pct'});
end