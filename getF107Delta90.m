function T = getF107Delta90(startDate, endDate)
%GETF107DELTA90  One % uncertainty per month using Δ90 = (H95-L5)/(2F)*100
%   T = getF107Delta90()
%   T = getF107Delta90(startDate, endDate)
%
% Returns a table with:
%   Month        : datetime (yyyy-MM)
%   F            : predicted smoothed monthly F10.7
%   L5, H95      : 5th and 95th percentile bounds
%   Delta90_pct  : ((H95 - L5) / (2*F)) * 100    <-- single % per month
%
% 
% Author: Rafail Panagiotidis
% The University of Manchester
% August 2025
%
%--- Copyright notice ---%
% Copyright (C) 2025 The University of Manchester

    % Inputs & timezone (Month is UTC; match filters to UTC to avoid warnings)
    if nargin < 1 || isempty(startDate), startDate = datetime(1900,1,1); end
    if nargin < 2 || isempty(endDate),   endDate   = datetime(9999,12,1); end
    startDate = dateshift(startDate,'start','month'); startDate.TimeZone = 'UTC';
    endDate   = dateshift(endDate,  'start','month'); endDate.TimeZone   = 'UTC';

    % Fetch JSON
    url = 'https://services.swpc.noaa.gov/json/solar-cycle/predicted-solar-cycle.json';
    S   = webread(url, weboptions('Timeout',30));

    % --- Field extraction (handles dash/underscore variants) ---
    tags = getCellField(S, {'time_tag','time-tag'});
    F    = getNumField (S, {'predicted_f10_7','predicted_f10_7'}); % keep both for safety
    L5   = getNumField (S, {'low_f10_7','low_f10_7','low_f10_7'}); % older MATLAB sometimes repeats
    H95  = getNumField (S, {'high_f10_7','high_f10_7','high_f10_7'});

    % Parse months
    Month = datetime(tags, 'InputFormat','yyyy-MM', 'TimeZone','UTC');

    % Compute single % per month
    Delta90_pct = ((H95 - L5) ./ (2 .* F)) * 100;

    % Assemble, filter, sort
    T = table(Month, F, L5, H95, Delta90_pct);
    mask = (T.Month >= startDate) & (T.Month <= endDate);
    T = sortrows(T(mask,:), 'Month');
end

% === helpers ===
function v = getNumField(S, candidates)
    for k = 1:numel(candidates)
        if isfield(S, candidates{k})
            v = [S.(candidates{k})]';
            return
        end
    end
    error('Expected numeric field not found. Tried: %s', strjoin(candidates, ', '));
end

function c = getCellField(S, candidates)
    for k = 1:numel(candidates)
        if isfield(S, candidates{k})
            c = {S.(candidates{k})}';
            return
        end
    end
    error('Expected string/cell field not found. Tried: %s', strjoin(candidates, ', '));
end

