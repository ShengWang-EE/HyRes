function gas_demand_table = preprocess_UK_gas_demand_data(filename)
%PREPROCESS_UK_GAS_DEMAND_DATA
% Reorganize NTS daily off-take table into: 1 row per calendar day (by Applicable At date),
% 1 column per Data Item (typically 5 items).
% Rules:
%   - Grouping day = dateshift(Applicable At, 'start','day') (ignore HH:MM)
%   - For each day, keep ONLY the latest Applicable At "batch" (e.g., 05:16 beats 00:00)
%   - Within that batch, if duplicates for same (Day, Data Item) exist, keep the latest row
%   - Missing values are filled by linear interpolation over days (EndValues='nearest')
%   - Output which days/items were filled

%% ============ 0) Read ============

%% 0) Read (preserve original headers)
T = readtable(filename, 'VariableNamingRule','preserve');

%% 1) Extract columns
appAt = T.("Applicable At");
item  = string(T.("Data Item"));
val   = T.Value;

% Parse datetime if needed
if ~isdatetime(appAt)
    appAt = datetime(appAt, 'InputFormat','dd/MM/yyyy HH:mm', 'Locale','en_GB');
end

% Numeric value
if ~isnumeric(val)
    val = str2double(string(val));
end
val = double(val);

% Calendar day based on Applicable At
day = dateshift(appAt, 'start','day');

TT = table(day, appAt, item, val, 'VariableNames', {'Day','AppAt','Item','Value'});

% Drop rows missing essentials
TT = TT(~isnat(TT.AppAt) & strlength(TT.Item) > 0, :);

%% 2) Keep latest record for each (Day, Item)
% Sort so latest AppAt is first within each (Day, Item)
TT = sortrows(TT, {'Day','Item','AppAt'}, {'ascend','ascend','descend'});

% Unique by (Day, Item) => keep first (latest AppAt)
[~, ia] = unique(TT(:, {'Day','Item'}), 'rows', 'stable');
TT = TT(ia, :);

%% 3) (Optional) restrict to the 5 expected items, in a fixed order
items5 = [ ...
    "NTS Energy Offtaken, Interconnector Exports Total", ...
    "NTS Energy Offtaken, Industrial Offtake Total", ...
    "NTS Energy Offtaken, LDZ Offtake Total", ...
    "NTS Energy Offtaken, Powerstations Total", ...
    "NTS Energy Offtaken, Storage Injection Total" ...
];

TT = TT(ismember(TT.Item, items5), :);

%% 4) Pivot to wide (1 row per day)
W = unstack(TT, 'Value', 'Item','GroupingVariables', 'Day');   % columns will be MATLAB-ized variable names
W = sortrows(W, 'Day');

%% 5) Ensure all 5 item columns exist (even if missing in data)
% Map original item strings -> the actual variable names created by unstack
% The created names are in W.Properties.VariableNames (MATLAB made them valid).
wantedVars = matlab.lang.makeValidName(items5);
for k = 1:numel(wantedVars)
    v = wantedVars(k);
    if ~ismember(v, W.Properties.VariableNames)
        W.(v) = NaN(height(W), 1);
    end
end

% Reorder columns: Day + the 5 items
W = W(:, ["Day", wantedVars]);

%% 6) Log missing + fill
numCols = wantedVars;

missMask = false(height(W), numel(numCols));
for i = 1:numel(numCols)
    c = numCols(i);
    x = W.(c);
    missMask(:,i) = isnan(x);
    W.(c) = fillmissing(x, 'linear', 'EndValues','nearest');
end

T_daily = W;

%% 7) Output filled days & log
rowsFixed = any(missMask, 2);
filledDays = T_daily.Day(rowsFixed);

[rowIdx, colIdx] = find(missMask);
filledLog = table( ...
    T_daily.Day(rowIdx), ...
    string(numCols(colIdx))', ...
    'VariableNames', {'Day','DataItemVar'} );




gas_demand_table = T_daily;

end