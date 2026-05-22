function ZWD = excessphasezenith(prof, height)
% EXCESSPHASEZENITH  Compute Zenith Wet Delay from a wet refractivity profile
%
% Integrates a vertical wet refractivity profile (in N-units) over a
% height grid to produce the Zenith Wet Delay (ZWD) in metres using a
% midpoint-rule sum.
%
% CONVENTION: prof must contain one value per *layer*, not per boundary.
%   length(prof) == length(height) == numel(diff(model.levels_TOMO))
%
% Standard formula:
%   ZWD = 1e-6 * sum( Nw_i * dh_i )
% where Nw_i is the layer-average wet refractivity (sampled at midpoint)
% and dh_i is the layer thickness in metres.
%
% INPUT
%   prof   - [nLayers x 1] wet refractivity at layer midpoints (N-units)
%   height - [nLayers x 1] layer thicknesses in km
%            (= diff(model.levels_TOMO)./1000, so length = num_levels-1)
%
% OUTPUT
%   ZWD    - Zenith Wet Delay in metres
%
% See also: stackedNwcalc, intersetionNw, weightObs

% Convert layer thicknesses from km to metres
dh_m = height(:) * 1e3;   % [nH x 1], metres

% Ensure profile is a column vector of the same length
prof = prof(:);

if numel(prof) ~= numel(dh_m)
    error('excessphasezenith: profile length (%d) does not match height grid length (%d)', ...
          numel(prof), numel(dh_m));
end

% Replace NaN/Inf with zero so they contribute nothing to the integral
prof(~isfinite(prof)) = 0;

% Zenith Wet Delay  [m]
ZWD = 1e-6 * sum(prof .* dh_m);

end
