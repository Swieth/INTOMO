function ZWD = excessphasezenith(prof, height)
% EXCESSPHASEZENITH  Compute Zenith Wet Delay from a wet refractivity profile
%
% Integrates a vertical wet refractivity profile (in N-units) over a
% height grid to produce the Zenith Wet Delay (ZWD) in metres.
%
% Standard formula:
%   ZWD = 1e-6 * sum( Nw_i * dh_i )
% where dh_i is the layer thickness in metres.
%
% INPUT
%   prof   - [nH x 1] wet refractivity profile (N-units, dimensionless x1e6)
%   height - [nH x 1] layer thicknesses in km (from diff(model.levels_TOMO)./1000)
%
% OUTPUT
%   ZWD    - Zenith Wet Delay in metres
%
% See also: stackedNwcalc, weightObs

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
