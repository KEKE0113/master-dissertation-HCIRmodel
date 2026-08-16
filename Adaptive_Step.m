function h = Adaptive_Step(state, speed, level, tol, hmax)
% ADAPTIVE_STEP Determines the adaptive time step based on the local drift.

h = min(hmax, tol ./ (speed^2 * abs(level - state)));

end
