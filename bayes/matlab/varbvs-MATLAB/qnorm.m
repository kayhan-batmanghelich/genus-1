% qnorm(x,A) returns the quadratic norm of vector x with respect to positive
% definite matrix A. For a definition of the quadratic norm, see p. 635 of
% Convex Optimization (2004) by Boyd & Vandenberghe.
function y = qnorm (x, A)

  % Part of the varbvs package, https://github.com/pcarbo/varbvs
  %
  % Copyright (C) 2012-2017, Peter Carbonetto
  %
  % This program is free software: you can redistribute it under the
  % terms of the GNU General Public License; either version 3 of the
  % License, or (at your option) any later version.
  %
  % This program is distributed in the hope that it will be useful, but
  % WITHOUT ANY WARRANY; without even the implied warranty of
  % MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  % General Public License for more details.
  %
  x = x(:);
  y = sqrt(x'*A*x);
  