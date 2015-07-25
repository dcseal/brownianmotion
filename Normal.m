function f = Normal( x )
% Compute the cumulative distribution function for a Normal random variable.
%
% This routine computes
%
%   NORMAL(x) = 1/sqrt(2*pi) * \int_{-inf}^x exp(-t^2/2) dt.
%
% Note that Normal( -inf ) = 0, Normal(0) = 1/2, and Normal( inf ) = 1.
%
    f = 0.5*(1.+ erf( x/sqrt(2.) ) );
end
