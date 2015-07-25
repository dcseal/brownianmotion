% This script produces plots for the paper, 
% "Merton's Model with Stochastic Recovery", A. Cohen, N. Costanzino,
% (to be submitted).
%
% TODO - save the seed value for the random number generator.  This is
% important to reproduce the plots.
%
% Author: David Seal, Michigan State University, Fall 2014

% Set the constants for this problem:
f     = 1.0;
sigA  = 0.25;  
sigR  = 0.25;
rho   = linspace( 0.0, 1.0 );
r     = 0.05; 
A0    = 1.25; 
R0    = 1.2;
N     = 1.0;
T     = 10.0;  % Final time

% Derived constants
gamma = rho.*(sigR / sigA );

t   = 0;                    % Elapsed time
tau = T - t;                % Time to maturity

% Asset price (constant for this problem)
At = 0*t + A0;

% Recovery price (constant for this problem)
Rt = 0*t + R0;

% Eqn. (2.4, and 2.5)
d0 = ( log( At/N ) + (r-0.5*sigA^2)*tau ) ./ ( sigA*sqrt(tau) );
d1 = ( log( At/N ) + (r+0.5*sigA^2)*tau ) ./ ( sigA*sqrt(tau) );

% Derived value (Eqn. ... )
dgam  = d0 + gamma.*sigA*sqrt(tau);

% Eqns. (2.3, , 2.8, 2.10)
B_MERt   = N*exp( -r*tau ) .* Normal( d0 ) + f * At .* Normal( -d1 );
S_MERt   = -log( Normal(d0) + (f*At)./(N*exp(-r*tau)).*Normal(-d1 ) )./tau;
LGD_MERt = 1.0-exp(r*tau).*( f*At/N ).* ( Normal( -d1 )./Normal( -d0 ) );

% Eqn. (3.15)
Bt = N*exp(-r*tau) .* Normal( d0 ) + Rt .* Normal( -dgam );

% Eqn. (3.20)
LGDt     = 1.0-exp(r*tau).*( Rt/N   ).* ( Normal(-dgam) ./ Normal( -d0 ) );

% Credit spread (Eqn. 3.21)
St   =  (1./tau).*log( N*exp(-r*tau ) ./ Bt );

d0B   = B_MERt - Bt;
d0S   = S_MERt - St;
d0LGD = LGD_MERt - LGDt;

% Alternatively, we can use Eqn. 4.1 and 4.2
Theta = Normal( -d1+(1.0-gamma).*sigA*sqrt(tau) ) - Normal( -d1 );
dB    = (f*At-Rt).*Normal( -d1 ) - Theta .* Rt;
dS    = -log( (N*exp(-r*tau) .* Normal( d0 ) + f*At.*Normal(-d1))./(N*exp(-r*tau).*Normal(d0) + Rt.*Normal(-dgam) ) )./tau;
dLGD  = ( (Rt-f*At).*Normal(-d1) + Rt.*Theta ) ./ ( N*exp(-r*tau).*Normal(-d0) );

% Eqn. (5.13)
%Rswap = ( exp(r*(tau)) .* Rt / N ) .* ( Normal(-dgam) ./ Normal( -d0 ) );

% ------------------------------- %
% Compute the absolute difference
% ------------------------------- %

figure(1);
clf;
p1 = plot(  rho,   dB, '--b' );
set( p1, 'LineWidth', 3 );
hold on;

p2 = plot(  rho,   dS, '-g' );
set( p2, 'LineWidth', 3 );

p3 = plot(  rho, dLGD, '-xk' );
set( p3, 'LineWidth', 1 );

l1 = legend('\Delta B_{t,T}', '\Delta S_{t,T}', '\Delta LGD_{t,T}' );
set(l1, 'FontSize', 16 );
set(l1, 'Location', 'Best' );
hold off;

xlabel('Correlation ( \rho )', 'Fontsize', 14 );
ylabel('Price per Unit', 'Fontsize', 14 );

t1 = title('Absolute Difference');
set( t1, 'Fontsize', 16 );

% ------------------------------- %
% Compute the relative difference
% ------------------------------- %

%   figure(2);
%   clf;
%   p1 = plot(  tau,   dB./B_MERt, '--b' );
%   set( p1, 'LineWidth', 3 );
%   hold on;

%   p2 = plot(  tau,   dS./S_MERt, '-g' );
%   set( p2, 'LineWidth', 3 );

%   p3 = plot(  tau, dLGD./LGD_MERt, '-xk' );
%   set( p3, 'LineWidth', 1 );

%   l1 = legend('\Delta B_{t,T} / B_{t,T}^{MER}', '\Delta S_{t,T} / S_{t,T}^{MER}', '\Delta LGD_{t,T} / LGD^{MER}_{t,T}' );
%   set(l1, 'FontSize', 16 );
%   set(l1, 'Location', 'Best' );
%   hold off;

%   xlabel('Time to maturity (\tau )',  'Fontsize', 14 );
%   ylabel('Price per Unit', 'Fontsize', 14 );

%   t1 = title('Relative Difference');
%   set( t1, 'Fontsize', 16 );
