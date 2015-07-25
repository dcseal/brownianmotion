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
N     = 1.0;
sigA  = 0.25;  
sigR  = 0.25;
rho   = 0.25; 
r     = 0.05; 
A0    = 1.25; 
R0    = 1.2;
T     = 10.0;  % Final time

% Derived constants
gamma = rho *(sigR / sigA );

NumSteps = 1000;  % Number of steps used for Brownian path

NumSimulations = 20;       % Number of simulations to run

for n1=1:NumSimulations

    % Simulate two separate Brownian motion paths
    [ZtA,t ] = BrownianPath( T, NumSteps );
    [Wt, t2] = BrownianPath( T, NumSteps );
    ZtR = rho * ZtA + sqrt( 1 - rho^2 )*Wt;

    % Waiting time
    tau = T-t;

    % Asset price (Eqn. 3.5)
    At = A0*exp( (r-0.5*sigA^2)*t + sigA*ZtA );
    At = 0*At + A0;

    % Recovery price (Eqn. 3.5)
    Rt = R0*exp( (r-0.5*sigR^2)*t2 + sigR*ZtR );
    Rt = 0*Rt + R0;

    % Eqn. (2.5)
    d0    = ( log( At / N ) + (r-0.5*sigA^2)*tau ) ./ ( sigA*sqrt(tau) );
    dgam  = d0 + gamma*sigA*sqrt(T-t);

    % Eqn. (3.15)
    Bt = N*exp(-r*(T-t)) .* Normal( d0 ) + Rt .* Normal( -dgam );

    % Eqn. (5.13)
    Rswap = ( exp(r*(T-t)) .* Rt / N ) .* ( Normal(-dgam) ./ Normal( -d0 ) );

    % Credit spread (Eqn. 3.21)
    %St = -(1./(T-t)).*( Normal(d0) + ( exp(r*(T-t)).*Rt/N ) .* Normal(-dgam) );
    St = (1./(T-t)).*log( N*exp(-r*(T-t) ) ./ Bt );

    figure(1);
    clf;
    p1 = plot( t, At, '-r' );

    hold on;

    p2 = plot( t2, Rt, '-b' );
    p3 = plot(  t, Bt, '-g' );
    p4 = plot(  t, St, '-k' );

    l1 = legend('A_t', 'R_t', 'B_{t,T}', 'S_{t,T}' );
    set(l1, 'FontSize', 16 );
    % set(l1, 'Location', 'NorthWest' );
    set(l1, 'Location', 'Best' );
    hold off;

    xlabel('Time',  'Fontsize', 14 );
    ylabel('Price per Unit', 'Fontsize', 14 );

    t1 = title('Asset vs. Recovery vs. SRM Bond Price');
    set( t1, 'Fontsize', 16 );

    fname = strcat('AssetVsRecovery', num2str(n1, '%02d' ) );
    fname = strcat(fname, '.eps' );

    % Save the pretty picture!
    print(1, '-depsc', fname );

end
