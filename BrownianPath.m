function [W,t] = BrownianPath( T, N, seed_num, plt_flag )
%BROWNIANPATH.   Brownian path simulation.
%
% Usage: [W] = BrownianPath(T,N,seed_num,plt_flag).
%
% This script was pulled from slides available at
% http://www.math.ksu.edu/~dmaldona/research/slides/SDE_FTSS.pdf
%
% Input
% -----
%
%       T           : Final time
%       N           : Number of time steps
%       seed        : Seed for random number generator
%       plt_flg     : Plotting flag.  If any argument is passed in, plot results.
%
% Returns
% -------
%
%       W  : (Discrete) Brownian motion vector.
%            size( W ) = (1,N).
%            

% set the state of random number generator
% Save the seed value in case we wish to reproduce these results!
if( nargin > 2 & seed_num > 0 )
    rng(seed_num);
end

dt = T/N;                   % Discrete time steps
dW = sqrt(dt)*randn(1,N);   % increments
W  = cumsum(dW);            % cumulative sum
t  = linspace( 0, T, N );

% Plot the solution
if( nargin > 3 )

    M = max(abs(W));            % Used for setting axes
    plot([0:dt:T],[0,W],'b')
    grid on
    xlabel('t','FontSize',16)
    ylabel('W(t)','FontSize',16,'Rotation',0)
    title(['Discretized Brownian path from time 0 to ' int2str(T) ' with ' int2str(N) ' steps']) 
    axis([0 T -M-0.5 M+0.5]);

end

end
