function G = grids_even(O,P)

% G = grids_even(O,P)
%   Constructs an evenly-spaced discretized state space
% Inputs:
%	O :   Structure of user-specified options
%   P :   Structure of parameters
% Output:
%   G :   Structure of grid points

%--------------------------------------------------------------------------
% Set up grids by discretizing the state space
%--------------------------------------------------------------------------
% Define grid points for interest rate
G.r_grid = linspace(O.rbound(1),O.rbound(2),O.r_pts);
G.MP_grid = linspace(O.MPbound(1),O.MPbound(2),O.MP_pts);

% Define grid points for preference shock
G.b_grid = linspace(O.bbound(1),O.bbound(2),O.b_pts);

% Define grid points for productivity shock
G.a_grid = linspace(O.abound(1),O.abound(2),O.a_pts);

%--------------------------------------------------------------------------
% Weights for numerical integration from truncated normal
%--------------------------------------------------------------------------
% [e_nodes,G.e_weight] = ghquad(O.e_pts);
% G.e_nodes = (2^.5) * P.sigma_g * e_nodes;
% [u_nodes,G.u_weight] = ghquad(O.u_pts);
% G.u_nodes = (2^.5) * P.sigma_z * u_nodes;
% 
% [v_nodes,G.v_weight] = ghquad(O.v_pts);
% G.v_nodes = (2^.5) * P.sigma_r * v_nodes;

G.e_nodes =0;
G.u_nodes = 0;
G.v_nodes = 0;
G.w_nodes = 0;
%--------------------------------------------------------------------------
% Construct grid arrays
%--------------------------------------------------------------------------

[G.r_gr G.b_gr G.a_gr G.MP_gr] = ndgrid(G.r_grid, G.b_grid, G.a_grid, G.MP_grid);

G.nodes = numel(G.r_gr);
G.griddim = size(G.r_gr);

