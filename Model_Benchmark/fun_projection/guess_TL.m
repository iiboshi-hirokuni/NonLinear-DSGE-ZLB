function pf = guess_TL(P,S,G)

% pf = guess(P,S,G) 
%   Sets the initial policy functions
% Inputs:
%   P : structure of parameters
%   S : structure of steady state values 
%   G : structure of grids
% Outputs:
%  pf : structure of policy functions

%----------------------------------------------------------------------
% Log-linear solution - ZLB not imposed
%----------------------------------------------------------------------
V = variables;
[T,T0,eu] = linmodel(P,S,V);

%%  check whether solutions exist or not.
% disp(['eu:    ' mat2str(eu)]) % "0" is OK.

% Transform discretized state space to percent deviation from steady state     
r_gr_per = (G.r_gr - S.R)./1;   
b_gr_per = (G.b_gr - 0)./1;
a_gr_per = (G.a_gr - 0)./1;
MP_gr_per = (G.MP_gr - 0) ./1;

% Calculate linear policy functions on discretized state space    
linpf_y = zeros(G.griddim);
linpf_pi = zeros(G.griddim);
linpf_r = zeros(G.griddim);
linpf_y_star = zeros(G.griddim);

state = [r_gr_per(:) b_gr_per(:) a_gr_per(:) MP_gr_per(:)  ]';    
TT = [T T0];

% output(Y_t) - potential(Y*_t) 
linpf_y(:) = TT(V.y,[V.r V.z_b V.mu_a V.nvar+V.eps_MP ])*state ...  % + T0(V.y,V.eps_MP)* [MP_gr_per(:)]';  %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA
            -TT(V.y_star,[V.r V.z_b V.mu_a V.nvar+V.eps_MP ])*state;                                    %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA

linpf_pi(:) = TT(V.pi,[V.r V.z_b V.mu_a V.nvar+V.eps_MP ])*state;  %+T0(V.pi,V.eps_MP)*MP_gr_per(:);    %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA
linpf_r(:) = TT(V.r,[V.r V.z_b V.mu_a V.nvar+V.eps_MP ])*state;    %+T0(V.pi,V.eps_MP)*MP_gr_per(:);    %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA


% linpf_y(:) = TT(V.y,[V.r V.nvar+V.eps_b V.nvar+V.eps_a V.nvar+V.eps_MP ])*state ...  % + T0(V.y,V.eps_MP)* [MP_gr_per(:)]';
%             -TT(V.y_star,[V.r V.nvar+V.eps_b V.nvar+V.eps_a V.nvar+V.eps_MP ])*state;  
% 
% linpf_pi(:) = TT(V.pi,[V.r V.nvar+V.eps_b V.nvar+V.eps_a V.nvar+V.eps_MP ])*state;  %+T0(V.pi,V.eps_MP)*MP_gr_per(:);
% linpf_r(:) = TT(V.r,[V.r V.nvar+V.eps_b V.nvar+V.eps_a V.nvar+V.eps_MP ])*state;    %+T0(V.pi,V.eps_MP)*MP_gr_per(:);

% T0
% 
% T


% Convert back to levels
pf.y =  S.y_star*(1 + linpf_y);
pf.y_star =  S.y_star*(1 + linpf_y_star);
pf.pi = S.pi_star + linpf_pi;
pf.r  = S.R  + linpf_r;

% pf.pi = S.pi_star*(1 + linpf_pi);
% pf.r  = S.R *(1 + linpf_r);




