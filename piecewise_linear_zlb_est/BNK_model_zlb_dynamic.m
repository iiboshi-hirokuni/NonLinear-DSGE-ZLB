function [residual, g1, g2, g3] = BNK_model_zlb_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(12, 1);
T34 = (params(8)-1)/(params(3)*params(4))*(params(7)+params(9)/(1-params(10)*exp((-params(5)))));
T47 = (1-params(10)*exp((-params(5))))/params(9);
T82 = (params(7)+params(9)/(1-params(10)*exp((-params(5)))))^(-1);
T86 = exp((-params(5)))*params(10)*params(9)*T82/(1-params(10)*exp((-params(5))));
T104 = (params(9)+exp((-params(5)))*params(10)*params(9)*params(7)*T82/(1-params(10)*exp((-params(5)))))*params(15);
lhs =y(5);
rhs =params(22)*params(6)*exp((1-params(9))*params(5))*y(17)+T34*(y(6)-y(7));
residual(1)= lhs-rhs;
lhs =y(6)-y(7);
rhs =params(21)*(y(18)-y(19)-params(10)*exp((-params(5)))*(y(6)-y(7))-T47*(y(8)/(params(4)*params(2))-y(17)/params(4)-y(10)/params(2)));
residual(2)= lhs-rhs;
lhs =y(9)/(params(4)*params(2));
rhs =params(17)*y(2)/(params(4)*params(2))+(1-params(17))*(y(5)*params(13)/params(4)+(y(6)-y(7))*params(14))+x(it_, 3);
residual(3)= lhs-rhs;
lhs =y(8);
rhs =params(1);
residual(4)= lhs-rhs;
lhs =y(10)/params(2);
rhs =y(7)*params(7)*(1/params(21)-T86)+(1-params(21)*params(16))/params(21)*y(12)+T104*y(11);
residual(5)= lhs-rhs;
lhs =y(7);
rhs =T86*(y(1)-y(11));
residual(6)= lhs-rhs;
lhs =y(11);
rhs =params(15)*y(3)+x(it_, 1);
residual(7)= lhs-rhs;
lhs =y(12);
rhs =params(16)*y(4)+x(it_, 2);
residual(8)= lhs-rhs;
lhs =y(13);
rhs =y(6)+y(11)-y(1);
residual(9)= lhs-rhs;
lhs =y(14);
rhs =y(5)+params(4)-1;
residual(10)= lhs-rhs;
lhs =y(15);
rhs =y(8);
residual(11)= lhs-rhs;
lhs =y(16);
rhs =y(9);
residual(12)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(12, 22);

  %
  % Jacobian matrix
  %

  g1(1,5)=1;
  g1(1,17)=(-(params(22)*params(6)*exp((1-params(9))*params(5))));
  g1(1,6)=(-T34);
  g1(1,7)=T34;
  g1(2,17)=(-(params(21)*(-(T47*(-(1/params(4)))))));
  g1(2,6)=1-params(21)*(-(params(10)*exp((-params(5)))));
  g1(2,18)=(-params(21));
  g1(2,7)=(-1)-params(10)*exp((-params(5)))*params(21);
  g1(2,19)=params(21);
  g1(2,8)=(-(params(21)*(-(T47*1/(params(4)*params(2))))));
  g1(2,10)=(-(params(21)*(-(T47*(-(1/params(2)))))));
  g1(3,5)=(-((1-params(17))*params(13)/params(4)));
  g1(3,6)=(-((1-params(17))*params(14)));
  g1(3,7)=(-((1-params(17))*(-params(14))));
  g1(3,2)=(-(params(17)*1/(params(4)*params(2))));
  g1(3,9)=1/(params(4)*params(2));
  g1(3,22)=(-1);
  g1(4,8)=1;
  g1(5,7)=(-(params(7)*(1/params(21)-T86)));
  g1(5,10)=1/params(2);
  g1(5,11)=(-T104);
  g1(5,12)=(-((1-params(21)*params(16))/params(21)));
  g1(6,1)=(-T86);
  g1(6,7)=1;
  g1(6,11)=T86;
  g1(7,3)=(-params(15));
  g1(7,11)=1;
  g1(7,20)=(-1);
  g1(8,4)=(-params(16));
  g1(8,12)=1;
  g1(8,21)=(-1);
  g1(9,1)=1;
  g1(9,6)=(-1);
  g1(9,11)=(-1);
  g1(9,13)=1;
  g1(10,5)=(-1);
  g1(10,14)=1;
  g1(11,8)=(-1);
  g1(11,15)=1;
  g1(12,9)=(-1);
  g1(12,16)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],12,484);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],12,10648);
end
end
