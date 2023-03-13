function [residual, g1, g2] = BNK_model_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 12, 1);

%
% Model equations
%

T33 = (params(8)-1)/(params(3)*params(4))*(params(7)+params(9)/(1-params(10)*exp((-params(5)))));
T43 = (1-params(10)*exp((-params(5))))/params(9);
T75 = (params(7)+params(9)/(1-params(10)*exp((-params(5)))))^(-1);
T79 = exp((-params(5)))*params(10)*params(9)*T75/(1-params(10)*exp((-params(5))));
T97 = (params(9)+exp((-params(5)))*params(10)*params(9)*params(7)*T75/(1-params(10)*exp((-params(5)))))*params(15);
lhs =y(1);
rhs =y(1)*params(22)*params(6)*exp((1-params(9))*params(5))+T33*(y(2)-y(3));
residual(1)= lhs-rhs;
lhs =y(2)-y(3);
rhs =params(21)*(y(2)-y(3)-params(10)*exp((-params(5)))*(y(2)-y(3))-T43*(y(4)/(params(4)*params(2))-y(1)/params(4)-y(6)/params(2)));
residual(2)= lhs-rhs;
lhs =y(5)/(params(4)*params(2));
rhs =y(4)/(params(4)*params(2))*params(17)+(1-params(17))*(y(1)*params(13)/params(4)+(y(2)-y(3))*params(14))+x(3);
residual(3)= lhs-rhs;
lhs =y(4);
rhs =y(5);
residual(4)= lhs-rhs;
lhs =y(6)/params(2);
rhs =y(3)*params(7)*(1/params(21)-T79)+(1-params(21)*params(16))/params(21)*y(8)+T97*y(7);
residual(5)= lhs-rhs;
lhs =y(3);
rhs =T79*(y(2)-y(7));
residual(6)= lhs-rhs;
lhs =y(7);
rhs =params(15)*y(7)+x(1);
residual(7)= lhs-rhs;
lhs =y(8);
rhs =params(16)*y(8)+x(2);
residual(8)= lhs-rhs;
lhs =y(9);
rhs =y(2)+y(7)-y(2);
residual(9)= lhs-rhs;
lhs =y(10);
rhs =y(1)+params(4)-1;
residual(10)= lhs-rhs;
lhs =y(11);
rhs =y(4);
residual(11)= lhs-rhs;
lhs =y(12);
rhs =y(5);
residual(12)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(12, 12);

  %
  % Jacobian matrix
  %

  g1(1,1)=1-params(22)*params(6)*exp((1-params(9))*params(5));
  g1(1,2)=(-T33);
  g1(1,3)=T33;
  g1(2,1)=(-(params(21)*(-(T43*(-(1/params(4)))))));
  g1(2,2)=1-(1-params(10)*exp((-params(5))))*params(21);
  g1(2,3)=(-1)-params(21)*((-1)-(-(params(10)*exp((-params(5))))));
  g1(2,4)=(-(params(21)*(-(T43*1/(params(4)*params(2))))));
  g1(2,6)=(-(params(21)*(-(T43*(-(1/params(2)))))));
  g1(3,1)=(-((1-params(17))*params(13)/params(4)));
  g1(3,2)=(-((1-params(17))*params(14)));
  g1(3,3)=(-((1-params(17))*(-params(14))));
  g1(3,4)=(-(params(17)*1/(params(4)*params(2))));
  g1(3,5)=1/(params(4)*params(2));
  g1(4,4)=1;
  g1(4,5)=(-1);
  g1(5,3)=(-(params(7)*(1/params(21)-T79)));
  g1(5,6)=1/params(2);
  g1(5,7)=(-T97);
  g1(5,8)=(-((1-params(21)*params(16))/params(21)));
  g1(6,2)=(-T79);
  g1(6,3)=1;
  g1(6,7)=T79;
  g1(7,7)=1-params(15);
  g1(8,8)=1-params(16);
  g1(9,7)=(-1);
  g1(9,9)=1;
  g1(10,1)=(-1);
  g1(10,10)=1;
  g1(11,4)=(-1);
  g1(11,11)=1;
  g1(12,5)=(-1);
  g1(12,12)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],12,144);
end
end
