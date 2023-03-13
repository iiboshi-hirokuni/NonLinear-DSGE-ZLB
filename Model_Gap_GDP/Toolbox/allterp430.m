function [o1,o2,o3] = allterp430(  x1, x2, x3, x4,...
                                   x1i, x2i, x3i, x4i,...
                                   pf1 ,pf2, pf3)

% allterp332 linear inter/extrapolation (3 states, 3 policies, 2 stoch comp)
% Inputs:
%   x*      :   Grid
%   x*i     :   Point to evaluate
%   pf*     :   Policy function
% Outputs:
%   o*      :   Interpolated/extrapolated values of dimension x4ipts

% Grid lengths
nx1 = length(x1);
nx2 = length(x2);
nx3 = length(x3);
nx4 = length(x4);

% Number of stochastic realizations
% x1ipts = length(x1i);
% x2ipts = length(x2i);
% x3ipts = length(x3i);
x1ipts = 1;
x2ipts = 1;
x3ipts = 1;
x4ipts = 1;

% Preallocate output
o1 = zeros(x1ipts,x2ipts,x3ipts,x4ipts);
o2 = zeros(x1ipts,x2ipts,x3ipts,x4ipts);
o3 = zeros(x1ipts,x2ipts,x3ipts,x4ipts);

% for i3 = 1:x3ipts
%   for i2 = 1:x2ipts 
%      for i1 = 1:x1ipts   
i1 =1;
i2 =1;
i3 = 1;
i4 =1;

        s1 = x1(2) - x1(1);
        x1i_min = x1i(1) - x1(1);                     
        loc1 = real(min(nx1-1,max(1,floor(x1i_min/s1) + 1)));
        
        s2 = x2(2) - x2(1);
        x2i_min = x2i(1) - x2(1);
        loc2 = real(min(nx2-1,max(1,floor(x2i_min/s2) + 1)));

        s3 = x3(2) - x3(1);
        x3i_min = x3i(1) - x3(1);
        loc3 = real(min(nx3-1,max(1,floor(x3i_min/s3) + 1)));
        
        s4 = x4(2) - x4(1);
        x4i_min = x4i(1) - x4(1);
        loc4 = real(min(nx4-1,max(1,floor(x4i_min/s4) + 1)));
        
        xi = [x1i(1) x2i(1) x3i(1) x4i(1)];
        xi_left = [x1(loc1) x2(loc2) x3(loc3) x4(loc4)];
        xi_right = [x1(loc1+1) x2(loc2+1) x3(loc3+1) x4(loc4+1) ];

        w_2 = (xi - xi_left)./ (xi_right - xi_left);
        w_1 = 1 - w_2;
        w1 = [w_1(1) w_2(1)];
        w2 = [w_1(2) w_2(2)];
        w3 = [w_1(3) w_2(3)];
        w4 = [w_1(4) w_2(4)];

        for m4 = 0:1
          for m3 = 0:1
              for m2 = 0:1
                  for m1 = 0:1
                      
                     o1(i1,i2,i3,i4) = o1(i1,i2,i3,i4) + ...
                            w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*...
                            pf1(loc1+m1,loc2+m2,loc3+m3,loc4+m4);
                 
                      o2(i1,i2,i3,i4) = o2(i1,i2,i3,i4) + ...
                            w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*...
                            pf2(loc1+m1,loc2+m2,loc3+m3,loc4+m4);

                     o3(i1,i2,i3,i4) = o3(i1,i2,i3,i4) + ...
                            w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*...
                            pf3(loc1+m1,loc2+m2,loc3+m3,loc4+m4);
                  end
              end
          end    
        end
        
        
%      end
%   end
% end