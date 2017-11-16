% value function iteration interpolation
v_tol = 1;
while v_tol >.0001;
   % CONSTRUCT TOTAL RETURN FUNCTION
   v_mat = ret + beta * ...
       repmat(permute(PI * v_guess, [3 2 1]), [num_a 1 1]);
   
   % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
   [vfn, pol_indx] = max(v_mat, [], 2);
   % only need to consider the interval (pol_ind-1, pol_ind+1) since
   % the max must be in this range
   vfn = permute(vfn, [3 1 2]);
   
   v_tol = abs(max(v_guess(:) - vfn(:)));
   
   v_guess = vfn; %update value functions
end;
pol_indx = permute(pol_indx, [3 1 2]);
% aa = a_lo:.25:a_hi;
% vv = spline(a,vfn,aa);
% plot(a,vfn,'o',aa,vv)

% x = [0 1 2.5 3.6 5 7 8.1 10];
% y = sin(x);
% xx = 0:.25:10;
% yy = spline(x,y,xx);
% plot(x,y,'o',xx,yy)