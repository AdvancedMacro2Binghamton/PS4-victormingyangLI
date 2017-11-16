% value function iteration (grid search)
v_tol = 1;
while v_tol >.0001;
   % CONSTRUCT TOTAL RETURN FUNCTION
   v_mat = ret + beta * ...
       repmat(permute(PI * v_guess, [3 2 1]), [num_a 1 1]);
   
   % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
   [vfn, pol_indx] = max(v_mat, [], 2);
   vfn = permute(vfn, [3 1 2]);
   
   v_tol = abs(max(v_guess(:) - vfn(:)));
   
   v_guess = vfn; % update value functions
end;
pol_indx = permute(pol_indx, [3 1 2]);