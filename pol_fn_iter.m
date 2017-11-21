% policy function iteration by Prof. Kuhn
k = 30;
v_tol = 1;
while v_tol > 1e-6;
   % construct total return function
   v_mat = ret + beta * repmat(permute(PI * v_guess, [3 2 1]), [num_a 1 1]);
   % choose highest value (associated with a' choice)
   [vfn, pol_indx] = max(v_mat, [], 2);
   vfn = permute(vfn, [3 1 2]);
   pol_indx = permute(pol_indx, [3 1 2]);
   v_tol = abs(max(v_guess(:) - vfn(:)));
   v_guess = vfn; % update value functions
   % construct Q matrix from policy index and PI
   Q_mat = makeQmatrix(pol_indx, PI);
   % construc return vector and value vector
   pol_fn = a(pol_indx);
   u_mat = bsxfun(@minus, r * a, pol_fn);
   u_mat = bsxfun(@plus, u_mat, z' * w);
   u_mat = (u_mat .^ (1 - sigma)) ./ (1 - sigma);
   u_vec = u_mat(:);
   w_vec = v_guess(:);
   % PFI
   for ii = 1:k
       w_vec_new = u_vec + beta * Q_mat * w_vec;
       w_vec = w_vec_new;
   end
   v_guess = reshape(w_vec, num_z, num_a);
end;