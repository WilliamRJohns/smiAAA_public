function [e_rms, e_rel, H2] = comp_error(f_act, f_approx)
  diff = abs(f_act - f_approx);
  % Get noramlizing constant N = Ne * Ns
  N = prod(size(f_act), 'all');
  e_rms = sqrt(sum(diff.^2,'all') / N);
  e_rel = sum(diff ./ abs(f_act), 'all') * 100 / N;
  H2=sqrt(norm(diff,'fro')^2/norm(f_act,'fro')^2);
 