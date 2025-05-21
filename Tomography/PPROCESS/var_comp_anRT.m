function [P_ll, P_c, P_c_bottom, P_c_top, R_SWD_k, output] = var_comp_anRT(P_ll, P_c, P_c_bottom, P_c_top, A_ll, A_c, Pplus, A_c_top, A_c_bottom, res_SWD, res_Nw, res_Nw_top, res_Nw_bottom, v_top, v_bottom, output, i, iter)
%variance component analysis

Cov_ll_adj = A_ll*Pplus*A_ll';
r_ll = diag(1-Cov_ll_adj*P_ll);

Cov_c_adj = A_c*Pplus*A_c';
r_c = diag(1-Cov_c_adj*P_c);

Cov_c_bottom_adj = A_c_bottom*Pplus*A_c_bottom';
r_c_bottom = diag(1-Cov_c_bottom_adj*P_c_bottom);

Cov_c_top_adj = A_c_top*Pplus*A_c_top';
r_c_top = diag(1-Cov_c_top_adj*P_c_top);

output.sigma_ll(i).sigma_ll(iter) = res_SWD'*P_ll*res_SWD/sum(r_ll);
output.sigma_c(i).sigma_c(iter) = res_Nw'*P_c*res_Nw/sum(r_c);
output.sigma_c_bottom(i).sigma_c_bottom(iter) = res_Nw_bottom'*P_c_bottom*res_Nw_bottom/sum(r_c_bottom);
output.sigma_c_top(i).sigma_c_top(iter) = res_Nw_top'*P_c_top*res_Nw_top/sum(r_c_top);

P_ll = P_ll./output.sigma_ll(i).sigma_ll(iter);
%P_c = P_c./output.sigma_c(i).sigma_c(iter);
P_c_bottom = P_c_bottom./output.sigma_c_bottom(i).sigma_c_bottom(iter);
P_c_top = P_c_top./output.sigma_c_top(i).sigma_c_top(iter);
%P_c = pinv(full(Cov_c));
P_c(v_bottom,v_bottom) = P_c_bottom;
P_c(v_top,v_top) = P_c_top;

R_ll = sqrt(diag(pinv(P_ll)));
R_c = sqrt(diag(pinv(P_c)));
%R_c_bottom = sqrt(diag(pinv(P_c_bottom)));

R_SWD_k = [R_ll;R_c];