function [Rtarget_sc_u, Rtarget_sc_uu] = compute_modular_response_second(zc_scaled, zv_scaled, kinetic_law, N, W, ind_ext, alpha_A, alpha_I, alpha_M, v, A, u, c, h, v_plus_fallback, v_minus_fallback, CS_un, CJ_un, target_type)

% [Rtarget_sc_u, Rtarget_sc_uu] = compute_modular_response_second(zc_scaled, zv_scaled, kinetic_law, N, W, ind_ext, alpha_A, alpha_I, alpha_M, v, A, u, c, h, v_plus_fallback, v_minus_fallback, CS_un, CJ_un, target_type)
%
% Efficient calculation of scaled first- and second-order response coefficients (with respect to
% enzyme levels u) of a target function z(c,v), without computing full second-order tensors
%
%   Rtarget_sc_u   (1 x nr)   d ln z / d ln u_l
%   Rtarget_sc_uu  (nr x nr)  d^2 ln z / d ln u_l d ln u_j
%
% target_type (default 'linear'):
%   'linear'     z = zc_scaled' * c + zv_scaled' * v
%                (zc_scaled, zv_scaled are the coefficients, i.e. the unscaled derivatives dz/dc and dz/dv)
%   'loglinear'  ln z = zc_scaled' * ln c + zv_scaled' * ln v
%                (zc_scaled, zv_scaled are the scaled derivatives d ln z/d ln c and d ln z/d ln v)
% For a single concentration or flux (unit vector zc_scaled or zv_scaled), both give the same result.
%
% v_plus_fallback and v_minus_fallback must be given whenever v contains zero values
%
% Formulae (y = [c; v], R^y scaled response coefficients):
%   linear:     R^z_lj = (u_l u_j / z) d^2 z/du_l du_j + delta_lj R^z_l - R^z_l R^z_j
%   loglinear:  R^z_lj = sum_p zhat_p R^p_lj
% where d^2 z/du du is obtained from the scaled second-order elasticities of the reactions, weighted by
% a = d z / d ln v (via control coefficients).
%
% Change 2026-11-09: before, the result was correct only for a single concentration or flux.
%   The help text said "the LOGARITHM of z must be a LINEAR function of the LOGARITHMIC concentrations
%   and fluxes", but the code treated zc_scaled, zv_scaled as coefficients of a linear function, and
%   it used the term - sum_p (y_p zhat_p/z) R^p_l R^p_j instead of - R^z_l R^z_j (correct only if z
%   depends on one variable). Now both target types are computed correctly (tested in
%   mnt/tests/test_modular_elasticities.m).

% upper limit for elasticities and other stuff (to avoid numerical trouble)

eval(default('v_plus_fallback','1', 'v_minus_fallback','1','target_type','''linear'''));

maxv = 10^20;

[nm,nr] = size(N);
if ~length(zc_scaled), zc_scaled = 0 * c; end
if ~length(zv_scaled), zv_scaled = 0 * v; end

switch target_type
  case 'linear'
    % zy_p = y_p dz/dy_p  (so zy/z are the scaled derivatives zhat)
    zy = [c .* zc_scaled; v .* zv_scaled];
    z  = sum(zy);
  case 'loglinear'
    % zy = zhat; with z = 1, the formula below gives sum_p zhat_p R^p_lj
    zy = [zc_scaled; zv_scaled];
    z  = 1;
  otherwise
    error('Unknown target_type');
end

Elasticities = compute_modular_elasticities(...
    kinetic_law,N,W,ind_ext,alpha_A,alpha_I,alpha_M, ...
    v,A,u,c,h,v_plus_fallback,v_minus_fallback,0);

Ec_sc  = lim(Elasticities.sc_E_c, maxv);  % lim : dirty fix
Eu_sc  = lim(Elasticities.sc_E_u, maxv);  % lim : dirty fix
ext    = zeros(nm,1); ext(ind_ext) = 1;

if ~exist('CS_un','var') || isempty(CS_un)
  ext = zeros(size(N,1),1);
  ext(ind_ext) = 1;
  [CJ_un, CS_un] = control_coefficients(N, Elasticities.un_E_c, ext);
end

CS_sc  = diag(1./c) * CS_un * diag(v);
CJ_sc  = diag(1./v) * CJ_un * diag(v);
RSu_sc = CS_sc * Eu_sc;
RJu_sc = CJ_sc * Eu_sc;

a      = zy' * [CS_sc; CJ_sc];

E_target = compute_modular_elasticities_second(...
    a',kinetic_law, N, W, ind_ext, alpha_A, alpha_I, alpha_M, v, ...
    A, u, c, h, v_plus_fallback, v_minus_fallback);

bc     = a * Ec_sc;
bu     = a * Eu_sc;
Bcc    = E_target.sc_E_cc;
Bcu    = E_target.sc_E_cu;
Buc    = E_target.sc_E_cu';
Buu    = E_target.sc_E_uu;

% Limit values (it would be better to understand their reasons ..)

a      = lim(a,  maxv);
bc     = lim(bc, maxv);
bu     = lim(bu, maxv);
Bcc    = lim(Bcc,maxv);
Bcu    = lim(Bcu,maxv);
Buc    = lim(Buc,maxv);
Buu    = lim(Buu,maxv);
RSu_sc = lim(RSu_sc,maxv);
RJu_sc = lim(RJu_sc,maxv);

Rtarget_sc_u = 1/z * column([zy' * [RSu_sc; RJu_sc]])';

% u_l u_j sum_p (zy_p / y_p) d^2 y_p / du_l du_j   (for linear z: u_l u_j d^2 z/du_l du_j)
second_derivative_term = ...
    RSu_sc' * [Bcc + Ec_sc' * diag(a) * Ec_sc - diag(bc)] * RSu_sc ...
    + RSu_sc' * [Bcu + Ec_sc' * diag(a) * Eu_sc] ...
    + [Buc + Eu_sc' * diag(a) * Ec_sc] * RSu_sc ...
    + [Buu + Eu_sc' * diag(a) * Eu_sc - diag(bu)];

switch target_type
  case 'linear'
    Rtarget_sc_uu = 1/z * second_derivative_term + diag(Rtarget_sc_u) - Rtarget_sc_u' * Rtarget_sc_u;
  case 'loglinear'
    Rtarget_sc_uu = second_derivative_term ...
        - [RSu_sc; RJu_sc]' * diag(zy) * [RSu_sc; RJu_sc] ...
        + diag( zy' * [RSu_sc; RJu_sc] );
end


% ---------------------------------------------------
% limit values in X to absolute values below maxv
% replace missing or infinite values by 0

function X = lim(X,maxv);

  if find(abs(X)>maxv)
    warning('Extremely large values encountered; replaing them by threshold values');
  end

  X(~isfinite(X)) = 0;
  X(X>maxv)       = maxv;
  X(X<-maxv)      = -maxv;
