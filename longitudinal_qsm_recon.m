function [x0_ppm, dx_ppm, cost, cost_hist] = longitudinal_qsm_recon(varargin)
%LONGITUDINAL_QSM_RECON Joint longitudinal QSM reconstruction for two time points.
%
%   [x0_ppm, dx_ppm, cost, cost_hist] = longitudinal_qsm_recon(...)
%
% This function jointly reconstructs:
%   x0 : baseline susceptibility map
%   dx : temporal susceptibility change map
%
% Forward model:
%   RDF1 = D1 * x0
%   RDF2 = D2 * (x0 + dx)
%
% Cost function:
%   data fidelity
% + MEDI regularization on x0 using iMag1
% + MEDI regularization on (x0 + dx) using iMag2
% + lambda_dx * |dx| within the valid brain region
%
% Required inputs:
%   RDF1, RDF2      : local field maps in radians
%   Mask1, Mask2    : brain masks
%   voxel_size      : voxel size
%   delta_TE        : echo spacing in seconds
%   CF              : center frequency in Hz
%   B0_dir1, B0_dir2: B0 direction unit vectors
%
% Optional inputs:
%   iMag1, iMag2    : magnitude images for structural regularization
%   iMag            : fallback magnitude image for backward compatibility
%   VesselMask      : mask where temporal sparsity is disabled
%   Mask_CSF        : CSF mask for referencing final outputs
%   x0_init_ppm     : initial x0 in ppm (use conventional medi)
%   dx_init_ppm     : initial dx in ppm (use conventional medi)
%
% Notes:
% - Temporal sparsity is applied only within the intersection of Mask1 and Mask2.
% - Sparsity can be disabled inside VesselMask.
% - Separate magnitude-guided structural weights are used for x0 and x0 + dx.
% - Optional CSF referencing can be applied to the final output.
%
% External dependencies:
%   dipole_kernel, sphere_kernel, SMV, dataterm_mask,
%   gradient_mask, reg_MEDI, fgrad, bdiv

%% Parse inputs
RDF1  = get_option(varargin, 'RDF1');
RDF2  = get_option(varargin, 'RDF2');

Mask1 = logical(get_option(varargin, 'Mask1'));
Mask2 = logical(get_option(varargin, 'Mask2'));

vox = get_option(varargin, 'voxel_size');
dTE = get_option(varargin, 'delta_TE');
CF  = get_option(varargin, 'CF');

H1 = get_option(varargin, 'B0_dir1');
H2 = get_option(varargin, 'B0_dir2');

iMag1 = get_option(varargin, 'iMag1', []);
iMag2 = get_option(varargin, 'iMag2', []);
iMag  = get_option(varargin, 'iMag', []);  % backward compatibility
if isempty(iMag1), iMag1 = iMag; end
if isempty(iMag2), iMag2 = iMag; end

lambda    = get_option(varargin, 'lambda', 2000);
lambda_dx = get_option(varargin, 'lambda_dx', 3);

outer_irls = get_option(varargin, 'outer_irls', 5);
cg_max_it  = get_option(varargin, 'cg_max_iter', 120);
cg_tol     = get_option(varargin, 'cg_tol', 1e-4);
linear     = get_option(varargin, 'linear', false);

use_smv = get_option(varargin, 'smv', false);
radius  = get_option(varargin, 'radius', 5);

VesselMask = get_option(varargin, 'VesselMask', []);
Mask_CSF   = get_option(varargin, 'Mask_CSF', []);

x0_init_ppm = get_option(varargin, 'x0_init_ppm', []);
dx_init_ppm = get_option(varargin, 'dx_init_ppm', []);

lambda_dx_in = get_option(varargin, 'lambda_dx_in', 1);

%% Define valid region for temporal sparsity
MaskU = logical(Mask1 & Mask2);

lambda_dx_eff = lambda_dx * lambda_dx_in * single(MaskU);
if ~isempty(VesselMask)
    lambda_dx_eff = lambda_dx_eff .* single(~VesselMask);
end
valid_mask = logical(lambda_dx_eff);

%% Dipole kernels and optional SMV filtering
D1 = dipole_kernel(size(RDF1), vox, H1);
D2 = dipole_kernel(size(RDF2), vox, H2);

if use_smv
    SphereK = single(sphere_kernel(size(RDF1), vox, radius));
    D1 = (1 - SphereK) .* D1;
    D2 = (1 - SphereK) .* D2;

    RDF1 = RDF1 - SMV(RDF1, SphereK);
    RDF2 = RDF2 - SMV(RDF2, SphereK);

    RDF1 = RDF1 .* Mask1;
    RDF2 = RDF2 .* Mask2;
end

Dconv1 = @(x) real(ifftn(D1 .* fftn(x)));
Dconv2 = @(x) real(ifftn(D2 .* fftn(x)));

%% Data weighting
m1 = dataterm_mask(1, ones(size(RDF1)), Mask1);
m2 = dataterm_mask(1, ones(size(RDF2)), Mask2);
W1 = abs(m1).^2;
W2 = abs(m2).^2;

%% MEDI structural weights
opts1.grad = @fgrad;
opts1.div = @bdiv;
opts1.percentage = 0.9;
opts1.voxel_size = vox;
opts1.B0_dir = H1;
opts1.matrix_size = size(RDF1);
opts1.wG = gradient_mask(1, iMag1, Mask1, opts1.grad, vox, opts1.percentage);

opts2.grad = @fgrad;
opts2.div = @bdiv;
opts2.percentage = 0.9;
opts2.voxel_size = vox;
opts2.B0_dir = H2;
opts2.matrix_size = size(RDF2);
opts2.wG = gradient_mask(1, iMag2, Mask2, opts2.grad, vox, opts2.percentage);

%% Unit conversion
scale = (2 * pi * dTE * CF) / 1e6;

if isempty(x0_init_ppm)
    x0 = zeros(size(RDF1), 'single');
else
    x0 = single(x0_init_ppm) * scale;
end

if isempty(dx_init_ppm)
    dx = zeros(size(RDF1), 'single');
else
    dx = single(dx_init_ppm) * scale;
end

%% History
cost_hist = struct( ...
    'data1', nan(1, outer_irls), ...
    'data2', nan(1, outer_irls), ...
    'TVx0', nan(1, outer_irls), ...
    'TVjoint', nan(1, outer_irls), ...
    'dxL1', nan(1, outer_irls), ...
    'total', nan(1, outer_irls));

%% Outer IRLS loop
for out = 1:outer_irls
    t_start = tic;

    reg0 = reg_MEDI(opts1, x0, @(z) z);
    xsum = x0 + dx;
    regJ = reg_MEDI(opts2, xsum, @(z) z);

    % IRL1 surrogate weights for dx
    w_dx_val = 1 ./ (abs(dx) + 1e-6);
    w_dx_val(~valid_mask) = 0;

    % Data term
    A00_data = @(z) Dconv1(W1 .* Dconv1(z)) + Dconv2(W2 .* Dconv2(z));
    Add_data = @(z) Dconv2(W2 .* Dconv2(z));
    A0d_data = Add_data;
    Ad0_data = Add_data;

    if linear
        r1 = Dconv1(x0) - RDF1;
        r2 = Dconv2(x0 + dx) - RDF2;

        b0_data = 2 * lambda * (Dconv1(W1 .* r1) + Dconv2(W2 .* r2));
        bd_data = 2 * lambda * Dconv2(W2 .* r2);
    else
        phi1 = Dconv1(x0);
        phi2 = Dconv2(x0 + dx);

        w1  = m1 .* exp(1i * phi1);
        b01 = m1 .* exp(1i * RDF1);

        w2  = m2 .* exp(1i * phi2);
        b02 = m2 .* exp(1i * RDF2);

        r1_term = real(conj(w1) .* (-1i) .* (w1 - b01));
        r2_term = real(conj(w2) .* (-1i) .* (w2 - b02));

        b0_data = 2 * lambda * (Dconv1(r1_term) + Dconv2(r2_term));
        bd_data = 2 * lambda * Dconv2(r2_term);
    end

    % Normal equations
    A00 = @(z) 2 * lambda * A00_data(z) + reg0.opHop(z) + regJ.opHop(z);
    Add = @(z) 2 * lambda * Add_data(z) + ...
               valid_mask .* (lambda_dx_eff .* (w_dx_val .* z)) + regJ.opHop(z);
    A0d = @(z) 2 * lambda * A0d_data(z) + regJ.opHop(z);
    Ad0 = @(z) 2 * lambda * Ad0_data(z) + regJ.opHop(z);

    b0 = b0_data + reg0.opHop(x0) + regJ.opHop(xsum);
    bd = bd_data + valid_mask .* (lambda_dx_eff .* (w_dx_val .* dx)) + regJ.opHop(xsum);

    [d0, dd] = pcg_block( ...
        @(z0, zd) deal(A00(z0) + Ad0(zd), A0d(z0) + Add(zd)), ...
        -b0, -bd, ...
        zeros(size(RDF1), 'single'), zeros(size(RDF1), 'single'), ...
        cg_tol, cg_max_it);

    x0 = x0 + d0;
    dx = dx + dd;

    % True cost
    if linear
        t1 = m1 .* (Dconv1(x0) - RDF1);
        t2 = m2 .* (Dconv2(x0 + dx) - RDF2);
    else
        w1  = m1 .* exp(1i * Dconv1(x0));
        b01 = m1 .* exp(1i * RDF1);

        w2  = m2 .* exp(1i * Dconv2(x0 + dx));
        b02 = m2 .* exp(1i * RDF2);

        t1 = w1 - b01;
        t2 = w2 - b02;
    end

    dx_vec  = dx(valid_mask);
    lam_vec = lambda_dx_eff(valid_mask);

    cost_hist.data1(out)   = norm(t1(:), 2);
    cost_hist.data2(out)   = norm(t2(:), 2);
    cost_hist.TVx0(out)    = reg0.cost(x0);
    cost_hist.TVjoint(out) = regJ.cost(xsum);
    cost_hist.dxL1(out)    = sum(lam_vec(:) .* abs(dx_vec(:)));
    cost_hist.total(out)   = cost_hist.data1(out)^2 + cost_hist.data2(out)^2 + ...
                             cost_hist.TVx0(out) + cost_hist.TVjoint(out) + ...
                             cost_hist.dxL1(out);

    fprintf(['OUT %d/%d | Time=%.1fs | TVx0=%.3e | TV(x0+dx)=%.3e | ', ...
             'dxL1=%.3e | TOTAL=%.3e\n'], ...
             out, outer_irls, toc(t_start), ...
             cost_hist.TVx0(out), cost_hist.TVjoint(out), ...
             cost_hist.dxL1(out), cost_hist.total(out));
end

%% Final cost
dx_vec  = dx(valid_mask);
lam_vec = lambda_dx_eff(valid_mask);

reg0 = reg_MEDI(opts1, x0, @(z) z);
regJ = reg_MEDI(opts2, x0 + dx, @(z) z);

cost.data1   = norm(t1(:), 2);
cost.data2   = norm(t2(:), 2);
cost.TVx0    = reg0.cost(x0);
cost.TVjoint = regJ.cost(x0 + dx);
cost.dxL1    = sum(lam_vec(:) .* abs(dx_vec(:)));
cost.total   = cost.data1^2 + cost.data2^2 + cost.TVx0 + cost.TVjoint + cost.dxL1;

%% Convert to ppm
x0_ppm = (x0 / scale) .* MaskU;
dx_ppm = dx / scale;

if ~isempty(Mask_CSF)
    mu0 = mean(x0_ppm(Mask_CSF));
    mud = mean(dx_ppm(Mask_CSF));
    x0_ppm = x0_ppm - mu0;
    dx_ppm = dx_ppm - mud;
end

end

function value = get_option(vargs, name, default)
if nargin < 3
    default = [];
end

value = default;
n = numel(vargs);

for k = 1:2:n
    if (ischar(vargs{k}) || isstring(vargs{k})) && strcmpi(vargs{k}, name)
        if k + 1 <= n
            value = vargs{k + 1};
        else
            warning('Option "%s" has no value. Using default.', name);
        end
        return;
    end
end
end

function [x0, xd] = pcg_block(Aop, b0, bd, x0, xd, tol, maxit)
[Ap0, Apd] = Aop(x0, xd);
r0 = b0 - Ap0;
rd = bd - Apd;
p0 = r0;
pd = rd;

rr_old = dotc(r0, r0) + dotc(rd, rd);

for it = 1:maxit
    [Ap0, Apd] = Aop(p0, pd);
    denom = dotc(p0, Ap0) + dotc(pd, Apd);
    alpha = rr_old / max(denom, eps);

    x0 = x0 + alpha * p0;
    xd = xd + alpha * pd;

    r0 = r0 - alpha * Ap0;
    rd = rd - alpha * Apd;

    rr_new = dotc(r0, r0) + dotc(rd, rd);

    if sqrt(rr_new) < tol * sqrt(numel(x0) + numel(xd))
        break;
    end

    beta = rr_new / max(rr_old, eps);
    p0 = r0 + beta * p0;
    pd = rd + beta * pd;
    rr_old = rr_new;
end
end

function value = dotc(a, b)
value = real(sum(conj(a(:)) .* b(:)));
end