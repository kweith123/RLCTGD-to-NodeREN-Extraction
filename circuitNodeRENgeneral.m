%% =========================================================
% RLCTGD -> F -> resistor fold-in -> Hughes M -> Schur reduction -> (A,B,C,D)
%
% Iteratively reduces M to Hughes form (M13=M31=M33=0, M32=-M23^T)
% per the proof of Theorem 5 in Hughes 2017.
%
% Three reduction cases:
%   (a) (M33)_kk ~= 0          -> Schur out b-index k
%   (b) M33_diag = 0, off-diag -> Schur out two b-indices i,j
%   (M13) M33 = 0, M13 ~= 0    -> Schur out (port i, b-index j); flips i/o role of port i
%% =========================================================

clear; clc;

%% =========================================================
% Symbolic + dimensions
%% =========================================================
syms iq ip id v1 v2 i3 ia ib v4 v5 v6 real

nP  = 2;     nD  = 1;     nu  = nP + nD;

nLa = 0;  
nCa = 2;          % a-part storage (-> states)
nLb = 0; 
nCb = 1;          % b-part storage (candidate for elimination)

nxa = nLa + nCa;
nxb = nLb + nCb;
nx  = nxa + nxb;

nr  = 2;  ng  = 3;  nw  = nr + ng;

%% =========================================================
% F construction (KCL/KVL cutset-loop)
% xin order: [iLa ; vCa ; vLb ; iCb]
% x-block of zout: [vLa ; iCa ; iLb ; vCb]
%% =========================================================
uin = [iq ; ip ; id ];
xin = [v1 ; v2; i3];
win = [ia ; ib ; v4 ; v5; v6];
zin = [uin ; xin ; win];


i4 = iq; %Testing flipping port/diode currents
i5 = -ip;

i6 = id;

va = v1;
vb = v2;

i1 = iq - i3 - ia;
i2 = ip + i3 - ib - id;

vq = v1 + v4;
vp = v2 - v5;

v3 = v1 - v2;
vd = -v2 + v6; %FLIP THE DIODE POLARITY IN ZOUT!!!!!



zout = [vq ; vp ; vd ; -i1 ; -i2 ; -v3 ; -va ; -vb ; -i4 ; -i5 ; -i6];

F = jacobian(zout, zin);

iu = 1:nu;
ix = nu + (1:nx);
iw = nu + nx + (1:nw);

F11 = F(iu,iu); F12 = F(iu,ix); F13 = F(iu,iw);
F21 = F(ix,iu); F22 = F(ix,ix); F23 = F(ix,iw);
F31 = F(iw,iu); F32 = F(iw,ix); F33 = F(iw,iw);

%% =========================================================
% Fold w block: F_red on (e ; e1a, e1b) -> (y ; r1a, r1b)
%% =========================================================
if nr > 0, syms r [nr 1] real, Rblk = diag(r); else, Rblk = sym(zeros(0)); end
if ng > 0, syms g [ng 1] real, Gblk = diag(g); else, Gblk = sym(zeros(0)); end

Gamma = blkdiag(Rblk, Gblk);
Sigma_w = blkdiag(eye(nr), -eye(ng));
Gamma_s = simplify(Sigma_w * Gamma);


K  = simplify(inv(Gamma + F33));
Fa = simplify(F11 - F13*K*F31);
Fb = simplify(F12 - F13*K*F32);
Fc = simplify(F21 - F23*K*F31);
Fd = simplify(F22 - F23*K*F32);
F_red = [Fa Fb ; Fc Fd];

%% =========================================================
% Coordinate change to Hughes M
%% =========================================================
% a-storage: natural states


S1   = blkdiag(-eye(nLa), eye(nCa));

S2   = blkdiag(-eye(nLb), eye(nCb));

%Solve F_red,H using permuatation P
%% =========================================================
% General b-storage permutation to Hughes order
%% =========================================================

if nxb > 0
    if nLb == 0 || nCb == 0
        Pb = sym(eye(nxb));
    else
        Pb = [sym(zeros(nLb,nCb)), sym(eye(nLb));
              sym(eye(nCb)),       sym(zeros(nCb,nLb))];
    end

    P_H = blkdiag(sym(eye(nu)), sym(eye(nxa)), Pb);

    F_red = simplify(P_H * F_red * P_H.');
end



T = blkdiag(eye(nu), S1, -S2);
M = simplify(T.' * F_red * T);

%% =========================================================
% Storage matrices
%% =========================================================
if nLa>0, syms La [nLa 1] real positive; LamLa = diag(La); else, LamLa = sym(zeros(0)); end
if nCa>0, syms Ca [nCa 1] real positive; LamCa = diag(Ca); else, LamCa = sym(zeros(0)); end
if nLb>0, syms Lb [nLb 1] real positive; LamLb = diag(Lb); else, LamLb = sym(zeros(0)); end
if nCb>0, syms Cb [nCb 1] real positive; LamCb = diag(Cb); else, LamCb = sym(zeros(0)); end
Lam1 = blkdiag(LamLa, LamCa);
Lam2 = blkdiag(LamLb, LamCb);
% Override for the existing test circuit:
if nLa==1 && nCa==1 && nxb==0
    Lam1 = blkdiag(L*eye(nLa), C*eye(nCa));
end

%% =========================================================
% ---- ITERATIVE SCHUR REDUCTION TO HUGHES FORM ----
%% =========================================================
[M, S2, Lam2, port_swaps, log] = hughes_reduce_keep_b_storage(M, nu, nxa, S2, Lam2);

fprintf('Reduction log:\n');
for s = 1:numel(log), fprintf('  %s\n', log{s}); end
if ~isempty(port_swaps)
    warning(['M13 trick was invoked at port(s) [%s]: input/output roles ', ...
             'flipped at those ports. Re-classify u/y/z/zt accordingly.'], ...
             num2str(port_swaps));
end

disp('Reduced M ='); disp(simplify(M));

%% =========================================================
% Partition reduced M
%% =========================================================
ia = nu + (1:nxa);
ib = nu + nxa + (1:nxb);

M11 = M(1:nu , 1:nu);
M12 = M(1:nu , ia);
M21 = M(ia   , 1:nu);
M22 = M(ia   , ia);
M23 = M(ia   , ib);

%% =========================================================
% Hughes-hatted blocks + Omega
%% =========================================================
M12h = simplify(M12 * S1);
M21h = simplify(S1  * M21);
M22h = simplify(S1  * M22 * S1);
if nxb > 0
    M23h = simplify(S1 * M23 * S2);
else
    M23h = sym(zeros(nxa,0));
end
Omega = simplify(Lam1 + M23h * Lam2 * M23h.');

%% =========================================================
% Theorem 5 state-space
%% =========================================================
A    = simplify(-Omega \ M22h);
B    = simplify(-Omega \ M21h);
Cmat = simplify(M12h);
Dmat = simplify(M11);

% Algebraic e1b realisation (only when nxb > 0)
if nxb > 0
    E_u = simplify( Lam2 * M23h.' * (Omega \ M21h) );
    E_x = simplify( Lam2 * M23h.' * (Omega \ M22h) );
    disp('e1b = E_u*u + E_x*x'); disp(E_u); disp(E_x);
end

disp('Omega ='); disp(Omega);
disp('A =');     disp(A);
disp('B =');     disp(B);
disp('C =');     disp(Cmat);
disp('D =');     disp(Dmat);

%% =========================================================
% Diode split + NodeREN compact form
%% =========================================================
nu_ext = nP;  nz = nD;  ny = nP;

Bu  = B(:, 1:nu_ext);          Bz  = B(:, nu_ext+1:end);
Cu  = Cmat(1:ny, :);           Cz  = Cmat(ny+1:end, :);
Duu = Dmat(1:ny, 1:nu_ext);    Duz = Dmat(1:ny, nu_ext+1:end);
Dzu = Dmat(ny+1:end, 1:nu_ext);Dzz = Dmat(ny+1:end, nu_ext+1:end);

B1  = Bz;   B2  = Bu;
C1  = Cz;   C2  = Cu;
D11 = Dzz;  D12 = Dzu;
D21 = Duz;  D22 = Duu;

Mcompact = [A   B1   B2;
            C1  D11  D12;
            C2  D21  D22];

disp('Mcompact = [A B1 B2 ; C1 D11 D12 ; C2 D21 D22] ='); disp(Mcompact);


%% =========================================================
% ===================== HELPERS ===========================
%% =========================================================

function [M, S2, Lam2, port_swaps, log] = hughes_reduce_keep_b_storage(M, nu, nxa, S2, Lam2)
% Converts raw M into Hughes hatted form while KEEPING b-storage alive.
%
% Important:
%   This function does NOT delete e1b/r1b or drop rows/cols of Lam2.
%   It only removes the direct feedthrough blocks M13, M31, M33
%   by an output-side Hughes transformation.
%
% Starting form:
%
%   [ r    ]   [ M11 M12 M13 ] [ e    ]
%   [ r1a  ] = [ M21 M22 M23 ] [ e1a  ]
%   [ r1b  ]   [ M31 M32 M33 ] [ e1b  ]
%
% If M33 is invertible, define transformed outputs:
%
%   r_h   = r   - M13 M33^{-1} r1b
%   r1a_h = r1a
%   r1b_h = r1b - M31 e - M33 e1b
%
% This gives:
%
%   [ r_h   ]   [ M11h M12h 0 ] [ e    ]
%   [ r1a_h ] = [ M21  M22  M23 ] [ e1a ]
%   [ r1b_h ]   [ 0    M32  0 ] [ e1b ]
%
% Therefore M23 and Lam2 survive, so
%
%   Omega = Lam1 + M23h Lam2 M23h'
%
% still contains the b-storage contribution.

    port_swaps = [];
    log        = {};

    nxb = size(S2,1);
    if nxb == 0
        return;
    end

    ie = 1:nu;
    ia = nu + (1:nxa);
    ib = nu + nxa + (1:nxb);

    M11 = M(ie, ie);  M12 = M(ie, ia);  M13 = M(ie, ib);
    M21 = M(ia, ie);  M22 = M(ia, ia);  M23 = M(ia, ib);
    M31 = M(ib, ie);  M32 = M(ib, ia);  M33 = M(ib, ib);

    % Already in Hughes form.
    if is_zero_block(M13) && is_zero_block(M31) && is_zero_block(M33)
        log{end+1} = 'M already in Hughes zero-block form.';
        return;
    end

    % Case: M33 is nonzero and invertible.
    % Do NOT delete b-storage. Canonicalise the blocks instead.
    if ~is_zero_block(M33)
        if is_sym_zero(det(M33))
            error(['M33 is singular. Do not Schur-drop b-storage. ', ...
                   'Repartition a/b storage variables or implement the full ', ...
                   'unimodular polynomial-ring Hughes reduction.']);
        end

        M11h = simplify(M11 - M13 * (M33 \ M31));
        M12h = simplify(M12 - M13 * (M33 \ M32));
        M13h = sym(zeros(size(M13)));

        M21h = M21;
        M22h = M22;
        M23h = M23;

        M31h = sym(zeros(size(M31)));
        M32h = M32;
        M33h = sym(zeros(size(M33)));

        M = simplify([
            M11h, M12h, M13h;
            M21h, M22h, M23h;
            M31h, M32h, M33h
        ]);

        log{end+1} = 'canonicalised nonzero M33 while keeping S2 and Lam2.';
        return;
    end

    % If M33 = 0 but M13/M31 are nonzero, this is not safe to fix by
    % deleting b-storage. It requires the Hughes port-swap/unimodular step.
    if ~is_zero_block(M13) || ~is_zero_block(M31)
        error(['M33 = 0 but M13 or M31 is nonzero. ', ...
               'Do not drop Lam2. This requires a proper Hughes port-swap ', ...
               'or unimodular transformation, not the old schur_step/drop_idx.']);
    end
end

function M_new = schur_step(M, elim)
    n    = size(M,1);
    keep = setdiff(1:n, elim);
    A    = M(keep, keep);
    B    = M(keep, elim);
    C    = M(elim, keep);
    D    = M(elim, elim);
    M_new = simplify(A - B * (D \ C));
end

function [S2, Lam2] = drop_idx(S2, Lam2, ks)
    ks = sort(ks, 'descend');
    for k = ks
        S2(k,:)   = []; S2(:,k)   = [];
        Lam2(k,:) = []; Lam2(:,k) = [];
    end
end

function tf = is_sym_zero(x)
    tf = isequal(simplify(x), sym(0));
end

function tf = is_zero_block(X)
    if isempty(X), tf = true; return; end
    tf = isequal(simplify(X), sym(zeros(size(X))));
end

syms h;

Rh = inv(eye(size(A)) - h*A);
Hh = simplify(D11 + h*C1*Rh*B1);

disp('Hh =');
disp(Hh);

disp('Off diagonal coupling =');
disp(simplify(Hh - diag(diag(Hh))));