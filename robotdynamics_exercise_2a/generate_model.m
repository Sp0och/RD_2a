%% Define Generalized Coordintates & Model Parameters

% Initialize generalized coordinates
gen_cor = generate_gen_cor;

% Model parameters
params = generate_params;

%% Generate Forward Kinematics

kin = generate_kin(gen_cor);

%% Generate Forward Differential Kinematics

jac = generate_jac(gen_cor, kin, params);

%% Generate Equations of Motion

eom = generate_eom(gen_cor, kin, params, jac);

%%
% EOF

function eom = generate_eom(gen_cor,kin,dyn,jac)

%%setup
phi = gen_cor.phi;
dphi = gen_cor.dphi;

T_Ik = kin.T_Ik;
R_Ik = kin.R_Ik;

k_I_s = dyn.k_I_s;
m = dyn.m;
I_g_acc = dyn.I_g_acc;
k_rk_s = dyn.k_rk_s;

I_Jp_s = jac.I_Jp_s;
I_Jr = jac.I_Jr_s;

eom.M = system(zeros(6,6));
eom.b = system(zeros(6,1));
eom.g = system(zeros(6,1));
eom.hamiltonian = system(zeros(1,1));


%%
fprintf('Computing Mass Matrix')
M = sym(zeros(6,6));
for k = 1:length(phi)
   M = M + I_Jp_s{k}' * m{k} * I_Jp_s{k} ...
       + I_Jr{k}' * R_Ik{k} * k_I_s{k} * R_Ik{k}' * I_Jr{k};
end

fprintf('Simplifying for Performance')
for k = 1:length(phi)
    for h = k:length(phi)
        m_kh = simplify(M(k,h));
        if(h == k)
            M(h,k) = simplify(m_kh);
        else
            M(h,k) = simplify(m_kh);
            M(k,h) = simplify(m_kh);
        end
    end
end
fprintf('done\n')


%gravity part
fprintf('calculating gravity vector')
g = system(6,1);
for k = 1:length(phi)
    I_F_gk = m{k}*I_g_acc;
    g = g - I_Jp_s{k}' * I_F_gk;
end
fprintf('simplifying...')
g = simplify(g);
fprintf('done\n')



fprintf('calculating b')
b = system(6,1);
for k = 1:length(phi)
    b = b + I_Jp_s{k}' * m{k} * simplify(dAdt(I_Jp_s{k},phi,dphi)) *dphi...
        + I_Jr{k}' * R_Ik{k} * k_I_s{k} * R_Ik{k}' * simplify(dAdt(I_Jr{k},phi,dphi)) * dphi...
        + cross(I_Jr{k}' * I_Jr{k}*dphi, R_Ik{k} * k_I_s{k} * R_Ik{k}' * I_Jr{k}*dphi);
end
fprintf('done\n')


fprintf('computing total energy')
kin_energy = 0.5*dphi'*M*dphi;
pot_energy = sym(0);
for i = 1:length(phi)
   pot_energy  = pot_energy - [eye(3),zeros(3,1)] * T_Ik{k} * [k_rk_s{k},1] * m{k} * I_g_acc;
end
hamiltonian = kin_energy + pot_energy;
fprintf('simplifying...')
hamiltonian = simplify(hamiltonian);

fname = mfilename;
fpath = mfilename('fullpath');dpath = strrep(fpath, fname, '');
fprintf('Generating eom scripts... ');
fprintf('M... ');
matlabFunction(M, 'vars', fphig, 'file', strcat(dpath,'/M fun'));
fprintf('g... ');
matlabFunction(g, 'vars', fphig, 'file', strcat(dpath,'/g fun'));
fprintf('b... ');
matlabFunction(b, 'vars', fphi, dphig, 'file', strcat(dpath,'/b fun'));
fprintf('hamiltonian... ');
matlabFunction(hamiltonian, 'vars', fphi, dphig, 'file', ...
strcat(dpath,'/hamiltonian_fun'));
fprintf('done!nn');

eom.M = M;
eom.b = b;
eom.g = g;
eom.hamiltonian = hamiltonian;
eom.enPot = enPot;
eom.enKin = enKin;

end