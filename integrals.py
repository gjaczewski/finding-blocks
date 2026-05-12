from pyscf import gto, scf, ao2mo, x2c
import numpy as np
mol = gto.M(atom = 'O 0.000000 0.000000 0.000000',basis='sto-3g',spin=2,charge=0,verbose=0)
kin = mol.intor('int1e_kin')
vnuc = mol.intor('int1e_nuc')
overlap = mol.intor('int1e_ovlp')
eri = mol.intor('int2e')
relativistic = True
mf = scf.RHF(mol)
mf.kernel()
print(mol.nelectron)
print(mol.nao)

C = mf.mo_coeff   
n_mo=C.shape[1]
h_ao = kin + vnuc + 0j

h_mo = C.conj().T @ h_ao @ C 


if relativistic:
    mf_soc = scf.GHF(mol).x2c()
    hcore_x2c = mf_soc.get_hcore()
    n_ao = mol.nao_nr()
    h_aa = hcore_x2c[:n_ao, :n_ao]
    h_bb = hcore_x2c[n_ao:, n_ao:]
    h_ab = hcore_x2c[:n_ao, n_ao:]
    h_ba = hcore_x2c[n_ao:, :n_ao]
    h_soc_mo_ab = C.T @ h_ab @ C
    h_soc_mo_ba = C.T @ h_ba @ C
    h_soc_mo_aa = C.T @ h_aa @ C
    h_soc_mo_bb = C.T @ h_bb @ C


    h_upup = h_soc_mo_aa
    h_dndn = h_soc_mo_bb

    h_updn = h_soc_mo_ab
    h_dnup = h_soc_mo_ba
else:
    h_upup = h_mo 
    h_dndn = h_mo
    h_updn = 0*h_mo
    h_dnup = 0*h_mo

real_h_upup = h_upup.real
imag_h_upup = h_upup.imag
out = np.column_stack((real_h_upup, imag_h_upup))
np.savetxt("h_upup.txt",out,delimiter=",")
with open('h_upup.txt', 'w') as f:
    for row in h_upup:
        line = ' '.join(f"({z.real},{z.imag})" for z in row)
        f.write(line + '\n')

real_h_dndn = h_dndn.real
imag_h_dndn = h_dndn.imag
out = np.column_stack((real_h_dndn, imag_h_dndn))
np.savetxt("h_dndn.txt",out,delimiter=",")
with open('h_dndn.txt', 'w') as f:
    for row in h_dndn:
        line = ' '.join(f"({z.real},{z.imag})" for z in row)
        f.write(line + '\n')

real_h_updn = h_updn.real
imag_h_updn = h_updn.imag
out = np.column_stack((real_h_updn, imag_h_updn))
np.savetxt("h_updn.txt",out,delimiter=",")
with open('h_updn.txt', 'w') as f:
    for row in h_updn:
        line = ' '.join(f"({z.real},{z.imag})" for z in row)
        f.write(line + '\n')

real_h_dnup = h_dnup.real
imag_h_dnup = h_dnup.imag
out = np.column_stack((real_h_dnup, imag_h_dnup))
np.savetxt("h_dnup.txt",out,delimiter=",")
with open('h_dnup.txt', 'w') as f:
    for row in h_dnup:
        line = ' '.join(f"({z.real},{z.imag})" for z in row)
        f.write(line + '\n')

eri_mo = ao2mo.kernel(mol, C)

norb = C.shape[1]


eri_mo = ao2mo.restore(1, eri_mo, norb)

eri_phys = eri_mo.transpose(0,2,1,3)
energies=mf.mo_energy
E_nuc = mol.energy_nuc()


np.savetxt("mo_energies.txt", energies,delimiter=",")
n = h_mo.shape[0]
eri_flat=np.zeros([norb**2,norb**2])
p = 0
for i in range(0,norb):
    for j in range(0,norb):
        q = 0
        for k in range(0,norb):
            for  l in range(0,norb):
                eri_flat[p,q]=eri_phys[i,j,k,l]
                q = q + 1
        p = p + 1
eri_flat = eri_flat + 0j
real_eri_flat=eri_flat.real
imag_eri_flat=eri_flat.imag
out = np.column_stack((real_eri_flat, imag_eri_flat))
np.savetxt("eri_mo.txt",out,delimiter=",")
with open('eri_mo.txt', 'w') as f:
    for row in eri_flat:
        line = ' '.join(f"({z.real},{z.imag})" for z in row)
        f.write(line + '\n')
