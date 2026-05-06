from pyscf import gto, scf, ao2mo
import numpy as np
mol =gto.Mole()
mol.build(atom='O 0.000000  0.000000  0.000000; H 0.000000 -0.757000  0.587000; H 0.000000  0.757000  0.587000', basis='sto-3g', charge=0, spin=0)
kin = mol.intor('int1e_kin')
vnuc = mol.intor('int1e_nuc')
overlap = mol.intor('int1e_ovlp')
eri = mol.intor('int2e')

mf = scf.RHF(mol)
mf.kernel()

C = mf.mo_coeff   

h_ao = kin + vnuc

h_mo = C.T @ h_ao @ C

eri_mo = ao2mo.kernel(mol, C)

norb = C.shape[1]
print(norb)

eri_mo = ao2mo.restore(1, eri_mo, norb)

eri_phys = eri_mo.transpose(0,2,1,3)
energies=mf.mo_energy
E_nuc = mol.energy_nuc()
h_mo = h_mo + 0j

real_h_mo = h_mo.real
imag_h_mo = h_mo.imag
out = np.column_stack((real_h_mo, imag_h_mo))
np.savetxt("h_mo.txt",out,delimiter=",")
with open('h_mo.txt', 'w') as f:
    for row in h_mo:
        line = ' '.join(f"({z.real},{z.imag})" for z in row)
        f.write(line + '\n')
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
