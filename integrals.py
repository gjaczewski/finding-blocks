from pyscf import gto, scf, ao2mo
import numpy as np
mol =gto.Mole()
mol.build(atom='H 0 0 0; H 0 0 1.1', basis='sto-3g')
kin = mol.intor('int1e_kin')
vnuc = mol.intor('int1e_nuc')
overlap = mol.intor('int1e_ovlp')
eri = mol.intor('int2e')

mf = scf.RHF(mol)
mf.kernel()

C = mf.mo_coeff   # macierz AO → MO

h_ao = kin + vnuc

h_mo = C.T @ h_ao @ C

eri_mo = ao2mo.kernel(mol, C)

norb = C.shape[1]
print(f"norb={norb}")

eri_mo = ao2mo.restore(1, eri_mo, norb)

eri_phys = eri_mo.transpose(0,2,1,3)
energies=mf.mo_energy
E_nuc = mol.energy_nuc()

np.savetxt("h_mo.txt", h_mo)
np.savetxt("mo_energies.txt", energies)
n = h_mo.shape[0]
eri_flat = eri_phys.reshape(n*n, n*n)
np.savetxt("eri_mo.txt", eri_flat)

print(eri_phys)
print(eri_flat)