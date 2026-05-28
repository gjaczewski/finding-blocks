from pyscf import gto, scf, ao2mo, x2c
import numpy as np
relativistic = True
mol = gto.Mole()
mol.atom = 'H 0.0 0.0 0.0; H 0 0 0.74'
mol.basis = 'unc-sto-3g'
mol.spin = 0
mol.charge = 0
#mol.verbose = 4 # Żeby widzieć szczegóły iteracji
mol.build()
nao = mol.nao 
print(nao)
if relativistic:
    mf_scalar = scf.RHF(mol).x2c()
    mf_scalar.kernel()
    C_spatial = mf_scalar.mo_coeff 
    energies = mf_scalar.mo_energy
    mf_2c = scf.GHF(mol).x2c()
    hcore_ao_2c = mf_2c.get_hcore()
    C_spin = np.block([[C_spatial, np.zeros_like(C_spatial)],[np.zeros_like(C_spatial), C_spatial]])
    hcore_mo_2c = C_spin.conj().T @ hcore_ao_2c @ C_spin
    hcore_mo_2c = hcore_mo_2c.reshape(2, nao, 2, nao).transpose(1, 0, 3, 2)
    eri_spatial_chem = ao2mo.kernel(mol, C_spatial)
    eri_spatial_chem = ao2mo.restore(1, eri_spatial_chem, nao)
    eri_phys = eri_spatial_chem.transpose(0, 2, 1, 3)
    h_upup = hcore_mo_2c[:,0,:,0]
    h_dndn = hcore_mo_2c[:,1,:,1]
    h_updn = hcore_mo_2c[:,0,:,1]
    h_dnup = hcore_mo_2c[:,1,:,0]
else:
    kin = mol.intor('int1e_kin')
    vnuc = mol.intor('int1e_nuc')
    overlap = mol.intor('int1e_ovlp')
    eri = mol.intor('int2e')
    mf = scf.RHF(mol)
    mf.kernel()
    energies=mf.mo_energy
    C = mf.mo_coeff  
    h_ao = kin + vnuc + 0j
    h_mo = C.conj().T @ h_ao @ C 
    eri_mo = ao2mo.kernel(mol, C)
    eri_mo = ao2mo.restore(1, eri_mo, nao)
    eri_phys = eri_mo.transpose(0,2,1,3)
    h_upup = h_mo 
    h_dndn = h_mo
    h_updn = 0*h_mo
    h_dnup = 0*h_mo

real_h_upup = h_upup.real
imag_h_upup = h_upup.imag
out = np.column_stack((real_h_upup, imag_h_upup))
np.savetxt("integrals/h_upup.txt",out,delimiter=",")
with open('integrals/h_upup.txt', 'w') as f:
    for row in h_upup:
        line = ' '.join(f"({z.real},{z.imag})" for z in row)
        f.write(line + '\n')

real_h_dndn = h_dndn.real
imag_h_dndn = h_dndn.imag
out = np.column_stack((real_h_dndn, imag_h_dndn))
np.savetxt("integrals/h_dndn.txt",out,delimiter=",")
with open('integrals/h_dndn.txt', 'w') as f:
    for row in h_dndn:
        line = ' '.join(f"({z.real},{z.imag})" for z in row)
        f.write(line + '\n')

real_h_updn = h_updn.real
imag_h_updn = h_updn.imag
out = np.column_stack((real_h_updn, imag_h_updn))
np.savetxt("integrals/h_updn.txt",out,delimiter=",")
with open('integrals/h_updn.txt', 'w') as f:
    for row in h_updn:
        line = ' '.join(f"({z.real},{z.imag})" for z in row)
        f.write(line + '\n')

real_h_dnup = h_dnup.real
imag_h_dnup = h_dnup.imag
out = np.column_stack((real_h_dnup, imag_h_dnup))
np.savetxt("integrals/h_dnup.txt",out,delimiter=",")
with open('integrals/h_dnup.txt', 'w') as f:
    for row in h_dnup:
        line = ' '.join(f"({z.real},{z.imag})" for z in row)
        f.write(line + '\n')


np.savetxt("integrals/mo_energies.txt", energies,delimiter=",")
eri_flat=np.zeros([nao**2,nao**2])
p = 0
for i in range(0,nao):
    for j in range(0,nao):
        q = 0
        for k in range(0,nao):
            for  l in range(0,nao):
                eri_flat[p,q]=eri_phys[i,j,k,l]
                q = q + 1
        p = p + 1

        
eri_flat = eri_flat + 0j
real_eri_flat=eri_flat.real
imag_eri_flat=eri_flat.imag
out = np.column_stack((real_eri_flat, imag_eri_flat))
np.savetxt("integrals/eri_mo.txt",out,delimiter=",")
with open('integrals/eri_mo.txt', 'w') as f:
    for row in eri_flat:
        line = ' '.join(f"({z.real},{z.imag})" for z in row)
        f.write(line + '\n')