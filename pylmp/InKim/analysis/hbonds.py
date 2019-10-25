import numpy as np
import bgftools as bt

def pair_energy(donors_aNo, acceptors_aNo):


def hbonds(mybgf, selection = ""):
    '''
    This function calculates Hydrogen bonds(hbonds) between O(donor) and O(acceptor), especially for water molecules.
    Input: 
        - bgf.BgfFile mybgf
        - string selection: a region to search donor and acceptor atoms in mybgf
    Output:
        - int n_sel_atoms: number of hbond-able atoms in the selection region
        - list hbonds: self-descriptive
    '''

    # variables
    pbc = mybgf.CRYSTX[:3]
    d_crit = 3.5; a_crit = 30.0 # Chandler's criteria

    if selection:
        A = [atom for atom in mybgf.a if "O" in atom.ffType and eval(selection)]
        D = [atom for atom in mybgf.a if "O" in atom.ffType and eval(selection)]
    else:
        A = [atom for atom in mybgf.a if "O" in atom.ffType]
        D = [atom for atom in mybgf.a if "O" in atom.ffType]

    if not len(A) or not len(D):
        nu.warn("There are no atoms which can make hbonds (O atoms so far)!")
        return

    # calculate hbonds
    hbonds = []; 

    for d_atom in D:
        d = np.array([d_atom.x, d_atom.y, d_atom.z])    # donor coord
        neigh_anos = bt.get_neighbors_aNo(A, d, r=d_crit, pbc=pbc, k=7)
        donors = [d_atom.aNo] + d_atom.CONECT

        for ano in neigh_anos:
            a_atom = mybgf.getAtom(ano)
            a = np.array([a_atom.x, a_atom.y, a_atom.z])    # acceptor coord
            acceptors = [a_atom.aNo] + a_atom.CONECT

            for ano in d_atom.CONECT:
                h_atom = mybgf.getAtom(ano)
                h = np.array([h_atom.x, h_atom.y, h_atom.z])
                u = h - d; v = a - d; 
                theta = np.dot(u, v) / norm(u) / norm(v); theta = np.degrees(arccos(theta))
                if theta < a_crit:  # HBond exists
                    dist = nu.pbc_dist(a, d, pbc)
                    dist_ah = nu.pbc_dist(d, h, pbc)

                    # E_vdw
                    sigma_r = O_sigma / dist; sigma_r_6 = sigma_r**6; sigma_r_12 = sigma_r**12
                    E_vdw = 4.0 * O_epsilon * (sigma_r_12 - sigma_r_6); # E_vdw in kcal/mol

                    # E_coul
                    E_coul = 0.0
                    for i, j in itertools.product(donors, acceptors):
                        atom1 = mybgf.getAtom(i)
                        atom2 = mybgf.getAtom(j)
                        a1 = [atom1.x, atom1.y, atom1.z]
                        a2 = [atom2.x, atom2.y, atom2.z]
                        dist_ij = nu.pbc_dist(a1, a2, pbc)
                        E_coul += 332.06371 * atom1.charge * atom2.charge / dist_ij # E_coul in kcal/mol

                    # E_hbond
                    E_hbond = E_coul + E_vdw  # E_hbond = E_vdw + E_coul

                    # update for v4
                    # angle between H-O-H plane and O..O vector
                    #H1 = mybgf.getAtom(d_atom.CONECT[0]); h1 = [H1.x, H1.y, H1.z]    # H1
                    #H2 = mybgf.getAtom(d_atom.CONECT[1]); h2 = [H2.x, H2.y, H2.z]    # H2
                    #p = d - h1; q = d - h2; n = np.cross(p, q)  # normal vector of the H1-O-H2 plane
                    #m = a - d;  # O..O vector
                    #alpha = np.dot(n, m) / norm(n) / norm(m); alpha = np.degrees(arcsin(alpha)) # angle between H-O-H plane and O..O vector

                    #hbonds.append([d_atom.aNo, a_atom.aNo, d, a, dist, theta, [E_coul, E_vdw, E_hbond], dist_ah, alpha])   # v4
                    hbonds.append([d_atom.aNo, a_atom.aNo, dist, theta, E_hbond])   # v5

    return hbonds

