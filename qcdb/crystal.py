"""
This class will include: 

1) generation of unique n-mers 
2) generation of crystal from cif with distance parameters

Inputs can be: 

1) An xyz file of a crystal
2) A cif file for a crystal

Expects Angstrom input
"""
import numpy as np
#generic functions defined outside of class
cov_rad = {   'H' : 0.37, 'C' : 0.77, 'O' : 0.73, 'N' : 0.75, 'F' : 0.71,
'P' : 1.10, 'S' : 1.03, 'Cl': 0.99, 'Br': 1.14, 'I' : 1.33, 'He': 0.30,
'Ne': 0.84, 'Ar': 1.00, 'Li': 1.02, 'Be': 0.27, 'B' : 0.88, 'Na': 1.02,
'Mg': 0.72, 'Al': 1.30, 'Si': 1.18, 'K' : 1.38, 'Ca': 1.00, 'Sc': 0.75,
'Ti': 0.86, 'V' : 0.79, 'Cr': 0.73, 'Mn': 0.67, 'Fe': 0.61, 'Co': 0.64,
'Ni': 0.55, 'Cu': 0.46, 'Zn': 0.60, 'Ga': 1.22, 'Ge': 1.22, 'As': 1.22,
'Se': 1.17, 'Br': 1.14, 'Kr': 1.03}

temp_symbol = ["X", "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE", "NA", "MG",
"AL", "SI", "P", "S", "CL", "AR", "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO",
"NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y", "ZR", "NB",
"MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I", "XE", "CS",
"BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM",
"YB", "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI",
"PO", "AT", "RN", "FR", "RA", "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK",
"CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG",
"UUB", "UUT", "UUQ", "UUP", "UUH", "UUS", "UUO"]

z_num = {}
for i,el in enumerate(temp_symbol): z_num[el] = i

def dist(a,b):
    x = (a[0]-b[0])**2
    y = (a[1]-b[1])**2
    z = (a[2]-b[2])**2
    dist = (x+y+z)**0.5
    return dist

def bound(a_el,a,b_el,b):
    if cov_rad[a_el]+cov_rad[b_el] >=  dist(a,b): return True
    else: return False

class crystal():

    def __init__(self, mol): 
        if type(mol) is str and mol.endswith('.xyz'):
            self.fil = mol
            with open(mol,'r') as ofil:
                n = int(next(ofil).split('\n')[0])
                next(ofil)
                els = []
                geom = np.array([]).reshape(0,3)
                for l in range(n):
                    lin = next(ofil)
                    ls = lin.split()
                    els.append(ls[0])
                    coords = [float(x) for x in ls[1:]]
                    geom = np.vstack([geom,coords])
                mol = [els,geom]
        self.mol = mol
        self.frags = []

    def print_out(self):
        els = self.mol[0]
        print(len(els))
        print("In Angstrom")
        coords = self.mol[1]
        for i in range(len(els)): 
            el = els[i]
            c = coords[i]
            print("%2s %7.3f %7.3f %7.3f" % (el,c[0],c[1],c[2]))
        print("")

    def write_out(self,fil_name):
        with open(fil_name,'w') as fil:
            els = self.mol[0]
            fil.write("%d\n" % (len(els)))
            fil.write("In Angstrom\n")
            coords = self.mol[1]
            for i in range(len(els)): 
                el = els[i]
                c = coords[i]
                fil.write("%2s %7.3f %7.3f %7.3f\n" % (el,c[0],c[1],c[2]))

    def bfs(self):
        els = self.mol[0]
        coords = self.mol[1]
        lc = len(coords)
        bonds = []
        #find all bonded atoms
        for i in range(1,lc):
            for j in range(i):
                if not bound(els[i],coords[i],els[j],coords[j]): continue
                bonds.append(set([i,j]))
        #Construct subsets such that each pair has a null intersection (BFS)
        minds = []
        while bonds != []:
            mol = bonds[0]
            bonds.remove(bonds[0])
            intersect = True
            while intersect: 
                intersect = False
                remove = []
                for i in bonds:
                    if i & mol == set([]): continue
                    for j in i: mol.add(j)
                    intersect = True
                    remove.append(i)
                for i in remove: bonds.remove(i)
            minds.append(mol)
        #build mol out of separated inds
        mols = []
        for mol in minds:
            mels = []
            mcoords = []
            for i in mol:
                mels.append(els[i])
                mcoords.append(coords[i])
            mols.append(crystal([mels,mcoords]))
        self.frags = mols

    #function that moves mol to align with self
    # 0 is Q
    # 1 is P
    def rmsd(self, mol):
        #self
        els0 = self.mol[0]  
        coords0 = self.mol[1] 
        #shift to center
        coords0 = coords0 - coords0.mean(axis=0)
        #mol
        els1 = mol.mol[0]  
        coords1 = mol.mol[1]  
        coords1 = coords1 - coords1.mean(axis=0)
        
        #covariance matrix
        H = np.dot(coords1.T, coords0)
        V, s, Wt = np.linalg.svd(H)
        d = np.linalg.det(np.dot(V,Wt))
        I = np.identity(3)
        I[2,2] = d 
        R = np.dot(np.dot(V,I),Wt)
        coords1 = np.dot(coords1,R)
        ans = np.linalg.norm(coords0 - coords1) / np.sqrt(coords0.shape[0])
        return ans
    
    def extract_frags(self, frag_nums):
        els = []
        coords = np.array([]).reshape(0,3)
        for f in frag_nums:
            els = els + self.frags[f].mol[0]    
            coords = np.vstack([coords,self.frags[f].mol[1]])
        return crystal([els,coords])    

    def nuclear_repulsion_energy(self):
        e = 0.
        for a1 in range(len(self.mol[0])):
            for a2 in range(a1):
                z1 = z_num[self.mol[0][a1]]
                z2 = z_num[self.mol[0][a2]]
                d = 1.88973*dist(self.mol[1][a1],self.mol[1][a2])
                e += z1 * z2 / d
        return e 

    def distance_matrix(self):
        natoms = len(self.mol[0])
        mat = np.zeros((natoms,natoms))
        for a1 in range(natoms):
            for a2 in range(a1):
                mat[a1][a2] = dist(self.mol[1][a1],self.mol[1][a2])
        return mat

    def get_nmers(self, N):
        """ Returns a dictionary of lists that contain crystal objects for
            all unique n-mers from n=1-N. Key values correspond to the value of 
            n.

        >>> H2O_20.get_nmers(1)
        {1: [<__main__.crystal object at 0x7fc5caa3b4e0>, <__main__.crystal object at 0x7fc5caa3b3c8>]}
        """

        from itertools import combinations
        
        self.bfs()
        
        nfrags = len(self.frags)
        if N > nfrags: 
            raise Exception("BFS only found %d fragments. Cannot find %d-mer with %d fragments." % (nfrags,N,nfrags)) 

        inds = [x for x in range(nfrags)]
        
        #dictionary with lists for each class of n-mer
        unique = {}
        for i in range(N): unique[i+1] = []

        #for each type of n-mer
        for i in range(N):
            #add each unique n-mer
            nres = {}
            for combo in combinations(inds,i+1):
                mol = self.extract_frags(combo)
                nre = mol.nuclear_repulsion_energy()
                #tweak to whatever tightness we'd likee
                nre = round(nre,2)                    
                try: nres[nre].append(mol)
                except: nres[nre] = [mol]
            for nre in nres.keys():
                lnre = len(nres[nre])
                #last element not checked for duplication 
                unique[i+1].append(nres[nre][-1])
                #nothing to compare too if only one 
                if lnre == 1: continue
                for m1 in range(lnre-1):
                    duplicate = False
                    for m2 in range(m1+1,lnre):
                        mol1 = nres[nre][m1] 
                        mol2 = nres[nre][m2]
                        diff = mol1.rmsd(mol2) 
                        #make keyword for this option
                        if diff < 0.01: 
                            duplicate = True
                            break
                    if not duplicate: unique[i+1].append(mol1)
        return unique 
 
sulf1 = crystal("s7-1_Angstroms.xyz")
sulf2 = crystal("s7-2_Angstroms.xyz")

s1dm = sulf1.distance_matrix()
s2dm = sulf2.distance_matrix()
print(np.allclose(s1dm,s2dm, atol = 1e-2))
print(sulf1.rmsd(sulf2))
print(sulf1.nuclear_repulsion_energy())
print(sulf2.nuclear_repulsion_energy())


#a = crystal("138K_SMALL_FINAL.xyz")
#a.bfs()
#for mon in a.frags: print(mon.distance_matrix())
#for i,mon in enumerate(monomers): mon.write_out("mon_"+str(i)+".xyz")
 
#a = crystal("sulfanilamide.xyz")
#a.bfs()
#frags = a.frags
#print(frags[0].distance_matrix(), frags[1].distance_matrix())
#monomers = a.get_nmers(1)[1]
#m1 = monomers[0]
#m2 = monomers[1]
#print(m1.align(m2))
#m1.write_out("sulf1.xyz")
#m2.write_out("sulf2.xyz")
#print(a.rmsd(a))
#benz0 = crystal("benz0.xyz")
#benz1 = crystal("benz1.xyz")

#b0dm = benz0.distance_matrix()
#b1dm = benz1.distance_matrix()

#print(b0dm,b1dm)

#mon0 = crystal("mon_0.xyz")
#mon1 = crystal("mon_1.xyz")
#mon2 = crystal("mon_2.xyz")
#mon3 = crystal("mon_3.xyz")
#mons = [mon0, mon1, mon2, mon3]
#print(mons[0].rmsd(mons[1]))
##for m in range(1,len(mons)):
##    for n in range(m):
##        print(m,n,mons[m].align(mons[n]))
