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

    def __init__(self, mol, fragment = True): 
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
        if not fragment: self.frags = []
        else: self.bfs()
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
            mols.append(crystal([mels,mcoords],fragment = False))
        self.frags = mols

    
 
a = crystal("sulfanilamide.xyz")
a.frags[0].write_out("TEST.xyz")
