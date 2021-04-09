import numpy as np
import random 
import sys

class bead:
    ''' class bead contains all info about bead in chain: global number, 
    the remaining valence, type of the bead and coordinates '''
    numbd=int()
    valence=int()
    typep=int()
    x=float()
    y=float()
    z=float()
    neighbors=[]
    def __lt__(self, other):
        return self.numbd < other.numbd
    
class bond:
    '''class bond contains all info about bond: which beads connected by this bond'''
    first=int()
    last=int()
    
class chain:
    '''class chain has two lists of beads and bonds and general info about system such as total number of particles, density and box size along each axis'''
    def __init__(self):
        self.bd=[]
        self.bnd=[]
        self.number_of_beads=int()
        self.number_of_bonds=int()
        self.density=float()
        self.xbox=float()
        self.ybox=float()
        self.zbox=float()

def read_rst (f, polymer):
    '''function to read the restart file to the class chain, get path to file and name of chain'''
    one=bead()
    sb=bond()
    bnd=False
    for i,line in enumerate(f):
        if i==0:
            head,tail=line.split()
            polymer.number_of_beads = int(head)
            polymer.density = float(tail)
        elif i==1:
            xbox, ybox, zbox=line.split()
            polymer.xbox=float(xbox)
            polymer.ybox=float(ybox)
            polymer.zbox=float(zbox)
        elif 'bonds:' in line:
            bnd=True
        elif i>1 & bnd==False:
            one=bead()
            numbd,valence,typep,x,y,z=line.split()
            one.numbd=int(numbd)
            one.valence=int(valence)
            one.typep=int(typep)
            one.x=float(x)
            one.y=float(y)
            one.z=float(z)
            if one.typep==1:# or one.typep==2:
                polymer.bd.append(one)
        elif 'angles' in line:
            print ('Reading finished')
        elif bnd:
            sb=bond()
            head,tail=line.split()
            sb.first=int(head)
            sb.last=int(tail)
            polymer.bnd.append(sb)

def readLmpRst (f, polymer):
    one=bead()
    sb=bond()
    flag_bead=False
    flag_bond=False
    for i,line in enumerate(f):
        if i>1:
            if not line.strip():
                continue
            if 'Atoms # bond' in line:
                flag_bond=False
                flag_bead=True
                continue
            if 'Bonds' in line:
                flag_bond=True
                flag_bead=False
                continue
            if 'Velocities' in line:
                flag_bond=False
                flag_bead=False
                continue
            if 'atoms' in line:
                head,tail = line.split()
                polymer.number_of_beads=int(head)
            elif 'bonds' in line:
                head,tail = line.split()
                polymer.number_of_bonds=int(head)
            elif 'xlo xhi' in line:
                q,w,e,r=line.split()
                polymer.xbox=float(w)
            elif 'ylo yhi' in line:
                q,w,e,r=line.split()
                polymer.ybox=float(w)
            elif 'zlo zhi' in line:
                q,w,e,r=line.split()
                polymer.zbox=float(w)
            elif flag_bead==True:
                one=bead()
                numbd,valence,typep,x,y,z,t1,t2,t3 = line.split()
                one.numbd=int(numbd)
                one.valence=int(valence)
                one.typep=int(typep)
                one.x=float(x)
                one.y=float(y)
                one.z=float(z)
                polymer.bd.append(one)
            elif flag_bond==True:
                sb=bond()
                num, typeb, head, tail=line.split()
                sb.first=int(head)
                sb.last=int(tail)
                polymer.bnd.append(sb)
				
def removePBC (polymer):
    '''revome periodic boundary conditions in case of single chain and numbers of beads corresspond to the global numbers and start from 1'''
    itx=0
    ity=0
    itz=0
    for i in range(len(polymer.bd)-1):
        if polymer.bd[i].x - itx * polymer.xbox - polymer.bd[i+1].x > polymer.xbox / 2:
            itx = itx + 1
            polymer.bd[i+1].x = polymer.bd[i+1].x + itx*polymer.xbox
        elif polymer.bd[i].x - itx * polymer.xbox - polymer.bd[i+1].x < -polymer.xbox / 2:
            itx = itx - 1
            polymer.bd[i+1].x = polymer.bd[i+1].x + itx*polymer.xbox
        else:
            polymer.bd[i+1].x = polymer.bd[i+1].x + itx*polymer.xbox
            
        if polymer.bd[i].y - ity * polymer.ybox - polymer.bd[i+1].y > polymer.ybox / 2:
            ity = ity + 1
            polymer.bd[i+1].y = polymer.bd[i+1].y + ity*polymer.ybox
        elif polymer.bd[i].y - ity * polymer.ybox - polymer.bd[i+1].y < -polymer.ybox / 2:
            ity = ity - 1
            polymer.bd[i+1].y = polymer.bd[i+1].y + ity*polymer.ybox
        else:
            polymer.bd[i+1].y = polymer.bd[i+1].y + ity*polymer.ybox
            
        if polymer.bd[i].z - itz * polymer.zbox - polymer.bd[i+1].z > polymer.zbox / 2:
            itz = itz + 1
            polymer.bd[i+1].z = polymer.bd[i+1].z + itz*polymer.zbox
        elif polymer.bd[i].z - itz * polymer.zbox - polymer.bd[i+1].z < -polymer.zbox / 2:
            itz = itz - 1
            polymer.bd[i+1].z = polymer.bd[i+1].z + itz*polymer.zbox
        else:
            polymer.bd[i+1].z = polymer.bd[i+1].z + itz * polymer.zbox

def distance (a, b):
    '''calculate distance between two beads'''
    return np.sqrt((a.x-b.x)**2+(a.y-b.y)**2+(a.z-b.z)**2)

def fill_contact_matr (polymer, matrix):
    '''create distance matrix for a system'''
    for i in range(len(polymer.bd)):
        for j in range(i,len(polymer.bd)):
            if distance(polymer.bd[i], polymer.bd[j])<0.7:              
                matrix[i][j]=1
                matrix[j][i]=1
            else:
                matrix[i][j]=0
                matrix[j][i]=0

def rst2ncc(matr, path):
    multiplicator=100000
    f=open(path, 'w')
    num=0
    for i in range(len(matr)):
        for j in range(i,len(matr)):
            if matr[i][j]==1:
                num+=1
                line='chr1\t' + str((i+1)*multiplicator-100)+'\t'+ str((i+1)*multiplicator+100)+'\t'+ str((i+1)*multiplicator-10)+'\t'+ str((i+1)*multiplicator+10)+'\t'+'+'+'\t'+'chr1\t' + str((j+1)*multiplicator-100)+'\t'+ str((j+1)*multiplicator+100)+'\t'+ str((j+1)*multiplicator-10)+'\t'+ str((j+1)*multiplicator+10)+'\t'+'+'+'\t' + str(num) + '\t' + str(123456)+'\t' + '0'
                f.write('%s\n'%line)
    f.close()
	
path=sys.argv[1]
poly=chain()
f=open(path)
if path[-4:] == '.dat':
	read_rst(f,poly)
else:
	readLmpRst(f,poly)
poly.bd.sort()
removePBC(poly)
matr=np.zeros((len(poly.bd),len(poly.bd)))
fill_contact_matr(poly,matr)
rst2ncc(matr, sys.argv[1] + '.ncc')
print ('One bead corresponds to 100Kb')