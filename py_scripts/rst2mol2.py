#!/bin/python
import sys
import numpy as np

class bead:
    ''' class bead contains all info about bead in chain: global number, 
    the remaining valence, type of the bead and coordinates '''
    def __init__(self):
        self.numbd = int()
        self.valence = int()
        self.typep = int()
        self.x = float()
        self.y = float()
        self.z = float()
        self.neighbors = []
        self.touched = bool()
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
    '''function to read the restart file to the class chain, it needs path to file and name of chain'''
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
            one.touched = False
            if one.typep==1:# or one.typep==2:
                polymer.bd.append(one)
        elif 'angles' in line:
            return
        elif bnd:
            sb=bond()
            head,tail = line.split()
            sb.first = min(int(head), int(tail))-1
            sb.last = max(int(head), int(tail))-1
            polymer.bnd.append(sb)
            
def removePBC (polymer):
    '''revome periodic boundary conditions in case of multiple chains'''
    itx=0
    ity=0
    itz=0
    buffer = [polymer.bnd[0]]
    polymer.bd[polymer.bnd[0].first].touched = True
    while len(buffer) > 0:
        for i in buffer:
            if polymer.bd[i.last].touched:
                i.first, i.last = i.last, i.first
            
            if polymer.bd[i.first].x - polymer.bd[i.last].x > polymer.xbox / 2:
                polymer.bd[i.last].x += polymer.xbox
            elif polymer.bd[i.first].x - polymer.bd[i.last].x < -polymer.xbox / 2:
                polymer.bd[i.last].x -= polymer.xbox

            if polymer.bd[i.first].y - polymer.bd[i.last].y > polymer.ybox / 2:
                polymer.bd[i.last].y += polymer.ybox
            elif polymer.bd[i.first].y - polymer.bd[i.last].y < -polymer.ybox / 2:
                polymer.bd[i.last].y -= polymer.ybox

            if polymer.bd[i.first].z - polymer.bd[i.last].z > polymer.zbox / 2:
                polymer.bd[i.last].z += polymer.zbox
            elif polymer.bd[i.first].z - polymer.bd[i.last].z < -polymer.zbox / 2:
                polymer.bd[i.last].z -= polymer.zbox
            
            polymer.bd[i.last].touched = True

            for j in polymer.bnd:
                if (j.first == i.last and (not polymer.bd[j.last].touched)) or (j.last == i.last and (not polymer.bd[j.first].touched)):
                    if (j not in buffer):
                        buffer.append(j)
            buffer.remove(i)

            
def writeMol2 (polymer, path):
    bstr='1  ala'
    the_file=open(path, 'w')
    the_file.write('@<TRIPOS>MOLECULE\n')
    the_file.write('mol_name\n')
    the_file.write('\t %d \t %d \t %s \t %s \t %s \n' %(len(polymer.bd), len(polymer.bnd), '0', '0', '0'))
    the_file.write('SMALL\n')
    the_file.write('USER_CHARGES\n')
    the_file.write('@<TRIPOS>ATOM\n')
    for i in range(len(polymer.bd)):
        ty='C'
        if polymer.bd[i].typep==2:
            ty='O'
        the_file.write('%d \t %s \t %f \t %f \t %f \t %s \t %s \t %f \n' %(i+1, ty, polymer.bd[i].x, polymer.bd[i].y, polymer.bd[i].z, ty, bstr, float(i)))
    the_file.write('@<TRIPOS>BOND\n')

    for i in range(len(polymer.bnd)):       
        the_file.write('%d \t %d \t %d \t %s \n' %(i+1, polymer.bnd[i].first+1, polymer.bnd[i].last+1, '1'))
    the_file.close()

def read_input():
    if len(sys.argv) < 5:
        sys.exit('Number of arguments isn\'t equal 4.')
    tr, arg1, arg2, arg3, arg4 = sys.argv
    if arg1=='-i':
        if arg3!='-o':
            sys.exit('-o option isn\'t found.')
        else:
            if 'mol2' in arg4:
                return arg2, arg4
            else:
                arg4 = arg4 + '.mol2'
                return arg2, arg4
    elif arg3=='-i':
        if arg1!='-o':
            sys.exit('-o option isn\'t found.')
        else:
            if 'mol2' in arg2:
                return arg4, arg2
            else:
                arg2 = arg2 + '.mol2'
                return arg4, arg2
    else:
        sys.exit('-i option isn\'t found.')


def main():
    inpath, outpath = read_input()
    poly=chain()
    f=open(inpath)
    read_rst(f, poly)
    poly.bd.sort()
    removePBC(poly)
    writeMol2(poly,outpath)
main()