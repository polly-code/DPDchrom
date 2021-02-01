import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

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
            if one.typep < 3:# or one.typep==2:
                polymer.bd.append(one)
        elif 'angles' in line:
            print ('Reading finished')
        elif bnd:
            sb=bond()
            head,tail=line.split()
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
            
def fill_distance_matr_chr4_200 (polymer, matrix, q):
    '''create distance matrix for chr4'''
    if q == 55:
        beg = 9291
    elif q == 25:
        beg = 11248
    elif q == 17:
        beg = 971
    if q == 48:
        beg = 5375
    #beg=9291#11248#971
    for i in range(beg, beg+764):
        for j in range(i+1,beg+764):
            matrix[i-beg][j-beg]=distance(polymer.bd[i], polymer.bd[j])
            matrix[j-beg][i-beg]=matrix[i-beg][j-beg]

def distance (a, b):
    '''calculate distance between two beads'''
    return np.sqrt((a.x-b.x)**2+(a.y-b.y)**2+(a.z-b.z)**2)

def calc_cr2 (a, b):
    '''calc critical index for two matrices'''
    return np.linalg.norm((a-b))/np.linalg.norm((a+b))

def check_imj_base_lvl():
'''calculate IMJ value for two symmetric random matrices'''
    imj=[]
    for _ in range(10):
        m1 = np.random.rand(1000,1000)
        m2 = np.random.rand(1000,1000)
    #    for i in range(10000):
    #        for j in range(i,10000):
    #            m1[i,j] = m1[j,i]
    #            m2[i,j] = m2[j,i]
        #imj.append(calc_cr2(m1,m2))
        print ((0.378-calc_cr2(m1, m1))/0.378*100)
    print(np.mean(imj), np.std(imj))

def calcl_contact_matrix(matrix, rcut=1.0):
    cm=np.zeros((len(matrix), len(matrix)))
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if matrix[i][j] < rcut:
                cm[i, j] = 1
    return cm

def shuffle_along_diag(m):
'''shuffle matrix along subdiagonal'''
    diag_sum_before = np.zeros(len(m))
    for i in range(len(m)):
        for j in range(i, len(m)):
            diag_sum_before[j-i] += m[i, j]
    
    for i in range(len(m)):
        arr=[]
        for j in range(0, len(m)-i):
            k=j+i
            arr.append(m[j, k])
        np.random.shuffle(arr)
        for j in range(0, len(m)-i):
            k=j+i
            m[j, k] = arr[j]
            m[k, j] = arr[j]
    
    diag_sum_after = np.zeros(len(m))
    for i in range(len(m)):
        for j in range(i, len(m)):
            diag_sum_after[j-i] += m[i, j]
    
    if (diag_sum_after.all() != diag_sum_before.all()):
        print ("Error! The sums of the subdiagonals are not equal!")

def main():
    df = pd.read_csv('path/to/NSN_combined_200kb.csv', sep=',', dtype={'chr1': str, 'chr2': str})
    df = df[df["chr1"]==df["chr2"]]
    df = df[df["chr1"]=='4']
    subtr_adapter = min(df['start1'].min(), df['start2'].min())
    df['start1'] = df['start1'] - subtr_adapter
    df['end1'] = df['end1'] - subtr_adapter
    df['start2'] = df['start2'] - subtr_adapter
    df['end2'] = df['end2'] - subtr_adapter
    a = (df['end1'].values + df['start1'].values)/400000
    b = (df['end2'].values + df['start2'].values)/400000
    c = df['count'].values
    cm_combined = np.zeros((764, 764))

    for it in range(len(a)):
        i=int(a[it]-0.5)
        j=int(b[it]-0.5)
        cm_combined[i,j] = c[it]#cm_combined[i,j] + 1
        cm_combined[j,i] = c[it]#cm_combined[j,i] + 1
    cm_combined = cm_combined / np.mean(cm_combined)
    for file1 in [48, 17, 25, 55]:
        path='path/to/'+str(file1)+'.dat'
        poly=chain()
        f=open(path)
        read_rst(f, poly)
        poly.bd.sort()
        removePBC(poly)
        rest = np.zeros((764,764))
        fill_distance_matr_chr4_200(poly, rest, file1)
        cmexp = calcl_contact_matrix(rest, rcut = 0.8)
        cmexp = cmexp / np.mean(cmexp)
        re = (0.378-calc_cr2(cm_combined, cmexp))/0.378*100
        print ('True', re)
        sh = []
        for _ in range(10):
            shuffle_along_diag(cmexp)
            sh.append((0.378-calc_cr2(cm_combined, cmexp))/0.378*100)
            #print ('Shuffled', sh)
        print(stats.ttest_ind([re], sh)[1])
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))
        bp1 = ax.boxplot([[re], sh], patch_artist=True, labels=['proper', 'shuffled'])
        colors = ['lightgreen', 'lightblue']
        for patch, color in zip(bp1['boxes'], colors):
            patch.set_facecolor(color)
        ax.set_title('Comparison of reconstructed contact map\nwith experimental one and\nwith shuffled ones for control')
        plt.savefig('path/to/'+str(file1)+'.pdf')
        plt.show()
        plt.clf()
        
def main2():
    reb = []
    sh = []
    for file1 in [48, 17, 25, 55]:
        path='path/to/'+str(file1)+'.dat'
        poly=chain()
        f=open(path)
        read_rst(f, poly)
        poly.bd.sort()
        removePBC(poly)
        rest = np.zeros((764,764))
        fill_distance_matr_chr4_200(poly, rest, file1)
        cmexp = calcl_contact_matrix(rest, rcut = 0.8)
        cmexp = cmexp / np.mean(cmexp)
        for file2 in [48, 17, 25, 55]:
            path2='path/to/'+str(file2)+'.dat'
            poly2=chain()
            f2=open(path2)
            read_rst(f2, poly2)
            poly2.bd.sort()
            removePBC(poly2)
            rest2 = np.zeros((764,764))
            fill_distance_matr_chr4_200(poly2, rest2, file2)
            cmexp2 = calcl_contact_matrix(rest2, rcut = 0.8)
            cmexp2 = cmexp2 / np.mean(cmexp2)
            
            re = (0.378-calc_cr2(cmexp2, cmexp))/0.378*100
            if file1 != file2:
                reb.append(re)
            print ('True', re)
            
            for _ in range(10):
                shuffle_along_diag(cmexp2)
                sh.append((0.378-calc_cr2(cmexp2, cmexp))/0.378*100)
            print ('Shuffled', sh[0])
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))
    
    bp1 = ax.boxplot([reb, sh], patch_artist=True, labels=['proper', 'shuffled'])
    print(stats.ttest_ind(reb, sh)[1])
    colors = ['lightgreen', 'lightblue']
    ax.set_title('Comparison of reconstructed contact map\nwith each other and\nwith shuffled ones for control')
    for patch, color in zip(bp1['boxes'], colors):
        patch.set_facecolor(color)
    plt.savefig('path/to/between_reconstructed.pdf')
    plt.show()

        