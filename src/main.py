import math,os,sys
import scipy as sp

data_dir = './data/'
pfi_file = ''

class Pfi():
    def __init__(self,pfi_data_file):
        self.pfi_x, self.pfi_y, self.pfi_r = sp.genfromtxt(pfi_data_file,usecols=(5,6,7),unpack=True)
        self.pfi_x /= 3600.
        self.pfi_y /= 3600.
        self.pfi_r /= 3600.
        return None
    def get_num_cobra(self):
        return len(self.pfi_x)
    def map_on_sky(self,pc=(0.0,0.0),pa=30.0):
        self.pc = pc
        self.pa = pa
        self.pfi_x_rot = sp.cos(sp.pi/180.*self.pa) * self.pfi_x - sp.sin(sp.pi/180.*self.pa) * self.pfi_y + self.pc[0]
        self.pfi_y_rot = sp.sin(sp.pi/180.*self.pa) * self.pfi_x + sp.cos(sp.pi/180.*self.pa) * self.pfi_y + self.pc[1]
        return self.pfi_x_rot,self.pfi_y_rot

class Targets():
    def __init__(self):
        return None
    def get_tgt(self,target_file):
        dat = sp.genfromtxt(target_file,usecols=(2,3,4),dtype=[('id','S10'),('ra','f8'),('dec','f8')])
        self.tgt_id, self.tgt_x, self.tgt_y  = dat['id'], dat['ra'], dat['dec']
        return 0
    def get_tgt_id(self):
        return self.tgt_id
    def get_tgt_field_center(self):
        return sp.median(self.tgt_x),sp.median(self.tgt_y)
    def map_on_sky(self):
        return self.tgt_x, self.tgt_y

class Assigner():
    def __init__(self,pfi_x,pfi_y,tgt_x,tgt_x):
        return None
    def get_pfi_map(self):
        self.pfi_x = pfi_x
        self.pfi_y = pfi_y
        return self.pfi_x,self.pfi_y

def main():
    pfi = Pfi('./data/fiberpos.el90.dat')
    tgt = Targets()
    tgt.get_tgt('./data/dbsim_test_cosmology.dat')

    pfi.map_on_sky(pc=(0.0,0.0),pa=30.0)
    tgt.map_on_sky()

    assgn = Assigner(pfi.pfi_x_rot,pfi.pfi_y_rot,tgt.tgt_x,tgt.tgt_y)
    assgn.get_pfi_map()

#    print pfi.get_num_cobra()
#    for i in range(pfi.get_num_cobra()):
#        print pfi.pfi_x[i],pfi.pfi_y[i],pfi.pfi_x_rot[i],pfi.pfi_y_rot[i]
#    for i in range(len(tgt.tgt_id)):
#        print tgt.tgt_id[i],tgt.tgt_x[i],tgt.tgt_y[i]

if __name__ == '__main__':
    main()
