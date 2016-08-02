import math,os,sys
import scipy as sp
import scipy.interpolate as it
import subprocess

home = '/work/PFS/pfs_obs_sim/'
data_dir = home + 'data/'
out_dir = home + 'out/'
ets_dir = home + 'ets_fiber_assigner/'
etc_dir = home + 'spt_ExposureTimeCalculator/'
wave_norm = 480.0

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
    def __init__(self,target_file):
        dat = sp.genfromtxt(target_file,usecols=(2,3,4,17,7,11),dtype=[('id','i4'),('ra','f8'),('dec','f8'),('texp','f8'),('priority','f8'),('mag','f8')])
        self.tgt_id, self.tgt_x, self.tgt_y, self.tgt_texp, self.tgt_prio, self.tgt_mag = dat['id'], dat['ra'], dat['dec'], dat['texp'], dat['priority'], dat['mag']
        self.tgt_z = sp.random.uniform(0.5,2.5,len(self.tgt_id))
        return None
    def gen_list_for_ets(self,filename):
        file = open(filename,'w')
        for i in range(len(self.tgt_id)):
            file.write('%10s   %9.5f  %9.5f  %7.1f  %2d  %6.3f  N/A  N/A\n' % ('ID%d'%(self.tgt_id[i]),self.tgt_x[i],self.tgt_y[i],self.tgt_texp[i],self.tgt_prio[i],self.tgt_mag[i]))
        file.close()
        return 0
    def get_tgt_dict(self):
        self.dict_tgt_x = {}
        self.dict_tgt_y = {}
        self.dict_tgt_mag = {}
        self.dict_tgt_z = {}
        for i in range(len(self.tgt_id)):
            self.dict_tgt_x[self.tgt_id[i]] = self.tgt_x[i]
            self.dict_tgt_y[self.tgt_id[i]] = self.tgt_y[i]
            self.dict_tgt_mag[self.tgt_id[i]] = self.tgt_mag[i]
            self.dict_tgt_z[self.tgt_id[i]] = self.tgt_z[i]
        return 0

    def get_tgt_id(self):
        return self.tgt_id
    def get_tgt_field_center(self):
        return sp.median(self.tgt_x),sp.median(self.tgt_y)
    def map_on_sky(self):
        return self.tgt_x, self.tgt_y

class Template():
    def __init__(self,filename):
        self.x, self.y = sp.genfromtxt(filename,unpack=True,usecols=(0,1))
        intrp = it.interp1d(self.x,self.y)
        self.mag_norm_temp = intrp(wave_norm)
        return None
    def get_norm_template(self,z,mag,outfile):
        self.z = z
        self.mag = mag
        self.x_norm = self.x * (1+self.z)
        self.y_norm = self.y + (self.mag - self.mag_norm_temp)
        file = open(outfile,'w')
        for i in range(len(self.x)):
            file.write('%.4e %.4e\n' % (self.x_norm[i],self.y_norm[i]))
        file.close()
        return 0

class Ets():
    def __init__(self):
        cmd = 'cd %s; make' % (ets_dir)
        subprocess.call(cmd, shell=True)
        return None
    def run_assigner(self,assigner):
        self.assigner = assigner
        self.fraction = '0.50'
        self.objlist = '%s/ets.list' % (out_dir)
        self.output = '%s/ets.out' % (out_dir)
        cmd = '%sets_demo assigner=%s input=%s fract=%s output=%s' % (ets_dir,self.assigner,self.objlist,self.fraction,self.output)
        proc = subprocess.call(cmd, shell=True)
        if proc == 0:
            print "Fiber Allocation Done!"
            C=0
            self.assigned_tgt = []
            self.assigned_fib = []
            self.assigned_tgt_ra = []
            self.assigned_tgt_dec = []
            for line in open(self.output,'r'):
                if C<=1:
                    a = line.split()
                    if a[0] == "Exposure":
                        C+=1
                    if a[0] != "Exposure" and a[0] != "Target":
                        self.assigned_tgt.append(int(a[0]))
                        self.assigned_fib.append(int(a[1]))
                        self.assigned_tgt_ra.append(float(a[2]))
                        self.assigned_tgt_dec.append(float(a[3]))
        else:
            print "Fiber Allocation Failed!"
        return 0

class Etc():
    def __init__(self):
        cmd = 'cd %s; make' % (etc_dir)
        proc = subprocess.call(cmd, shell=True)
        return None
    def get_sn_table(self):
        cmd = 'cd %s; python run_etc.py @run_etc.defaults --OUTFILE_NOISE=out/ref.noise.dat --OUTFILE_SNC=out/ref.snc.dat --OUTFILE_SNL=- --EXP_TIME=450 --EXP_NUM=2 --MAG_FILE=23.0 ' % (etc_dir)
        proc = subprocess.call(cmd, shell=True)
        if proc == 0:
            print "SN table is created!"
        return 0
    def get_simulated_spectra(self,mag_file,objId):
        cmd = 'cd %s; python gen_sim_spec.py @gen_sim_spec.defaults --EXP_NUM=2 --MAG_FILE=%s --etcFile=out/ref.snc.dat --nrealize=1 --outDir=%s --asciiTable=None --writeFits=True --objId=%d' % (etc_dir,mag_file,out_dir,objId)
        proc = subprocess.call(cmd, shell=True)
        return proc

def main():
    if os.path.exists(out_dir)==False:
        os.mkdir(out_dir)
#    pfi.map_on_sky(pc=(0.0,0.0),pa=30.0)
#    tgt.map_on_sky()
## taget list ##
    print "loading target list ..."
    tgt = Targets('%s/catalogue/dbsim_test_cosmology.dat' % (data_dir))
    tgt.gen_list_for_ets('%s/ets.list' % (out_dir))
    tgt.get_tgt_dict()
## ETS package ##
    print "initializing ETS package ..."
    ets = Ets()
    ets.run_assigner(assigner='naive')
    print "%d objects are assigned!" % (len(ets.assigned_tgt))
    file = open('test.dat','w')
    for i in range(len(ets.assigned_tgt)):
        file.write('%d %.5f %.5f %.5f %.5f\n' % (ets.assigned_tgt[i],ets.assigned_tgt_ra[i],ets.assigned_tgt_dec[i],tgt.dict_tgt_x[ets.assigned_tgt[i]],tgt.dict_tgt_y[ets.assigned_tgt[i]]))
    file.close()
## spectral template ##
    stmp = Template(data_dir+'/spec_template/ex_gal_sf.dat')
## ETC package ##
    print "initializing ETC package ..."
    etc = Etc()
#    etc.get_sn_table()
## Simulated spectra ##
    print "generating simulated spectra ..."
    for i in range(len(ets.assigned_tgt)):
        print "ObjId: %d" % (ets.assigned_tgt[i])
        stmp.get_norm_template(tgt.dict_tgt_z[ets.assigned_tgt[i]],tgt.dict_tgt_mag[ets.assigned_tgt[i]],out_dir+'mag.dat')
        etc.get_simulated_spectra(mag_file=out_dir+'mag.dat',objId=ets.assigned_tgt[i])

#    print pfi.get_num_cobra()
#    for i in range(pfi.get_num_cobra()):
#        print pfi.pfi_x[i],pfi.pfi_y[i],pfi.pfi_x_rot[i],pfi.pfi_y_rot[i]
#    for i in range(len(tgt.tgt_id)):
#        print tgt.tgt_id[i],tgt.tgt_x[i],tgt.tgt_y[i]

if __name__ == '__main__':
    main()
