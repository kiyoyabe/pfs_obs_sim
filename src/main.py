import math,os,sys
import scipy as sp
import scipy.interpolate as it
import subprocess
import psycopg2

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

class pfsDB():
    def __init__(self):
        return None
    def connect(self,name,passwd):
        self.name = name
        self.passwd = passwd
        self.conn = psycopg2.connect(
        host = "192.168.156.76",
        port = 5432,
        database=self.name,
        user="pfsdbadmin",
        password=self.passwd)
        print('connection to %s ... ok' % (self.name))
        self.cur = self.conn.cursor()
        return None

    def getTargetList(self,filename):
## load target data ##
        print "loading target list ..."
        self.cur.execute("SELECT \"targetId\",\"ra\",\"dec\",\"fiducialExptime\",\"priority\",\"fiberMag_g\",\"programId\" FROM \"Target\" WHERE (\"priority\" BETWEEN 1 AND 2) AND (\"fiberMag_g\" < 23.0);")
        self.dat = sp.array(self.cur.fetchall(),dtype=[('targetId','i4'),('ra','f8'),('dec','f8'),('fiducialExptime','i4'),('priority','i4'),('fiberMag_g','f8'),('programId','i4')])
        self.targetId = self.dat['targetId']
        self.ra = self.dat['ra']
        self.dec = self.dat['dec']
        self.fiducialExptime = self.dat['fiducialExptime']
        self.priority = self.dat['priority']
        self.fiberMag_g = self.dat['fiberMag_g']
        self.programId = self.dat['programId']
        self.zph = sp.random.uniform(0.5,2.5,len(self.targetId))
#        print self.targetId
## save data ##
        file = open(filename,'w')
        for i in range(len(self.targetId)):
            file.write('%10s   %9.5f  %9.5f  %7.1f  %2d  %6.3f  %3d  N/A\n' % ('ID%d'%(self.targetId[i]),self.ra[i],self.dec[i],self.fiducialExptime[i],self.priority[i],self.fiberMag_g[i],self.programId[i]))
        file.close()
#        self.targetId = dat[0]
#        print self.targetId
        return 0

    def getTargetDic(self):
        self.ra_dic = {}
        self.dec_dic = {}
        self.fiducialExptime_dic = {}
        self.priority_dic = {}
        self.fiberMag_g_dic = {}
        self.programId_dic = {}
        self.zph_dic = {}
        for i in range(len(self.targetId)):
            self.ra_dic[self.targetId[i]] = self.ra[i]
            self.dec_dic[self.targetId[i]] = self.dec[i]
            self.fiducialExptime_dic[self.targetId[i]] = self.fiducialExptime[i]
            self.priority_dic[self.targetId[i]] = self.priority[i]
            self.fiberMag_g_dic[self.targetId[i]] = self.fiberMag_g[i]
            self.programId_dic[self.targetId[i]] = self.programId[i]
            self.zph_dic[self.targetId[i]] = self.zph[i]
#            return self.ra_dic,self.dec_dic,self.fiducialExptime_dic,self.priority_dic,self.fiberMag_g_dic,self.programId_dic
        return 0

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

class pfsEts():
    def __init__(self):
        cmd = 'cd %s; make' % (ets_dir)
        subprocess.call(cmd, shell=True)
        return None
    def run_assigner(self,assigner,pc,pa):
        self.assigner = assigner
        self.pc = pc
        self.pa = pa
        self.fraction = '0.50'
        self.objlist = '%s/ets.list' % (out_dir)
        self.output = '%s/ets.out' % (out_dir)
        cmd = '%sets_demo assigner=%s input=%s fract=%s output=%s ra=%.5f dec=%.5f posang=%.1f' % (ets_dir,self.assigner,self.objlist,self.fraction,self.output,self.pc[0],self.pc[1],self.pa)
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

class pfsEtc():
    def __init__(self):
        cmd = 'cd %s; make' % (etc_dir)
        proc = subprocess.call(cmd, shell=True)
        return None
    def get_sn_table(self):
        cmd = 'cd %s; python run_etc.py @run_etc.defaults --OUTFILE_NOISE=out/ref.noise.dat --OUTFILE_SNC=out/ref.snc.dat --OUTFILE_SNL=- --EXP_TIME=450 --EXP_NUM=2 --MAG_FILE=23.0' % (etc_dir)
        proc = subprocess.call(cmd, shell=True)
        if proc == 0:
            print "SN table is created!"
        return 0
    def get_simulated_spectra(self,mag_file,objId):
        cmd = 'cd %s; python gen_sim_spec.py @gen_sim_spec.defaults --EXP_NUM=2 --MAG_FILE=%s --etcFile=out/ref.snc.dat --nrealize=1 --outDir=%s/fits/ --asciiTable=None --writeFits=True --objId=%d' % (etc_dir,mag_file,out_dir,objId)
        proc = subprocess.call(cmd, shell=True)
        return proc

def main():
    if os.path.exists(out_dir)==False:
        os.mkdir(out_dir)
######################
## connect to pfsDB ##
######################
    print "Enter DB password:"
    dbpasswd = raw_input()
    db = pfsDB()
    db.connect(name='pfsdbsim',passwd=dbpasswd)
    db.getTargetList(filename='%s/ets.list' % (out_dir))
    db.getTargetDic()
#################
## ETS package ##
#################
    tile = 11
    pc = (33.84251,-4.58684)
    pa = 146.9
    print "initializing ETS package ..."
    ets = pfsEts()
    ets.run_assigner(assigner='naive',pc=pc,pa=pa)
    print "%d objects are assigned!" % (len(ets.assigned_tgt))
#    file = open('test.dat','w')
#    for i in range(len(ets.assigned_tgt)):
#        file.write('%d %.5f %.5f %.5f %.5f\n' % (ets.assigned_tgt[i],ets.assigned_tgt_ra[i],ets.assigned_tgt_dec[i],tgt.dict_tgt_x[ets.assigned_tgt[i]],tgt.dict_tgt_y[ets.assigned_tgt[i]]))
#    file.close()
#######################
## spectral template ##
#######################
    stmp = Template(data_dir+'/spec_template/ex_gal_sf.dat')
#################
## ETC package ##
#################
    print "initializing ETC package ..."
    etc = pfsEtc()
#    etc.get_sn_table()
## Simulated spectra ##
    print "generating simulated spectra ..."
    for i in range(len(ets.assigned_tgt)):
#        print "ObjId: %d" % (ets.assigned_tgt[i])
        stmp.get_norm_template(db.zph_dic[ets.assigned_tgt[i]],db.fiberMag_g_dic[ets.assigned_tgt[i]],out_dir+'mag.dat')
#        etc.get_simulated_spectra(mag_file=out_dir+'mag.dat',objId=ets.assigned_tgt[i])

if __name__ == '__main__':
    main()
