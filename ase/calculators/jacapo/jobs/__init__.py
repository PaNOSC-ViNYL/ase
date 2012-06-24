import os
from subprocess import *
from ase.calculators.jacapo import *
from torque.torque import *

class JobSubmitted(exceptions.Exception):
    pass

class JobInQueue(exceptions.Exception):
    pass

class JobDone(exceptions.Exception):
    pass

class UnknownJobStatus(exceptions.Exception):
    pass

dacapo_run = '''\
#!/usr/bin/env python
import torque
from ase.calculators.jacapo import *
fp = '%(NCFILE)s'
calc = Jacapo(fp)
calc.calculate()
calc.strip()
'''

bfgs_run = '''\
#!/usr/bin/env python
import os, torque
from ase.calculators.jacapo import *
from ase.optimize import BFGS
from ase.io import PickleTrajectory

ncfile = '%(NCFILE)s'
calc = Jacapo(ncfile)
atoms = calc.get_atoms()

opt = BFGS(atoms)
trajfile = ncfile.replace('.nc','.traj')
traj = PickleTrajectory(trajfile,'w',atoms)
opt.attach(traj)
opt.run(fmax=0.05)
#end '''


def queue_job(self,*args,**kwargs):
    '''monkey patch for calculate method to run a calculation script
    through the queue'''

    CWD = os.getcwd()

    wholencfile = self.get_nc()

    if not os.path.exists(wholencfile):
        self.write(wholencfile)

    basepath,NCFILE = os.path.split(wholencfile)
    basename,ext = os.path.splitext(NCFILE)
    TXTFILE = basename + '.txt'
    runningfile =  NCFILE + '.running'
    stopfile = NCFILE + '.stop'

    JOBOUTPUT = basename + '.output'
    JOBIDFILE = basename + '.jobid'

    if basepath != '':
        os.chdir(basepath)
        #now you are in teh directory where the ncfile is

    atoms = Jacapo.read_atoms(NCFILE)
    self = atoms.get_calculator()

    if self.get_status() == 'finished':
        #this means the job is done.
        #do some clean-up of unimportant files
        for jfile in [JOBIDFILE,
                      runningfile,
                      stopfile]:
            if os.path.exists(jfile):
                os.remove(jfile)
        #slave files from parallel runs
        import glob
        slvpattern = TXTFILE + '.slave*'
        for slvf in glob.glob(slvpattern):
            os.unlink(slvf)

        #exit so that we can move on.
        os.chdir(CWD)
        return True

    #Past here means we have to check if the calculation is running
    if os.path.exists(JOBIDFILE):
        JOBID = open(JOBIDFILE).readline()
    else:
        JOBID = None

    # get a queue object
    pbs = PBS()
    pbs.fastpoll()

    if JOBID is not None:

        fields = JOBID.split('.')
        jobnumber = fields[0]

        #the job has been submitted before, and we need
        #to find out what the status of the jobid is
        for job in pbs:

            if job['Job Id'] == jobnumber + '.gilgamesh':
                if job['job_state'] == 'R':
                    os.chdir(CWD)
                    raise JobRunning, job['job_state']
                elif job['job_state'] == 'Q':
                    os.chdir(CWD)
                    raise JobInQueue, job['job_state']
                elif  job['job_state'] == 'C':
                    os.chdir(CWD)
                    raise JobDone, job['job_state']
                else:
                    os.chdir(CWD)
                    raise UnknownJobStatus, job['job_state']
        # if you get here, the job is not in the queue anymore
        # getting here means the calculation was not finished, and is not in
        # the queue anymore. so we print the output. an exception should get raised.
        print open(JOBOUTPUT).readlines()

        if os.path.exists(TXTFILE): #sometimes a job may not even start and there is no txtfile
            #check output of Dacapo
            f = open(TXTFILE,'r')
            for line in f:
                if 'abort_calc' in line:
                    f.close()
                    os.chdir(CWD)
                    raise DacapoAborted(line,self)
                continue
            f.close()

            #check last line for proper finish
            if not ('clexit: exiting the program' in line
                    or 'PAR: msexit halting Master' in line):
                os.unlink(JOBIDFILE)
                os.chdir(CWD)
                print line
                raise DacapoNotFinished('Job is over, but Dacapo is not done',self)

            # something else must be wrong
            os.unlink(JOBIDFILE)
            print self.nc
            os.chdir(CWD)
            raise Exception,'something is wrong with your job!'
        else:
            os.chdir(CWD)
            print '%s does not exist!' % TXTFILE
            raise Exception,'%s does not exist!' % TXTFILE

    #Past here, we need to submit a job.
    p = Popen(['qsub',
               '-joe',
               '-N',
               "%s" % wholencfile,
               '-o',
               "%s" % JOBOUTPUT,
               '-l walltime=168:00:00'],
              stdin=PIPE,stdout=PIPE,stderr=PIPE)

    if not hasattr(self,'script'):
        #by default use dacapo_run
        script = dacapo_run % locals()
    else:
        script = eval(self.script) % locals()

    out,err = p.communicate(script)
    fields = out.split('.')
    f = open(JOBIDFILE,'w')
    f.write('.'.join(fields[0:2]))
    f.close()

    if len(err) > 0:
        print 'err = "%s"' % err

    os.chdir(CWD)
    raise JobSubmitted, out


if __name__ == '__main__':

    print eval('dacapo_run')
