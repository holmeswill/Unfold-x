# coding: utf-8
import numpy as np
import pylab as plt
import sys
from lxml import etree, objectify
import glob
import unittest


class TestProj(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.prefix = prefix = 'test'
        cls.nkpt   = nkpt   = 159
        cls.nbnd   = nbnd   = 44
        cls.cell_scale      = 4 # new cell has 4 atoms instead of 1, 4 times more bands, etc
        

        bc = np.loadtxt('BC/band.dat.gnu')
        cls.Ebc = bc[:,1].reshape((-1, cls.nkpt))

        E=np.zeros([nbnd,nkpt])
        P=np.zeros([nbnd,nkpt])
        
        for fname in glob.glob('SC/{prefix}.dat.save/*.xml'.format(prefix=prefix)):
            with open(fname,'r') as fobj:
                xml = fobj.read()
            data = objectify.fromstring(xml)
            ik = int(data.EIGENVALUES.attrib['ik'])
            ispin = int(data.EIGENVALUES.attrib['ispin'])
            E[:,ik-1]=[float(x) for x in data.EIGENVALUES.text.split()]
            P[:,ik-1]=[float(x) for x in data.PROJECTIONS.text.split()]
        
        
        with open('BC/{}.scf.out'.format(prefix),'r') as f:
            for line in f.readlines():
                if '(alat)' in line:
                    BCalat = float(line.split()[-2])
                if 'Fermi' in line:
                    EFBC = float(line.split()[-2])
                    break
        with open('SC/{}.scf.out'.format(prefix),'r') as f:
            for line in f.readlines():
                if '(alat)' in line:
                    SCalat = float(line.split()[-2])
                if 'Fermi' in line:
                    EFSC = float(line.split()[-2])
                    break
        cls.alatRatio = BCalat/SCalat
        cls.P = P
        cls.E = E


    def test_projection(self):
        nkpt = self.nkpt
        nbnd = self.nbnd
        cell_scale = self.cell_scale
        P=self.P
        E=self.E
        
        F=np.where(P>0.1, E, 0)
        for i in range(0,nbnd,cell_scale):
            S = np.zeros(nkpt)
            PE = np.zeros(nkpt)
            # skip incomplete bands
            for j in range(min(cell_scale,(nbnd-i)//cell_scale)):
                S  += P[j+i,:]
                PE += F[j+i,:]
            else:
                break # this will avoit testing incomplete bands
            np.testing.assert_allclose(S, 1., rtol=1e-5, atol=0)
            np.testing.assert_allclose(PE, Ebc[i,:], rtol=1e-5, atol=0)

    def test_path(self):
        nkpt = self.nkpt
        nbnd = self.nbnd

        data = np.fromfile('SC/{}.dat'.format(self.prefix))
        A=data.reshape([nkpt,-1,3])
        bc = np.loadtxt('BC/band.dat.gnu')
        
        np.testing.assert_allclose( A[:,0,0]*self.alatRatio, bc[:nkpt,0],atol=0.01)


if __name__ == '__main__':
    unittest.main()
