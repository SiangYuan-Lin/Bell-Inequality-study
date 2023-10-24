import numpy as np # linear algebra

def vec_unit(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    return x/r, y/r, z/r
vvec_unit = np.vectorize(vec_unit)

def vec_mag(x, y, z):
    return np.sqrt(x**2 + y**2 + z**2)
vvec_mag = np.vectorize(vec_mag)

def vec_mag_2(x, y, z):
    return (x**2 + y**2 + z**2)
vvec_mag_2 = np.vectorize(vec_mag)

def vec_cos(ax, ay, az, bx, by, bz):
    ra = np.sqrt(ax**2 + ay**2 + az**2)
    rb = np.sqrt(bx**2 + by**2 + bz**2)
    return (ax*bx + ay*by + az*bz)/(ra*rb)
vvec_cos = np.vectorize(vec_cos)

def inv_mass(E, px, py, pz):
    P2 = px**2 + py**2 + pz**2
    return np.sqrt(E**2 - P2)
vinv_mass = np.vectorize(inv_mass)

def inv_mass_2(E, px, py, pz):
    P2 = px**2 + py**2 + pz**2
    return (E**2 - P2)
vinv_mass_2 = np.vectorize(inv_mass_2)

class Vector:
    def __init__(self,X,Y,Z):
        R = vec_mag(X,Y,Z)
        self.x = X
        self.y = Y
        self.z = Z
        self.r = R
        self.r2= R**2
    def __rmul__(self, A):
        return Vector(A*self.x,A*self.y,A*self.z)
    def __mul__(self, A):
        return Vector(A*self.x,A*self.y,A*self.z)
    def __truediv__(self, A):
        return Vector(self.x/A,self.y/A,self.z/A)
    def __add__(self, A):
        return Vector(A.X() + self.x, A.Y() + self.y, A.Z() + self.z)
    def __radd__(self, A):
        return Vector(A.X() + self.x, A.Y() + self.y, A.Z() + self.z)
    def __sub__(self, A):
        return Vector(self.x - A.X(), self.y - A.Y(), self.z - A.Z())
    def __rsub__(self, A):
        return Vector(A.X() - self.x, A.Y() - self.y, A.Z() - self.z)
    def X(self):
        return self.x
    def Y(self):
        return self.y
    def Z(self):
        return self.z
    def R(self):
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)
    def R2(self):
        return (self.x**2 + self.y**2 + self.z**2)
    def unit(self):
        return Vector(self.x/self.r,self.y/self.r,self.z/self.r)
    def dot(self, V):
        return self.x*V.x + self.y*V.y + self.z*V.z
    def cos(self,V):
        return self.dot(V)/(self.r * V.r)
    def theta(self, V):
        return np.arccos(self.dot(V))
    def cross(self,V):
        vi =  self.y*V.z - self.z*V.y #123 + 132
        vj =  self.z*V.x - self.x*V.z #231 + 213 
        vk =  self.x*V.y - self.y*V.x #312 + 321
        return Vector(vi,vj,vk)
    def components(self):
        print(self.x, self.y, self.z)
        
class LorentzVector:
    p0 = 0
    p1 = 0
    p2 = 0
    p3 = 0
    #Four momentum elements are (px,py,pz,E)
    def __init__(self, p1,p2,p3,p0):#, mass): 
        self.e  = p0
        self.px = p1
        self.py = p2
        self.pz = p3
    def Px(self):
        return self.px
    def Py(self):
        return self.py
    def Pz(self):
        return self.pz
    def E(self):
        return self.e
    def SetPx(self,p1):
        return LorentzVector(     p1,self.py,self.pz,self.e)#, self.m)
    def SetPy(self,p2):
        return LorentzVector(self.px,     p2,self.pz,self.e)# ,self.m)
    def SetPz(self,p3):
        return LorentzVector(self.px,self.py,     p3,self.e)# ,self.m)
    def SetE( self,p0):
        return LorentzVector(self.px,self.py,self.pz,     p0)#,self.m)
    def SetP4(p1,p2,p3,p0):#, mass):
        return LorentzVector(p1,p2,p3,p0)#,mass)
    def P( self): 
        return np.sqrt(self.px**2 + self.py**2 + self.pz**2)
    def P2(self):
        return self.px**2 + self.py**2 + self.pz**2
    def P3(self):
        return Vector(self.px, self.py, self.pz)
    def M2(self):
        return (self.e**2 - self.P2())
    def Gamma(self):
        if(self.e == 0):
            if (self.P2() == 0): return 0 # to avoid Nan
            else: return 0
        elif (self.P2() >  self.e**2): return 0 # spacelike
        elif (self.P2() == self.e**2): return np.Inf # lightlike
        else: return 1/np.sqrt(1 - self.P2()/(self.e**2))
    def Beta(self):
        if(self.e == 0):
            if (self.P2() == 0): return 1
            else: return 0
        elif (self.M2() <= 0): return np.Nan
        else: return self.P()/self.e
    def BoostVector(self):
        return Vector(self.px/self.e,self.py/self.e,self.pz/self.e)
    def Boost(self,V):
        vp = V.dot(self.P3())
        gamma = 1/np.sqrt(1 - V.R2())
        gamma2 = (gamma - 1)/V.R2() if V.R2() > 0 else 0
        p1 = self.px + gamma2*vp*V.X() + gamma*V.X()*self.e
        p2 = self.py + gamma2*vp*V.Y() + gamma*V.Y()*self.e
        p3 = self.pz + gamma2*vp*V.Z() + gamma*V.Z()*self.e
        p0 = gamma*(vp + self.e)
        return LorentzVector(p1,p2,p3,p0) 
    def components(self):
        print(self.px,self.py,self.pz,self.e)
    def __add__(self, A):
        return LorentzVector(A.Px() + self.px, A.Py() + self.py, A.Pz() + self.pz, A.E() + self.e)
    def __radd__(self, A):
        return LorentzVector(self.px + A.Px(), self.py + A.Py(), self.pz + A.Pz(), self.e + A.E())
    def __sub__(self, A):
        return LorentzVector(A.Px() - self.px, A.Py() - self.py, A.Pz() - self.pz, A.E() - self.e)
    def __rsub__(self, A):
        return LorentzVector(self.px - A.Px(), self.py - A.Py(), self.pz - A.Pz(), self.e - A.E())

def four_momentum(particle, p4, idx):
    if p4 == None: p4=['E','px','py','pz']
    return LorentzVector(particle[p4[1]][idx],particle[p4[2]][idx],particle[p4[3]][idx],particle[p4[0]][idx])

def BoostToRest(Part, part_label, PartP4, Rest, frame_label, RestP4):
    # four-vector component label
    PartP4 = ['E','px','py','pz'] if PartP4 == None else PartP4 
    RestP4 = ['E','px','py','pz'] if RestP4 == None else RestP4
    Rest['v2'] = (Rest[RestP4[1]]/Rest[RestP4[0]])**2\
               + (Rest[RestP4[2]]/Rest[RestP4[0]])**2\
               + (Rest[RestP4[3]]/Rest[RestP4[0]])**2
    Rest['gamma'] = 1/np.sqrt(1 -  Rest['v2'])
    Rest['gamma2'] = (Rest['gamma'] - 1)/Rest['v2']
    # "part_label" names the particle to be boosted
    # e.g. 'lp' for lepton + ['vp_lp','g2_lp']
    vp = 'vp_' + part_label # inner product of boost vector and three-momentum
    Rest[vp] = Part[PartP4[1]]*Rest[RestP4[1]]/Rest[RestP4[0]]\
             + Part[PartP4[2]]*Rest[RestP4[2]]/Rest[RestP4[0]]\
             + Part[PartP4[3]]*Rest[RestP4[3]]/Rest[RestP4[0]]
    # "frame_label" names the rest frame to which "Part" is boosted to
    # e.g. 'H' for Higgs rest frame ['E_H','px_H','py_H','pz_H']
    Part[PartP4[1] + '_' +frame_label] = Part[PartP4[1]] + Rest['gamma2']*Rest[vp]*Rest[RestP4[1]]/Rest[RestP4[0]] - Rest['gamma']*Rest[RestP4[1]]*Part[PartP4[0]]/Rest[RestP4[0]]
    Part[PartP4[2] + '_' +frame_label] = Part[PartP4[2]] + Rest['gamma2']*Rest[vp]*Rest[RestP4[2]]/Rest[RestP4[0]] - Rest['gamma']*Rest[RestP4[2]]*Part[PartP4[0]]/Rest[RestP4[0]]
    Part[PartP4[3] + '_' +frame_label] = Part[PartP4[3]] + Rest['gamma2']*Rest[vp]*Rest[RestP4[3]]/Rest[RestP4[0]] - Rest['gamma']*Rest[RestP4[3]]*Part[PartP4[0]]/Rest[RestP4[0]]
    Part[PartP4[0] + '_' +frame_label] = Rest['gamma']*(Part[PartP4[0]] - Rest[vp])