import pyasf
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

class Sample(object):
    def __init__(self, substrate, *Layers):
        self.substrate=substrate
        self.Layers=Layers
        
    def set_Miller(self, Miller):
        self.substrate.set_Miller(Miller)
        self.substrate.calc_H_cartesian()
    
    def calc_g0_gH(self, Energy):
        self.substrate.calc_g0_gH(Energy)
        for layer in self.Layers:
            layer.g0=self.substrate.g0
            layer.gH=self.substrate.gH
    
    def add_Layer(self, Layer):
        self.Layers.append(Layer)
        
    def calc_layer_Miller(self):
        for layer in self.Layers:
            M=layer.structure.M.subs(layer.structure.subs).evalf()
            R=layer.R.subs(layer.structure.subs).evalf()
            RS=self.substrate.R.subs(self.substrate.structure.subs).evalf()
            H=RS*self.substrate.H.subs(self.substrate.structure.subs).evalf()
            Recq=M.T*R.T*H
            layer.Miller=sp.Matrix([int(round(element)) for element in Recq])
            u=H.normalized()
            k_in_unit = self.substrate.k_in_unit(sp.Symbol("theta"))
            thetalayer= pyasf.makefunc(-sp.asin(k_in_unit.dot(u)))
            layer.thetalayer=thetalayer

    def calc_reflectivity(self, theta, Energy, Polarization=1):
        self.calc_g0_gH(Energy)
        X0=self.substrate.calc_reflection_amplitude(theta, Energy)
        for layer in self.Layers:
            XR=layer.calc_reflection_amplitude(layer.thetalayer(theta), Energy)
            XT=layer.calc_transmission_amplitude(layer.thetalayer(theta), Energy)
            Xt=(XR-X0*(XR**2-XT**2))/(1-X0*XR)
            X0=Xt
        return X0
        
    def print_values(self, theta, Energy):
        print "wavelenght (A): ",12398./Energy, " Energy (eV): ", Energy
        print "Bragg reflection: "
        for layer in self.Layers:
            print layer.Miller
        print "Surface: ", self.substrate.v_perp
        phi=np.degrees(float(sp.acos(self.substrate.H.normalized().dot(self.substrate.w3)).subs(self.substrate.structure.subs)))
        print "Bragg plane angle to surface (degrees): ", phi
        thetaBragg=self.substrate.calc_Bragg_angle(Energy).subs(self.substrate.structure.subs).evalf()
        thetaBraggdegree=np.degrees(float(thetaBragg))
        print "Bragg angle (degrees): ", thetaBraggdegree
        inc_angle=np.degrees(-np.arcsin(float(self.substrate.k_in_unit(thetaBragg)[2])))
        print "Incident angle (degrees): ", inc_angle
        ex_angle=np.degrees(np.arcsin(float(self.substrate.k_sc_unit(thetaBragg)[2])))
        print "Exit angle (degrees): ", ex_angle
        b=self.substrate.g0/self.substrate.gH
        print "Asymmetry factor: ", b
        if self.substrate.w1.dot(self.substrate.H)==0:
            print "Symmetric"
        if (self.substrate.w1.cross(self.substrate.w3)).dot(self.substrate.H.normalized())==0:
            print "Coplanar"
            
class Epitaxial_Layer(object):
    def __init__(self, structure, thickness, Miller=None, R=None, v_par=None, v_perp=None):
        self.structure=structure
        self.thickness=thickness
        self.R = R
        self.Miller=Miller
        self.v_par=v_par
        self.v_perp=v_perp
    
    def calc_orientation(self, v_par=None, v_perp=None):
        """
            Define the Orientation Matrix of the Layer with respect to the
            sample system.
            
            Inputs:
                v_par : sequence of length 3
                    vector of the reciprocal lattice system which is pointing parallel
                    to the sample surface and parallel to the scattering plane
                v_perp : sequence of length 3
                    vector of the reciprocal lattice system which is pointing perpendicular
                    to the sample surface
        """
        if v_par==None:
            v_par=self.v_par
        if v_perp==None:
            v_perp=self.v_perp
        self.v_par=v_par
        self.v_perp=v_perp
        M=self.structure.M
        w1=M.T.inv()*v_par/((M.T.inv()*v_par).norm())
        w3=M.T.inv()*v_perp/((M.T.inv()*v_perp).norm())
        if w1.dot(w3)!=0:
            print("v_par and v_perp are not perpendicular")
        w2=w3.cross(w1)
        self.w1=w1
        self.w3=w3
        R=sp.Matrix([w1.T, w2.T, w3.T])
        self.R=R
        return R
    
    def calc_orientation_from_angle(self, psi, v_perp=None):
        if v_perp==None:
            v_perp=self.v_perp
        self.v_perp=v_perp
        M=self.structure.M
        w3=M.T.inv()*v_perp/((M.T.inv()*v_perp).norm())
        self.w3=w3
        b_rec=sp.Matrix([0,1,0])
        b=M.T.inv()*b_rec
        pv=b.cross(w3)
        if pv.norm()==0:
            c_rec=sp.Matrix([0,0,1])
            c=M.T.inv()*c_rec
            pv=w3.cross(c)
        p=pv.normalized()
        q=p.cross(w3)
        w1=p*sp.cos(psi)+q*sp.sin(psi)
        self.w1=w1
        self.v_par=M*w1
        w2=w3.cross(w1)
        R=sp.Matrix([w1.T, w2.T, w3.T])
        self.R=R
        return R
            
    def calc_parameters(self, theta, Energy):
        Miller=self.Miller
        struct=self.structure
        wavelenght=12398./Energy
        Volume=struct.V.subs(struct.subs).evalf()
        r_e=2.818e-5
        Gamma=r_e*wavelenght**2/(sp.pi*Volume)
        g0=self.g0
        gH=self.gH
        b=g0/gH
        C=1
        
        thetaBragg=self.calc_Bragg_angle(Energy).subs(struct.subs).evalf()
        FH=struct.DAFS(Energy, Miller)
        FHc=struct.DAFS(Energy, tuple([-i for i in Miller]))
        F0=struct.DAFS(Energy, (0,0,0))

        
        thetasym = sp.Symbol("theta", real=True)
        eta=(-b*(thetasym-thetaBragg)*sp.sin(2*thetaBragg)-Gamma*F0[0]*(1-b)/2)/(sp.sqrt(sp.Abs(b))*C*Gamma*sp.sqrt(FH[0]*FHc[0]))
        etafunc = pyasf.makefunc(eta)
        self.etafunc = etafunc
        etaval = etafunc(theta)
        self.etaval=etaval
        
        
        T=sp.pi*C*Gamma*sp.sqrt(FH[0]*FHc[0])*self.thickness/(wavelenght*sp.sqrt(abs(g0*gH)))

        self.T=T = complex(T.evalf(subs=struct.subs))
        
        alpha=T*np.sqrt(etaval**2-1)
        Q=np.sqrt(etaval**2-1)*np.cos(alpha)+1j*etaval*np.sin(alpha)
        self.alpha=alpha
        self.Q=Q
        
    def calc_reflection_amplitude(self, theta, Energy):
        self.calc_parameters(theta, Energy)
        alpha=self.alpha
        Q=self.Q
        XR=1j*np.sin(alpha)/Q
        return XR
    
    def calc_transmission_amplitude(self, theta, Energy):
        self.calc_parameters(theta, Energy)
        alpha=self.alpha
        Q=self.Q
        etaval=self.etaval
        XT=np.sqrt(etaval**2-1)/Q
        return XT
    
    def calc_H_cartesian(self):
        Miller=sp.Matrix(self.Miller)
        M=self.structure.M
        H=M.T.inv()*Miller
        self.H=H
        return H
        
    def calc_Bragg_angle(self, Energy, Miller=None):
        if Miller == None:
            Miller = self.Miller
        wavelength=12398./Energy
        self.structure.hkl(*Miller)
        Braggq= 2*sp.pi * self.structure.qfunc.dictcall(self.structure.subs)
        #thetaBragg=self.structure.theta_degrees()
        thetaBragg=sp.asin(Braggq*wavelength/(4*sp.pi))
        return thetaBragg
        
    def calc_wavevectors(self):
        #thetaBragg=self.calc_Bragg_angle(Energy).subs(struct.subs).evalf()
        theta = sp.Symbol("theta", real=True)
        w1=self.w1
        n=self.w3
        H=self.calc_H_cartesian()
        u=H/H.norm()
        wv=w1-(w1.dot(u))*u
        w=wv/wv.norm()
        k_in_unit=(self.R*(-sp.sin(theta)*u+sp.cos(theta)*w)).subs(self.structure.subs).evalf()
        k_sc_unit=(self.R*( sp.sin(theta)*u+sp.cos(theta)*w)).subs(self.structure.subs).evalf()
        self.k_in_unit=pyasf.makefunc(k_in_unit, "sympy")
        self.k_sc_unit=pyasf.makefunc(k_sc_unit, "sympy")
        self.q_unit=u
        
class Substrate(Epitaxial_Layer):
    def __init__(self, structure):
        super(Substrate, self).__init__(structure, np.inf)
    
    
    def set_Miller(self, Miller):
        self.Miller=Miller
        self.calc_wavevectors()
    
    def calc_g0_gH(self, Energy):
        thetaBragg=self.calc_Bragg_angle(Energy).subs(self.structure.subs).evalf()
        k_in_unit_B=self.k_in_unit(float(thetaBragg))
        k_sc_unit_B=self.k_sc_unit(float(thetaBragg))
        g0=k_in_unit_B[2]
        gH=k_sc_unit_B[2]
        self.g0=g0
        self.gH=gH   

    def calc_reflection_amplitude(self, theta, Energy):
        Miller=self.Miller
        struct=self.structure
        wavelenght=12398./Energy
        Volume=struct.V.subs(struct.subs).evalf()
        r_e=2.818e-5
        Gamma=r_e*wavelenght**2/(sp.pi*Volume)
        g0=self.g0
        gH=self.gH
        b=g0/gH
        C=1
        
        thetaBragg=self.calc_Bragg_angle(Energy).subs(struct.subs).evalf()
        FH=struct.DAFS(Energy, Miller)
        FHc=struct.DAFS(Energy, tuple([-i for i in Miller]))
        F0=struct.DAFS(Energy, (0,0,0))
        
        thetasym = sp.Symbol("theta", real=True)
        eta=(-b*(thetasym-thetaBragg)*sp.sin(2*thetaBragg)-Gamma*F0[0]*(1-b)/2)/(sp.sqrt(sp.Abs(b))*C*Gamma*sp.sqrt(FH[0]*FHc[0]))
        etafunc = pyasf.makefunc(eta)
        self.etafunc = etafunc
        etaval = etafunc(theta)
        s=-np.sign(etaval.real)
        X=etaval+s*np.sqrt(etaval**2-1)
        self.etaval=etaval
        return X
        