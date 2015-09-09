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
        for layer in self.Layers:
            M=layer.structure.M
            R=layer.R
            RM = (R*M).subs(layer.structure.subs).evalf()
            
            RS=self.substrate.R
            H=RS*self.substrate.H
            Recq=RM.T * H.subs(self.substrate.structure.subs).evalf()
            layer.Miller=sp.Matrix([int(round(element)) for element in Recq])
    
    def calc_g0_gH(self, Energy):
        self.substrate.calc_g0_gH(Energy)
        for layer in self.Layers:
            layer.g0=self.substrate.g0
            layer.gH=self.substrate.gH
    
        
    def calc_theta_layer(self):
        theta=sp.Symbol("theta")
        self.substrate.theta_layer_expr=theta
        for layer in self.Layers:
            M=layer.structure.M
            R=layer.R
            RM = (R*M)
            u=(RM.T.inv()*layer.Miller).normalized()
            layer.u=u
            k_in_unit = self.substrate.k_in_unit(theta)
            theta_layer_expr=-sp.asin(k_in_unit.dot(u))
            layer.theta_layer_expr=theta_layer_expr
            layer.theta_layer_func=pyasf.makefunc(theta_layer_expr, "numpy")

    def calc_reflectivity(self, theta, Energy, Polarization=1):
        self.calc_theta_layer()
        self.calc_g0_gH(Energy)
        self.substrate.calc_amplitudes(theta, Energy)
        X0=self.substrate.XR
        for layer in self.Layers:
            layer.calc_amplitudes(theta, Energy)
            XR=layer.XR
            XT=layer.XT
            Xt=(XR-X0*(XR**2-XT**2))/(1-X0*XR)
            X0=Xt
        return X0
        
    def print_values(self, theta, Energy):
        print "wavelength (A): ",12398./Energy, " Energy (eV): ", Energy
        print "Bragg reflection: "
        for layer in self.Layers:
            print layer.Miller
        print "Surface: ", self.substrate.v_perp
        phi=np.degrees(float(sp.acos(self.substrate.H.normalized().dot(self.substrate.w3)).subs(self.substrate.structure.subs)))
        print "Bragg plane angle to surface (degrees): ", phi
        thetaBragg=self.substrate.thetaBragg
        thetaBraggdegree=np.degrees(thetaBragg)
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
        self.S = self.structure.S
    

        
    
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
    
    
    def calc_parameters(self, Energy):
        Miller=self.Miller
        struct=self.structure
        wavelength=12398./Energy
        Volume=struct.V
        r_e=2.818e-5
        Gamma=r_e*wavelength**2/(sp.pi*Volume)
        g0=self.g0
        gH=self.gH
        b=g0/gH
        C=1
        
        thetaBragg=self.calc_Bragg_angle(Energy)
        FH=struct.DAFS(Energy, Miller)
        FHc=struct.DAFS(Energy, tuple([-i for i in Miller]))
        F0=struct.DAFS(Energy, (0,0,0))

        
        thetasym = sp.Symbol("theta", real=True)
        thicknesssym=sp.Symbol("thickness", real=True, positive=True)
        eta=(-b*(thetasym-thetaBragg)*sp.sin(2*thetaBragg)-Gamma*F0[0]*(1-b)/2)/(sp.sqrt(sp.Abs(b))*C*Gamma*sp.sqrt(FH[0]*FHc[0]))
        self.eta=eta
        
        T=sp.pi*C*Gamma*sp.sqrt(FH[0]*FHc[0])*thicknesssym/(wavelength*sp.sqrt(abs(g0*gH)))

        self.T=T
        
        alpha=T*sp.sqrt(eta**2-1)
        Q=sp.sqrt(eta**2-1)*sp.cos(alpha)+sp.I*eta*sp.sin(alpha)
        self.alpha=alpha
        self.Q=Q
        
    def calc_amplitudes(self, theta, Energy):
        theta_layer_sub=self.theta_layer_func.dictcall(dict(theta=theta, **self.structure.subs))
        self.calc_parameters(Energy)
        alpha=self.alpha
        eta=self.eta
        Q=self.Q
        XR = pyasf.makefunc(sp.I*sp.sin(alpha)/Q, "numpy")
        XT=pyasf.makefunc(sp.sqrt(eta**2-1)/Q, "numpy")
        thissubs=dict(theta=theta, thickness=self.thickness, **self.structure.subs)
        self.XR=XR.dictcall(thissubs)
        self.XT=XT.dictcall(thissubs)
    
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
        thissubs = dict(zip(self.structure.miller, Miller))
        Braggq= 2*sp.pi * self.structure.qfunc.dictcall(thissubs)
        thetaBragg=sp.asin(Braggq*wavelength/(4*sp.pi))
        return thetaBragg
        
    def calc_wavevectors(self):
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
        self.XT=0
    
    
    def set_Miller(self, Miller):
        self.Miller=Miller
        self.calc_wavevectors()
    
    def calc_g0_gH(self, Energy):
        thetaBragg=self.calc_Bragg_angle(Energy).subs(self.structure.subs).evalf()
        self.thetaBragg=float(thetaBragg)
        k_in_unit_B=self.k_in_unit(float(thetaBragg))
        k_sc_unit_B=self.k_sc_unit(float(thetaBragg))
        g0=k_in_unit_B[2]
        gH=k_sc_unit_B[2]
        self.g0=g0
        self.gH=gH   

    def calc_amplitudes(self, theta, Energy):
        Miller=self.Miller
        struct=self.structure
        wavelength=12398./Energy
        Volume=struct.V.subs(struct.subs).evalf()
        r_e=2.818e-5
        Gamma=r_e*wavelength**2/(sp.pi*Volume)
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
        self.XR=X

class Strained_Layer(Epitaxial_Layer):
    def __init__(self, structure, thickness_vector, **strainparam):
        for k in strainparam: # translate strings to symbols
            if isinstance(k, str) and hasattr(structure, k):
                sym = getattr(structure, k)
                strainparam[sym] = strainparam.pop(k)
        
        length=len(thickness_vector)-1
        self.length=length
        
        assert all([len(arr) == length for arr in strainparam.values()]), \
                "Thickness vector and all strain parameters should be of the same length."


        self.strain = strainparam
        super(Strained_Layer, self).__init__(structure, thickness_vector)
        
    def calc_amplitudes(self, theta, Energy):
        self.calc_parameters(Energy)
        Q=self.Q
        eta=self.eta
        alpha=self.alpha
        R=self.R
        M=self.structure.M
        Miller=self.Miller
        XR = pyasf.makefunc(sp.I*sp.sin(alpha)/Q, "numpy")
        XT = pyasf.makefunc(sp.sqrt(eta**2-1)/Q, "numpy")
        XR0=0
        XT0=1
        thissubs = self.structure.subs.copy()
        dt = np.diff(self.thickness)
        for i in xrange(self.length):
            for sym in self.strain:
                thissubs[sym] = self.structure.subs[sym] * (self.strain[sym][i] + 1)
            theta_layer_sub=self.theta_layer_func.dictcall(dict(theta=theta, **thissubs))
            thickness=dt[i]
            subs=dict(theta=theta_layer_sub,thickness=thickness, **thissubs)
            XR1=XR.dictcall(subs)
            XT1=XT.dictcall(subs)
            XT0*=XT1/(1-XR0*XR1)
            XR0=(XR1-XR0*(XR1**2-XT1**2))/(1-XR0*XR1)
        self.XR=XR0
        self.XT=XT0
