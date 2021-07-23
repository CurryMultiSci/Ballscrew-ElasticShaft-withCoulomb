import numpy as np

class Tribology():

    def SliceForceRatio(self, msmax):
        i = np.arange(msmax)
        x1 = -1 + 2 * i / msmax
        x2 = -1 + 2 * (i + 1) / msmax
        ForceRatio = 0.25 * (3 * (x2 - x1) - (x2**3 - x1**3))
        return ForceRatio

    def ReducedYoung(self, E0,m0,E1,m1):
        E = 2.0 / ((1.0 - m0 * m0) / E0 + (1.0 - m1 * m1) / E1)
        return E

    def ReducedMass(self, m0, m1):
        M = 1.0 / (1.0 / m0 + 1.0 / m1)
        return M

    def CompositeRoughness(self, s0,s1):
        s = np.sqrt(s0 * s0 + s1 * s1)
        return s
    
    def ForceRatio(self, lamb):
        if (lamb < 0):
            return 0.0;
        
        return 1.0 - np.exp(-1.8 * np.power(lamb, 1.2))

    def Tangent(self, mu0, us, ur, s):

        return mu0 * 2.0 / np.pi * np.arctan(s * us / ur)


    def BrewHamrock(self, Rx,Ry,dx,E):
        sum_rho = 1. / Rx + 1. / Ry
        Kk = 1.5277 + 0.6023 * np.log(Ry / Rx)
        Ek = 1.0003 + 0.5968 * Rx / Ry
        mu = 0.87959 * np.power(Ek, 1.0/3.0) * np.power((Ry / Rx), 0.4240)
        
        k_el = 1 / (1.0339 * np.power(Ry / Rx, 0.636))
        nu = np.power(2 * Ek * k_el / np.pi, 1.0/3.0);
        
        K2pimu = 2 * Kk / np.pi / mu;
        
        k = np.power(np.power(1.125 / E / E * sum_rho, 1.0/3.0) * K2pimu, -1.5);
        a = np.power(6 * Ek * k / E / sum_rho / np.pi / k_el / k_el , 1.0/3.0) * np.sqrt(dx)
        b = k_el * a
        
        return k, a, b
    
    def Tsuji(self,k,zeta,m,v,x):
        c = 2 * zeta * np.sqrt(1.5 * m * k)
        Fn = k * np.power(x, 1.5) + c * v * np.power(x, 0.25)
        return np.max([Fn, 1e-20])
    
    def surface_velocity(self, x0, p0, v0, w0):
        r = p0 - x0
        rw = np.cross(w0, r)
        return rw + v0
    
    def calc_Torque(self, x0, p0, F0):
        r = p0 - x0
        return - np.cross(F0, r)

    def Parabola(self, chi1, chi2, dx):
        c1 = 2 * (chi1 - chi2) * (chi1 + chi2 - chi1 * chi2 * dx)
        c2 = (chi1 * dx - 2) * (chi2 * dx - 2) * (2 * chi1 + 2 * chi2 - chi1 * chi2 * dx)
        return c1 / c2
    
'''

    def Houpert(Rx,a,b,w,E,um,eta,alpha,beta,k):

        p = 1.5 * w / (np.pi * a * b);
        kap = b / a
        beta_ = b * E / (np.pi * p * Rx);
        W = 2 * np.pi * p * b * b / (3 * kap * E * Rx * Rx)
        U = 2 * eta * um / E / Rx
        M = np.power(beta_, 0.25) * kap * W * np.power(U, -0.75)
        T = 2 * a / Rx * ((0.77 * np.power(beta_, -1 / 3) * np.power(kap, 0.12) * np.power(U, -1. / 12.) - 1.8 * np.power(beta_, -0.25)) / (1 + M / 6.6) + 1.8 * np.power(beta_, -0.25)) * np.power(U, 0.75)
        
        t = T * E * Rx * Rx * Rx

        return t



     AiharaR::calc
    (					
	     Rx,		
	     a,		
	     b,		
	     w,		
	     E,		
	     um,		
	     eta,		
	     alpha,	
	     beta,	
	     k		
    ) 
	     U = eta * um / E / Rx
	     G = alpha * E
	     W = w / E / Rx / Rx
	     L = eta * beta * um * um / k

	     M = 176 / (1 + 0.29 * np.power(L, 0.78)) / alpha * np.power(G * U, 0.658) * np.power(W, 0.31) * Rx * Rx * a
	    return M



     Fujiwara::calc
    (					
	     Rx,		
	     a,		
	     b,		
	     w,		
	     E,		
	     um,		
	     eta,		
	     alpha,	
	     beta,	
	     k		
    ) 
	     U = eta * um / E / Rx
	     G = alpha * E
	     W = w / E / Rx / Rx
	     kap = b / a
	     FHD = 44.6 * np.power(U, 0.694) * np.power(G, 0.961) * np.power(W, 0.401) * (1 - 0.962 * exp(-0.818 * kap))
	     Fr = FHD * Rx * Rx / alpha

	     M = Fr * Rx
	    return M



     RollingResistanceNothing::calc(Rx,  a,  b,  w,  E,  um,  eta,  alpha,  beta,  k) 
	    return 0.0;




     HamrockDowson::calc
    (					
	     w,		
	     E,		
	     Rx,		
	     Ry,		
	     alpha,	
	     eta,		
	     u,		
	     lm,		
	     bh		
    ) 
	     W = w / E / Rx / Rx
	     G = alpha * E
	     U = eta * u / E / Rx
	     k_ = 1.0339 * np.power(Ry / Rx, 0.636)

	     H = 2.69 * np.power(U, 0.67) * np.power(G, 0.53) * np.power(W, -0.067) * (1.0 - 0.61 * exp(-0.73 * k_))
	     h = H * Rx

	    h *=  this->Starvation(H, Rx, lm, bh)
	    return h;



     HamrockDowson::Starvation
    (H,		
	     Rx,		
	     lm,		
	     bh) 
	    if (lm <= -1.0)
		    return 1.0;

	    if (lm <= 0.0)
		    return fabs(lm)

	     f  = lm / bh;
	     fc = 3.06 * np.power(H * Rx * Rx / bh / bh, 0.58)
	    if (f > fc)
		    return 1.0;

	    return np.power(f / fc, 0.29)



     FilmThicknessNothing::calc(w,  E,  Rx,  Ry,  alpha,  eta,  u,  lm,  bh) 
	    return 0.0;



     ErtelGrubin
    (eta,		
	     beta,	
	     k,		
	     u) 
	     L = eta * beta * u * u / k
	     phi = 3.94 / (3.94 + np.power(L, 0.62))
	    return phi



     AiharaT::calc
    (					
	     eta0,	
	     p0,		
	     v0,		
	     u0		
    ) 
	    if (eta0 <= 0.0)
		    return 0.0;

	    if (p0 <= 0.0)
		    return 0.0;

	    if (v0 <= 0.0)
		    return 0.0;

	    if (u0 <= 0.0)
		    return 0.0;

	     s    = v0 / u0;
	     eta = Unit::Pas2cP(eta0)
	     u = Unit::ms2mms(u0)
	     pmax = Unit::Pa2kgfmm2(p0)
	     smax = 1.2448e7 * np.power(eta, -1.5) * np.power(u, -0.75) * np.power(pmax, -1.5)

	     mumax
	    if (pmax > 94.404)
		    mumax = 4.078e-2 * np.power(eta, 0.1) * np.power(u, -0.12) * np.power(pmax, 0.18)
	    else
		    mumax = 1.008e-4 * np.power(eta, 0.1) * np.power(u, -0.12) * np.power(pmax, 1.5)

	     mu = mumax * (1.0 - exp(-2.6 * s / smax))

	    return mu





     Tsuji::calc
    (k,	
	     zeta,
	     m,	
	     v,	
	     x) 
	     c = 2 * zeta * sqrt(1.5 * m * k)
	     Fn = k * np.power(x, 1.5) + c * v * np.power(x, 0.25)

	    return std::max(Fn, 1e-20)



     KelvinVoigt::calc
    (
	     k,	
	     zeta,
	     m,	
	     v,	
	     x	
    ) 
	     c = 2 * zeta * sqrt(m * k)
	     Fn = k * x + c * v

	    return std::max(Fn, 1e-20)


'''
