import numpy as np

class Spiral:

    def __init__(self, alp, l, r, eta, R):
        self.alp = alp;
        self.l   = l;
        self.r   = r;
        self.eta = eta;
        self.R   = R;

        self.l_2pi = l / (2 * np.pi);
        self.nd    = np.sqrt(self.r*self.r + self.l_2pi*self.l_2pi);
        self.l_nd  = self.l_2pi / self.nd;
        self.r_nd  = self.r     / self.nd;
        return;

    def to_eta(self, xyz):
        x      = xyz[0];
        y      = xyz[1];
        z      = xyz[2];
        phi    = x / self.nd / self.l_nd;
        beta   = phi + self.alp;
        s_beta = np.sin(beta);
        c_beta = np.cos(beta);
        y0     = c_beta * y + s_beta * z;
        z0     =-s_beta * y + c_beta * z;
        th     = z0 / (y0 + self.l_nd * self.l_2pi / self.r_nd);
        ze     = -th * self.l_2pi / self.r_nd;
        et     = self.r_nd * self.nd - y0 - y0 * th*th / 2 + self.l_nd * th * ze;
        eta  = np.array([th+phi, et, ze]);
        return eta;

    def to_xyz(self, eta):
        th     = eta[0];
        beta   = th + self.alp;
        s_beta = np.sin(beta);
        c_beta = np.cos(beta);
        q    = np.array([self.l_2pi * th, self.r * c_beta, self.r * s_beta]);
        n    = np.array([0.0, -c_beta, -s_beta]);
        b    = np.array([self.r_nd, self.l_nd * s_beta, -self.l_nd * c_beta]);
        xyz  = q + eta[1] * n + eta[2] * b;
        return xyz;

    def to_eta2(self, xyz):
        eta0 = self.to_eta(xyz);
        xyz0 = self.to_xyz(eta0);
        xyz1 = 2 * xyz - xyz0;
        eta1 = self.to_eta(xyz1);
        return eta1;
    
    def get_xyz2eta(self, theta):
        beta   = theta + self.alp;
        s_beta = np.sin(beta);
        c_beta = np.cos(beta);

        xyz2eta = np.array([
            [self.l_nd, -self.r_nd * s_beta, self.r_nd * c_beta,],
            [0.0, -c_beta, -s_beta,],
            [self.r_nd, self.l_nd * s_beta, -self.l_nd * c_beta]
            ])
        return xyz2eta;
    
    def get_surface(self, t, a, i):

        eta = np.zeros(3)
        eta[0] = t
        eta[1] = self.R[i] * np.cos(a) + self.eta[i, 0]
        eta[2] = self.R[i] * np.sin(a) + self.eta[i, 1]

        return self.to_xyz(eta)

    def get_mesh(self, th, alp, i):

        xyz = np.zeros(np.array([th.shape[0], alp.shape[0], 3]))

        for idt, t in enumerate(th):
            for ida, a in enumerate(alp):
                xyz[idt, ida, :] = self.get_surface(t, a, i)

        return xyz
    
    def get_rho(self, cos_alp, i):
        
        R = self.nd * self.nd / self.r - self.eta[i, 0]
        r_inv = 1.0 / self.R[i]
        
        return np.array([cos_alp / (R - self.R[i] * cos_alp), -r_inv])
        
    def get_contact(self, xyz, R, i):
        eta = self.to_eta2(xyz)[1:]
        eF  = eta - self.eta[i, :]
        print(eF)
        return 0


  
        