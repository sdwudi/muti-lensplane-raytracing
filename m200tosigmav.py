import numpy as np
from astropy.cosmology import Planck15 as cosmo
"set constant value"
Msun=1.98847*pow(10,30)#kg,cite wikipedia solar mass
c_v=2.99792458*pow(10,8)#m/s
G=6.6740831*pow(10,-11)#m^3/(kg*s^2)
Mpc = 30.85734528 * pow(10, 21)#unit:m
CH=1000/(30.85734528*pow(10,21))#unit:km/Mpc
tza=cosmo.angular_diameter_distance
h=np.array(cosmo.H(0))/100

"calculate sigmav using mass_200 and redshift,cite:mo,mao,white(1998)"
def cal(m,z):#sigma_v
    sigma_v=np.power(m*Msun*10*G*np.array(cosmo.H(z))*CH,1/3)/np.sqrt(2)#cite:mo,mao,white(1998)
    return sigma_v/1000    #单位km/s

#velocity dispersion= rotation velocity/sqrt(2)--(stt012)Chae(2010)


#r_200--rhi_crit
"calculate theta_200 using mass_200 and redshift,cite:mo,mao,white(1998)"
def M_2002thate_200(m,z):


    # rho_crit = 3 * np.array(cosmo.H(z) * CH)**2/8/np.pi/G
    #
    # r_3=m*Msun/800*3/rho_crit/np.pi
    #
    # r_200=np.power(r_3,1/3)

    r_200=cal(m,z)*np.sqrt(2)/(10*np.array(cosmo.H(z)))#cite:Mo,Mao,White(1998)
    theta_200=r_200/np.array(tza(z))*57.3*60*60
    return theta_200#unit:arcsec

def M_2002R_200(m,z):


    r_200=cal(m,z)*np.sqrt(2)/(10*np.array(cosmo.H(z)))#cite:Mo,Mao,White(1998)

    return r_200#unit:arcsec

"calculate m_star in Er et al.(2012).eq(12)"

def m_star(r,z):#c-m of nfw,from c.source file M_starcal.c(Planck15)
    omega_m = 0.3089
    rho_crit =3 * np.array(cosmo.H(z) * CH) ** 2 / 8 / np.pi / G
    rho_background=rho_crit*omega_m
    v=r**3*4/3*np.pi
    m=v*rho_background
    return m/Msun






# print(cal(10**12,0))

