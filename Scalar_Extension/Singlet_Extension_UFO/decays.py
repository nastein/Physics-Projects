# This file was automatically created by FeynRules 2.3.32
# Mathematica version: 12.0.0 for Mac OS X x86 (64-bit) (April 7, 2019)
# Date: Wed 4 Mar 2020 09:40:59


from object_library import all_decays, Decay
import particles as P


Decay_H = Decay(name = 'Decay_H',
                particle = P.H,
                partial_widths = {(P.h,P.h):'(cmath.sqrt(-4*Mh**2*MH**2 + MH**4)*(lam3**2*svev**2*cmath.cos(alpha)**6 + 12*lam*lam3*svev*vev*cmath.cos(alpha)**5*cmath.sin(alpha) - 4*lam3**2*svev*vev*cmath.cos(alpha)**5*cmath.sin(alpha) + 12*lam2*lam3*svev**2*cmath.cos(alpha)**4*cmath.sin(alpha)**2 - 4*lam3**2*svev**2*cmath.cos(alpha)**4*cmath.sin(alpha)**2 + 36*lam**2*vev**2*cmath.cos(alpha)**4*cmath.sin(alpha)**2 - 24*lam*lam3*vev**2*cmath.cos(alpha)**4*cmath.sin(alpha)**2 + 4*lam3**2*vev**2*cmath.cos(alpha)**4*cmath.sin(alpha)**2 + 72*lam*lam2*svev*vev*cmath.cos(alpha)**3*cmath.sin(alpha)**3 - 24*lam*lam3*svev*vev*cmath.cos(alpha)**3*cmath.sin(alpha)**3 - 24*lam2*lam3*svev*vev*cmath.cos(alpha)**3*cmath.sin(alpha)**3 + 10*lam3**2*svev*vev*cmath.cos(alpha)**3*cmath.sin(alpha)**3 + 36*lam2**2*svev**2*cmath.cos(alpha)**2*cmath.sin(alpha)**4 - 24*lam2*lam3*svev**2*cmath.cos(alpha)**2*cmath.sin(alpha)**4 + 4*lam3**2*svev**2*cmath.cos(alpha)**2*cmath.sin(alpha)**4 + 12*lam*lam3*vev**2*cmath.cos(alpha)**2*cmath.sin(alpha)**4 - 4*lam3**2*vev**2*cmath.cos(alpha)**2*cmath.sin(alpha)**4 + 12*lam2*lam3*svev*vev*cmath.cos(alpha)*cmath.sin(alpha)**5 - 4*lam3**2*svev*vev*cmath.cos(alpha)*cmath.sin(alpha)**5 + lam3**2*vev**2*cmath.sin(alpha)**6))/(32.*cmath.pi*abs(MH)**3)',
                                  (P.b,P.b__tilde__):'(cmath.sqrt(-4*MB**2*MH**2 + MH**4)*(-12*MB**2*yb**2*cmath.sin(alpha)**2 + 3*MH**2*yb**2*cmath.sin(alpha)**2))/(16.*cmath.pi*abs(MH)**3)',
                                  (P.ta__minus__,P.ta__plus__):'(cmath.sqrt(MH**4 - 4*MH**2*MTA**2)*(MH**2*ytau**2*cmath.sin(alpha)**2 - 4*MTA**2*ytau**2*cmath.sin(alpha)**2))/(16.*cmath.pi*abs(MH)**3)',
                                  (P.t,P.t__tilde__):'(cmath.sqrt(MH**4 - 4*MH**2*MT**2)*(3*MH**2*yt**2*cmath.sin(alpha)**2 - 12*MT**2*yt**2*cmath.sin(alpha)**2))/(16.*cmath.pi*abs(MH)**3)',
                                  (P.W__minus__,P.W__plus__):'(cmath.sqrt(MH**4 - 4*MH**2*MW**2)*((3*ee**4*vev**2*cmath.sin(alpha)**2)/(4.*sw**4) + (ee**4*MH**4*vev**2*cmath.sin(alpha)**2)/(16.*MW**4*sw**4) - (ee**4*MH**2*vev**2*cmath.sin(alpha)**2)/(4.*MW**2*sw**4)))/(16.*cmath.pi*abs(MH)**3)',
                                  (P.Z,P.Z):'(cmath.sqrt(MH**4 - 4*MH**2*MZ**2)*((9*ee**4*vev**2*cmath.sin(alpha)**2)/2. + (3*ee**4*MH**4*vev**2*cmath.sin(alpha)**2)/(8.*MZ**4) - (3*ee**4*MH**2*vev**2*cmath.sin(alpha)**2)/(2.*MZ**2) + (3*cw**4*ee**4*vev**2*cmath.sin(alpha)**2)/(4.*sw**4) + (cw**4*ee**4*MH**4*vev**2*cmath.sin(alpha)**2)/(16.*MZ**4*sw**4) - (cw**4*ee**4*MH**2*vev**2*cmath.sin(alpha)**2)/(4.*MZ**2*sw**4) + (3*cw**2*ee**4*vev**2*cmath.sin(alpha)**2)/sw**2 + (cw**2*ee**4*MH**4*vev**2*cmath.sin(alpha)**2)/(4.*MZ**4*sw**2) - (cw**2*ee**4*MH**2*vev**2*cmath.sin(alpha)**2)/(MZ**2*sw**2) + (3*ee**4*sw**2*vev**2*cmath.sin(alpha)**2)/cw**2 + (ee**4*MH**4*sw**2*vev**2*cmath.sin(alpha)**2)/(4.*cw**2*MZ**4) - (ee**4*MH**2*sw**2*vev**2*cmath.sin(alpha)**2)/(cw**2*MZ**2) + (3*ee**4*sw**4*vev**2*cmath.sin(alpha)**2)/(4.*cw**4) + (ee**4*MH**4*sw**4*vev**2*cmath.sin(alpha)**2)/(16.*cw**4*MZ**4) - (ee**4*MH**2*sw**4*vev**2*cmath.sin(alpha)**2)/(4.*cw**4*MZ**2)))/(32.*cmath.pi*abs(MH)**3)'})

Decay_h = Decay(name = 'Decay_h',
                particle = P.h,
                partial_widths = {(P.H,P.H):'(cmath.sqrt(Mh**4 - 4*Mh**2*MH**2)*(lam3**2*vev**2*cmath.cos(alpha)**6 - 12*lam2*lam3*svev*vev*cmath.cos(alpha)**5*cmath.sin(alpha) + 4*lam3**2*svev*vev*cmath.cos(alpha)**5*cmath.sin(alpha) + 36*lam2**2*svev**2*cmath.cos(alpha)**4*cmath.sin(alpha)**2 - 24*lam2*lam3*svev**2*cmath.cos(alpha)**4*cmath.sin(alpha)**2 + 4*lam3**2*svev**2*cmath.cos(alpha)**4*cmath.sin(alpha)**2 + 12*lam*lam3*vev**2*cmath.cos(alpha)**4*cmath.sin(alpha)**2 - 4*lam3**2*vev**2*cmath.cos(alpha)**4*cmath.sin(alpha)**2 - 72*lam*lam2*svev*vev*cmath.cos(alpha)**3*cmath.sin(alpha)**3 + 24*lam*lam3*svev*vev*cmath.cos(alpha)**3*cmath.sin(alpha)**3 + 24*lam2*lam3*svev*vev*cmath.cos(alpha)**3*cmath.sin(alpha)**3 - 10*lam3**2*svev*vev*cmath.cos(alpha)**3*cmath.sin(alpha)**3 + 12*lam2*lam3*svev**2*cmath.cos(alpha)**2*cmath.sin(alpha)**4 - 4*lam3**2*svev**2*cmath.cos(alpha)**2*cmath.sin(alpha)**4 + 36*lam**2*vev**2*cmath.cos(alpha)**2*cmath.sin(alpha)**4 - 24*lam*lam3*vev**2*cmath.cos(alpha)**2*cmath.sin(alpha)**4 + 4*lam3**2*vev**2*cmath.cos(alpha)**2*cmath.sin(alpha)**4 - 12*lam*lam3*svev*vev*cmath.cos(alpha)*cmath.sin(alpha)**5 + 4*lam3**2*svev*vev*cmath.cos(alpha)*cmath.sin(alpha)**5 + lam3**2*svev**2*cmath.sin(alpha)**6))/(32.*cmath.pi*abs(Mh)**3)',
                                  (P.b,P.b__tilde__):'((-12*MB**2*yb**2*cmath.cos(alpha)**2 + 3*Mh**2*yb**2*cmath.cos(alpha)**2)*cmath.sqrt(-4*MB**2*Mh**2 + Mh**4))/(16.*cmath.pi*abs(Mh)**3)',
                                  (P.ta__minus__,P.ta__plus__):'((Mh**2*ytau**2*cmath.cos(alpha)**2 - 4*MTA**2*ytau**2*cmath.cos(alpha)**2)*cmath.sqrt(Mh**4 - 4*Mh**2*MTA**2))/(16.*cmath.pi*abs(Mh)**3)',
                                  (P.t,P.t__tilde__):'((3*Mh**2*yt**2*cmath.cos(alpha)**2 - 12*MT**2*yt**2*cmath.cos(alpha)**2)*cmath.sqrt(Mh**4 - 4*Mh**2*MT**2))/(16.*cmath.pi*abs(Mh)**3)',
                                  (P.W__minus__,P.W__plus__):'(((3*ee**4*vev**2*cmath.cos(alpha)**2)/(4.*sw**4) + (ee**4*Mh**4*vev**2*cmath.cos(alpha)**2)/(16.*MW**4*sw**4) - (ee**4*Mh**2*vev**2*cmath.cos(alpha)**2)/(4.*MW**2*sw**4))*cmath.sqrt(Mh**4 - 4*Mh**2*MW**2))/(16.*cmath.pi*abs(Mh)**3)',
                                  (P.Z,P.Z):'(((9*ee**4*vev**2*cmath.cos(alpha)**2)/2. + (3*ee**4*Mh**4*vev**2*cmath.cos(alpha)**2)/(8.*MZ**4) - (3*ee**4*Mh**2*vev**2*cmath.cos(alpha)**2)/(2.*MZ**2) + (3*cw**4*ee**4*vev**2*cmath.cos(alpha)**2)/(4.*sw**4) + (cw**4*ee**4*Mh**4*vev**2*cmath.cos(alpha)**2)/(16.*MZ**4*sw**4) - (cw**4*ee**4*Mh**2*vev**2*cmath.cos(alpha)**2)/(4.*MZ**2*sw**4) + (3*cw**2*ee**4*vev**2*cmath.cos(alpha)**2)/sw**2 + (cw**2*ee**4*Mh**4*vev**2*cmath.cos(alpha)**2)/(4.*MZ**4*sw**2) - (cw**2*ee**4*Mh**2*vev**2*cmath.cos(alpha)**2)/(MZ**2*sw**2) + (3*ee**4*sw**2*vev**2*cmath.cos(alpha)**2)/cw**2 + (ee**4*Mh**4*sw**2*vev**2*cmath.cos(alpha)**2)/(4.*cw**2*MZ**4) - (ee**4*Mh**2*sw**2*vev**2*cmath.cos(alpha)**2)/(cw**2*MZ**2) + (3*ee**4*sw**4*vev**2*cmath.cos(alpha)**2)/(4.*cw**4) + (ee**4*Mh**4*sw**4*vev**2*cmath.cos(alpha)**2)/(16.*cw**4*MZ**4) - (ee**4*Mh**2*sw**4*vev**2*cmath.cos(alpha)**2)/(4.*cw**4*MZ**2))*cmath.sqrt(Mh**4 - 4*Mh**2*MZ**2))/(32.*cmath.pi*abs(Mh)**3)'})

Decay_Z = Decay(name = 'Decay_Z',
                particle = P.Z,
                partial_widths = {(P.W__minus__,P.W__plus__):'(((-12*cw**2*ee**2*MW**2)/sw**2 - (17*cw**2*ee**2*MZ**2)/sw**2 + (4*cw**2*ee**2*MZ**4)/(MW**2*sw**2) + (cw**2*ee**2*MZ**6)/(4.*MW**4*sw**2))*cmath.sqrt(-4*MW**2*MZ**2 + MZ**4))/(48.*cmath.pi*abs(MZ)**3)',
                                  (P.u,P.u__tilde__):'(MZ**2*(-(ee**2*MZ**2) + (3*cw**2*ee**2*MZ**2)/(2.*sw**2) + (17*ee**2*MZ**2*sw**2)/(6.*cw**2)))/(48.*cmath.pi*abs(MZ)**3)',
                                  (P.c,P.c__tilde__):'(MZ**2*(-(ee**2*MZ**2) + (3*cw**2*ee**2*MZ**2)/(2.*sw**2) + (17*ee**2*MZ**2*sw**2)/(6.*cw**2)))/(48.*cmath.pi*abs(MZ)**3)',
                                  (P.t,P.t__tilde__):'((-11*ee**2*MT**2 - ee**2*MZ**2 - (3*cw**2*ee**2*MT**2)/(2.*sw**2) + (3*cw**2*ee**2*MZ**2)/(2.*sw**2) + (7*ee**2*MT**2*sw**2)/(6.*cw**2) + (17*ee**2*MZ**2*sw**2)/(6.*cw**2))*cmath.sqrt(-4*MT**2*MZ**2 + MZ**4))/(48.*cmath.pi*abs(MZ)**3)',
                                  (P.d,P.d__tilde__):'(MZ**2*(ee**2*MZ**2 + (3*cw**2*ee**2*MZ**2)/(2.*sw**2) + (5*ee**2*MZ**2*sw**2)/(6.*cw**2)))/(48.*cmath.pi*abs(MZ)**3)',
                                  (P.s,P.s__tilde__):'(MZ**2*(ee**2*MZ**2 + (3*cw**2*ee**2*MZ**2)/(2.*sw**2) + (5*ee**2*MZ**2*sw**2)/(6.*cw**2)))/(48.*cmath.pi*abs(MZ)**3)',
                                  (P.b,P.b__tilde__):'((-7*ee**2*MB**2 + ee**2*MZ**2 - (3*cw**2*ee**2*MB**2)/(2.*sw**2) + (3*cw**2*ee**2*MZ**2)/(2.*sw**2) - (17*ee**2*MB**2*sw**2)/(6.*cw**2) + (5*ee**2*MZ**2*sw**2)/(6.*cw**2))*cmath.sqrt(-4*MB**2*MZ**2 + MZ**4))/(48.*cmath.pi*abs(MZ)**3)',
                                  (P.ve,P.ve__tilde__):'(MZ**2*(ee**2*MZ**2 + (cw**2*ee**2*MZ**2)/(2.*sw**2) + (ee**2*MZ**2*sw**2)/(2.*cw**2)))/(48.*cmath.pi*abs(MZ)**3)',
                                  (P.vm,P.vm__tilde__):'(MZ**2*(ee**2*MZ**2 + (cw**2*ee**2*MZ**2)/(2.*sw**2) + (ee**2*MZ**2*sw**2)/(2.*cw**2)))/(48.*cmath.pi*abs(MZ)**3)',
                                  (P.vt,P.vt__tilde__):'(MZ**2*(ee**2*MZ**2 + (cw**2*ee**2*MZ**2)/(2.*sw**2) + (ee**2*MZ**2*sw**2)/(2.*cw**2)))/(48.*cmath.pi*abs(MZ)**3)',
                                  (P.e__minus__,P.e__plus__):'(MZ**2*(-(ee**2*MZ**2) + (cw**2*ee**2*MZ**2)/(2.*sw**2) + (5*ee**2*MZ**2*sw**2)/(2.*cw**2)))/(48.*cmath.pi*abs(MZ)**3)',
                                  (P.mu__minus__,P.mu__plus__):'(MZ**2*(-(ee**2*MZ**2) + (cw**2*ee**2*MZ**2)/(2.*sw**2) + (5*ee**2*MZ**2*sw**2)/(2.*cw**2)))/(48.*cmath.pi*abs(MZ)**3)',
                                  (P.ta__minus__,P.ta__plus__):'((-5*ee**2*MTA**2 - ee**2*MZ**2 - (cw**2*ee**2*MTA**2)/(2.*sw**2) + (cw**2*ee**2*MZ**2)/(2.*sw**2) + (7*ee**2*MTA**2*sw**2)/(2.*cw**2) + (5*ee**2*MZ**2*sw**2)/(2.*cw**2))*cmath.sqrt(-4*MTA**2*MZ**2 + MZ**4))/(48.*cmath.pi*abs(MZ)**3)'})

Decay_W__plus__ = Decay(name = 'Decay_W__plus__',
                        particle = P.W__plus__,
                        partial_widths = {(P.u,P.d__tilde__):'(ee**2*MW**4)/(16.*cmath.pi*sw**2*abs(MW)**3)',
                                          (P.c,P.s__tilde__):'(ee**2*MW**4)/(16.*cmath.pi*sw**2*abs(MW)**3)',
                                          (P.t,P.b__tilde__):'(((-3*ee**2*MB**2)/(2.*sw**2) - (3*ee**2*MT**2)/(2.*sw**2) - (3*ee**2*MB**4)/(2.*MW**2*sw**2) + (3*ee**2*MB**2*MT**2)/(MW**2*sw**2) - (3*ee**2*MT**4)/(2.*MW**2*sw**2) + (3*ee**2*MW**2)/sw**2)*cmath.sqrt(MB**4 - 2*MB**2*MT**2 + MT**4 - 2*MB**2*MW**2 - 2*MT**2*MW**2 + MW**4))/(48.*cmath.pi*abs(MW)**3)',
                                          (P.ve,P.e__plus__):'(ee**2*MW**4)/(48.*cmath.pi*sw**2*abs(MW)**3)',
                                          (P.vm,P.mu__plus__):'(ee**2*MW**4)/(48.*cmath.pi*sw**2*abs(MW)**3)',
                                          (P.vt,P.ta__plus__):'((-MTA**2 + MW**2)*(-(ee**2*MTA**2)/(2.*sw**2) - (ee**2*MTA**4)/(2.*MW**2*sw**2) + (ee**2*MW**2)/sw**2))/(48.*cmath.pi*abs(MW)**3)'})

Decay_b = Decay(name = 'Decay_b',
                particle = P.b,
                partial_widths = {(P.W__minus__,P.t):'(((3*ee**2*MB**2)/(2.*sw**2) + (3*ee**2*MT**2)/(2.*sw**2) + (3*ee**2*MB**4)/(2.*MW**2*sw**2) - (3*ee**2*MB**2*MT**2)/(MW**2*sw**2) + (3*ee**2*MT**4)/(2.*MW**2*sw**2) - (3*ee**2*MW**2)/sw**2)*cmath.sqrt(MB**4 - 2*MB**2*MT**2 + MT**4 - 2*MB**2*MW**2 - 2*MT**2*MW**2 + MW**4))/(96.*cmath.pi*abs(MB)**3)'})

Decay_t = Decay(name = 'Decay_t',
                particle = P.t,
                partial_widths = {(P.W__plus__,P.b):'(((3*ee**2*MB**2)/(2.*sw**2) + (3*ee**2*MT**2)/(2.*sw**2) + (3*ee**2*MB**4)/(2.*MW**2*sw**2) - (3*ee**2*MB**2*MT**2)/(MW**2*sw**2) + (3*ee**2*MT**4)/(2.*MW**2*sw**2) - (3*ee**2*MW**2)/sw**2)*cmath.sqrt(MB**4 - 2*MB**2*MT**2 + MT**4 - 2*MB**2*MW**2 - 2*MT**2*MW**2 + MW**4))/(96.*cmath.pi*abs(MT)**3)'})

Decay_ta__minus__ = Decay(name = 'Decay_ta__minus__',
                          particle = P.ta__minus__,
                          partial_widths = {(P.W__minus__,P.vt):'((MTA**2 - MW**2)*((ee**2*MTA**2)/(2.*sw**2) + (ee**2*MTA**4)/(2.*MW**2*sw**2) - (ee**2*MW**2)/sw**2))/(32.*cmath.pi*abs(MTA)**3)'})

