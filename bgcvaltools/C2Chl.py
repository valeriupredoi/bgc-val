#
# Copyright 2014 Plymouth Marine Laboratory
#
# This file is part of the ukesm-validation library.
#
# ukesm-validation is free software: you can redistribute it and/or modify it
# under the terms of the Lesser GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, 
# or (at your option) any later version. 

# ukesm-validation is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU
# General Public License for more details.
# You should have received a copy of the Lesser GNU General
# Public License along with ukesm-validation. If not, see <http://www.gnu.org/licenses/>.
#
# Address:
# Plymouth Marine Laboratory
# Prospect Place, The Hoe
# Plymouth, PL1 3DH, UK
#
# Email:
# momm@pml.ac.ukesm
#

class C2Chl:
    """class for relationships of carbon to chlorophyll following
    Sathyendranath et al. 2009.
    The class includes the loglog functions of the from log10(C)=m+p*log10(Chl) and their transformation C=m*Chl^p as well as the carbon to chlorophyll ratio C2Chl=m*Chl^(p-1)."""
    def __init__(self,m,p):
        self.loglog= lambda x:self._loglog(x,m,p)
        self.C= lambda x:self._C(x,m,p)
        self.C2Chl= lambda x:self._C2Chl(x,m,p)

    _loglog = lambda self,lchl,m,p: m+p*lchl
    
    _C = lambda self,chl,m,p: 10**m*chl**p

    _C2Chl = lambda self,chl,m,p: 10**m*chl**(p-1)

    def __call__(self,lchl):
       return self.loglog(lchl)

#Sathyendranath HPLC (C_P):
SH = C2Chl(1.9,.65)
#Sathyendranath Turner (C_P):
ST = C2Chl(1.81,.63)
#Buck Turner:
B = C2Chl(1.92,.69)

#Sathyendranath HPLC (C_T):
SH_T = C2Chl(2.25,0.48)
#Sathyendranath Turner (C_T):
ST_T = C2Chl(2.19,0.45)

