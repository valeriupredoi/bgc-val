#
# Here, I load a few random samples, in order to test the robust stats diagrams,
#

import numpy as np
from matplotlib import pyplot
from bgcvaltools import robust, RobustStatistics,StatsDiagram

a = np.random.random_integers(0,100,50)
b = np.random.random_integers(0,500,50)
b = (a*100.+b-250)/100.




print "Robust:"
rstats = robust.StatsDiagram(a,b,1.)
print "Standard:"
sstats = StatsDiagram.StatsDiagram(a,b)


fig1 = pyplot.figure()
#fig.add_subplot(211)

rTarget=robust.TargetDiagram(rstats.gamma,rstats.E0, rstats.E, rstats.R, marker = 'o',s=150, label="Robust",)

fig2 = pyplot.figure()
sTarget=StatsDiagram.TargetDiagram(sstats.gamma,sstats.E0, sstats.R, label="Standard",)

fig3 = pyplot.figure()
pyplot.scatter(a,b)

fig1.show()
fig2.show()

fig3.show()
#self,gam,E0,E,rho,marker='o',s=40,antiCorrelation=False,*opts,**keys):

