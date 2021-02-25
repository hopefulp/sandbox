###### plot params following jobs
try:
    from Ni_CO2red import *
except:
    print("import module python/qchem")

class JobParams():
    def __init__(self, title=None, xtitle=None, ytitle=None):
        self.title=title
        self.xtitle=xtitle
        self.ytitle=ytitle


def get_jobtitle(job, title=None, xtitle=None, ytitle=None):

    if job == 'qcmo':
        if title==None:
            title="MO Plot"
        if xtitle==None:
            xtitle="Models"
        if ytitle==None:
            ytitle="Charges"
    elif job == 'ai':
        if title==None:
            title="filename"
        if xtitle==None:
            xtitle='t(ps)'
        if ytitle==None:
            ytitle='E(kJ/mol)'
    
    return JobParams(title, xtitle, ytitle)

#Ni_files=[m1.f3, m2.f3, m3.f3, m4.f3, m5.f3, m6.f3]
Ni_6x=['1-PP', '2-PP$^{Me}$P', '3-PN$^{Me}$P', '4-PNP', '5-Ni-CO$_2$-Fe', '6-Ni-CO$_2$-Fe$^{SM}$']
Ni_5x=['1-PP', '2-PP$^{Me}$P', "2`-PN$^{Me}$P", '3-PNP', '4-Ni-CO2-Fe']
Ni_4x=['1PP', '2PPP', '3PNP', '4PNP-bridged']



