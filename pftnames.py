"""	This is a dictionary of all the terms that you'd need to pick out.
	pftnames['Model Name']['Functional type']['currency']

"""


class AutoVivification(dict):
    """Implementation of perl's autovivification feature.
    	This class allows you to automate the creating of layered dictionaries.
    """
    def __getitem__(self, item):
        try: return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
            
            
pftnames = AutoVivification()
pftnames['ERSEM']['diatoms']['c'] = 'P1c'
pftnames['ERSEM']['diatoms']['n'] = 'P1n'
pftnames['ERSEM']['diatoms']['p'] = 'P1p'
pftnames['ERSEM']['diatoms']['f'] = 'P1f'
pftnames['ERSEM']['diatoms']['s'] = 'P1s'
pftnames['ERSEM']['diatoms']['chl'] = 'Chl1'

pftnames['ERSEM' ]['total']['chl'] = 'chl'
pftnames['MEDUSA']['total']['chl'] = 'CHL'

pftnames['MEDUSA']['diatoms']['n'] = ['']
pftnames['MEDUSA']['diatoms']['f'] = ['']
pftnames['MEDUSA']['diatoms']['s'] = ['']


# overloading different names with the same dict.
pftnames['E'] = pftnames['ERSEM']
pftnames['e'] = pftnames['ERSEM']
pftnames['ersem'] = pftnames['ERSEM']
pftnames['M'] = pftnames['MEDUSA']
pftnames['m'] = pftnames['MEDUSA']
pftnames['medusa'] = pftnames['MEDUSA']

