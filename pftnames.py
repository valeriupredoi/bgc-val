"""	This is a dictionary of all the terms that you'd need to pick out 

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
pftnames['ERSEM']['diatoms']['c'] = ['P1c']
pftnames['ERSEM']['diatoms']['n'] = ['P1c']
pftnames['ERSEM']['diatoms']['p'] = ['P1c']
pftnames['ERSEM']['diatoms']['f'] = ['P1c']
pftnames['ERSEM']['diatoms']['s'] = ['P1c']

pftnames['MEDUSA']['diatoms']['n'] = ['']
pftnames['MEDUSA']['diatoms']['f'] = ['']
pftnames['MEDUSA']['diatoms']['s'] = ['']


pftnames['E'] = pftnames['ERSEM']

pftnames['M'] = pftnames['MEDUSA']


