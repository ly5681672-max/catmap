from ase.symbols import string2symbols

abinitio_energies = {
        
        'C6H14_gas': -6444.26783629,
        'C2H4_gas': -2136.090987, 
        'C2H6_gas': -1028.45297997, 
        'C6H14_111':-341545.455848095, 
        'C6H13_111':-341529.614200962,
        'H_111': -335116.12846635,
        'C4H9_111':-339392.045398202,
        'C2H4_111':-339392.19,
        'C2H5_111':-337254.656307692,
        'C4H9-C2H4_111':-341525.218267321,
        'C2H5-C2H4_111':-339388.420864213,
        'C2H5-H_111':-337268.500243499,
        'slab_111': -335100.56850075,
        }

ref_dict = {}
ref_dict['H'] = (abinitio_energies['C2H6_gas'] - abinitio_energies['C2H4_gas'])*0.5
ref_dict['C'] = (abinitio_energies['C2H4_gas'] - 4*ref_dict['H'])*0.5
ref_dict['111'] = abinitio_energies['slab_111']

def get_formation_energies(energy_dict,ref_dict):
    formation_energies = {}
    for key in energy_dict.keys(): #iterate through keys
        E0 = energy_dict[key] #raw energy
        name,site = key.split('_') #split key into name/site
        if 'slab' not in name: #do not include empty site energy (0)
            if site == '111':
                E0 -= ref_dict[site] #subtract slab energy if adsorbed
            #remove - from transition-states
            formula = name.replace('-','')
            #get the composition as a list of atomic species
            composition = string2symbols(formula)
            #for each atomic species, subtract off the reference energy
            for atom in composition:
                E0 -= ref_dict[atom]
            #round to 3 decimals since this is the accuracy of DFT
            E0 = round(E0,3)
            formation_energies[key] = E0
    return formation_energies

formation_energies = get_formation_energies(abinitio_energies,ref_dict)

for key in formation_energies:
    print(str(key) + ' ' + str(formation_energies[key]))

frequency_dict = {
                
                'C6H14_gas': [],
                'C2H4_gas': [],
                'C2H6_gas': [],
                'C6H14_111':[], 
                'C6H13_111':[],
                'H_111': [],
                'C4H9_111':[],
                #'C2H4_111':[],
                'C2H5_111':[],
                'C4H9-C2H4_111':[],
                'C2H5-C2H4_111':[],
                'C2H5-H_111':[],
                }


def make_input_file(file_name,energy_dict,frequency_dict):

    #create a header
    header = '\t'.join(['surface_name','site_name',
                        'species_name','formation_energy',
                        'frequencies','reference'])

    lines = [] #list of lines in the output
    for key in energy_dict.keys(): #iterate through keys
        E = energy_dict[key] #raw energy
        name,site = key.split('_') #split key into name/site
        if 'slab' not in name: #do not include empty site energy (0)
            frequency = frequency_dict[key]
            if site == 'gas':
                surface = None
            else:
                surface = 'Rh'
            outline = [surface,site,name,E,frequency,'Input File Tutorial.']
            line = '\t'.join([str(w) for w in outline])
            lines.append(line)

    lines.sort() #The file is easier to read if sorted (optional)
    lines = [header] + lines #add header to top
    input_file = '\n'.join(lines) #Join the lines with a line break

    input = open(file_name,'w') #open the file name in write mode
    input.write(input_file) #write the text
    input.close() #close the file

    print('Successfully created input file')

file_name = 'energies.txt'
make_input_file(file_name,formation_energies,frequency_dict)

#Test that input is parsed correctly
from catmap.model import ReactionModel
from catmap.parsers import TableParser
rxm = ReactionModel()
#The following lines are normally assigned by the setup_file
#and are thus not usually necessary.
rxm.surface_names = ['MgO']
rxm.adsorbate_names = ('CO','C','O','H','CH','OH','CH2','CH3') 
rxm.transition_state_names = ('C-O','H-OH','H-C')
rxm.gas_names = ('CO_g','H2_g','CH4_g','H2O_g')
rxm.site_names = ('s',)
rxm.species_definitions = {'s':{'site_names':['111']}}
#Now we initialize a parser instance (also normally done by setup_file)
parser = TableParser(rxm)
parser.input_file = file_name
parser.parse()
#All structured data is stored in species_definitions; thus we can
#check that the parsing was successful by ensuring that all the
#data in the input file was collected in this dictionary.
for key in rxm.species_definitions:
    print(str(key) + ' ' + str(rxm.species_definitions[key]))
