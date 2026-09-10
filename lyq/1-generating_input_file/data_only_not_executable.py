from ase.symbols import string2symbols
import csv

#abinitio_energies 是
#1.绝对DFT能量，不考虑温度，后续需要频率 2.G吉布斯自由能，后续不需要频率
# 单位是eV C2H5_111是表面吸附物质 C4H9-C2H4_111是过渡态 'slab_111':Pt/MOX-1
abinitio_energies = {
        
        'CH2OHCHOHCH2OH_gas': -74.215128,
        'O2_gas': -9.268161,
        'H2O_gas': -14.204158,
        'CH2OHCHOHCOOH_gas':-74.668974,
        'CH2OHCHOHCH2OH_111':-888.536937,
        'O2_111': -816.31943,
        'H2O_11="111':-818.646247,
        'OH_111':-815.434485,
        'O_111': -811.318878,
        'CH2OHCHOHCH2O_111':-880.530884,
        'CH2OHCHOHCHO_111':-880.172371,
        'H_111':-816.443571,
        'CH2OHCHOHCHOOH_111':-887.550389,
        'CH2OHCHOHCOOH_111':-883.813514,
        'O2-2OH_111':-830.335421,
        'CH2OHCHOHCH2OH-CH2OHCHOHCH2O_111':-893.39887,
        'CH2OHCHOHCH2O-CH2OHCHOHCHO_111':-879.804246,
        'CH2OHCHOHCHO-CH2OHCHOHCHOOH_111':-887.517083,
        'CH2OHCHOHCHOOH-CH2OHCHOHCOOH_111':-886.830805,
        'slab_111': -804.838953,

        }
# 直接使用预计算好原子能量的数值，单位：eV
ref_dict = {}
ref_dict['O'] = -4.37
ref_dict['C'] = -9.28
ref_dict['H'] = -1.11
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

                'CH2OHCHOHCH2OH_gas':[],
                'O2_gas':[],
                'H2O_gas':[],
                'CH2OHCHOHCOOH_gas':[],
                'CH2OHCHOHCH2OH_111':[],
                'O2_111':[],
                'H2O_111':[],
                'OH_111':[],
                'O_111':[],
                'CH2OHCHOHCH2O_111':[],
                'CH2OHCHOHCHO_111':[],
                'H_111':[],
                'CH2OHCHOHCHOOH_111':[],
                'CH2OHCHOHCOOH_111':[],
                'O2-2OH_111':[],
                'CH2OHCHOHCH2OH-CH2OHCHOHCH2O_111':[],
                'CH2OHCHOHCH2O-CH2OHCHOHCHO_111':[],
                'CH2OHCHOHCHO-CH2OHCHOHCHOOH_111':[],
                'CH2OHCHOHCHOOH-CH2OHCHOHCOOH_111':[],
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
                surface = 'ceo2'
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

file_name = 'ceo2_energies.txt'
make_input_file(file_name,formation_energies,frequency_dict)

#Test that input is parsed correctly
from catmap.model import ReactionModel
from catmap.parsers import TableParser
rxm = ReactionModel()
#The following lines are normally assigned by the setup_file
#and are thus not usually necessary.
rxm.surface_names = ['ceo2']
rxm.adsorbate_names = ('CH2OHCHOHCH2OH','O2','H2O','OH','O','CH2OHCHOHCH2O','CH2OHCHOHCHO','H','CH2OHCHOHCHOOH','CH2OHCHOHCOOH')
rxm.transition_state_names = ('O-O','O-H3','C-H4','C-OH5','C-H6')
rxm.gas_names = ('CH2OHCHOHCH2OH_g','O2_g','H2O_g','CH2OHCHOHCOOH_g')
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
