from ase.symbols import string2symbols
import csv

#abinitio_energies 是
#1.绝对DFT能量，不考虑温度，后续需要频率 2.G吉布斯自由能，后续不需要频率
# 单位是eV C2H5_111是表面吸附物质 C4H9-C2H4_111是过渡态 'slab_111':Pt/MOX-1
abinitio_energies = {
        # ========== 气相物种 ==========
    'H2O2_gas': -18.122091,
    'C6H12_gas': -95.031223,
    'C6H12O_gas': -101.037771,
    'H2O_gas': -14.204457,
         # ========== 吸附物种 ==========
     'H2O2_111': -888.055199,
     'H_a_111': -873.9655615,  # H在位置1
     'H_b_111': -875.018362,  # H在位置2
     'OOH_111': -883.605234,
     'C6H12_111': -966.14116,
     'C6H12O_111': -972.2223845,
     'O_111': -874.6632055,
     'OH_111': -880.2182855,
     'H2O_111': -885.0369905,
         # ========== 过渡态 ==========
     'OOH-C6H12_111': -980.1874005,
     'H_b-O_111': -879.6075115,
        # ========== 裸板 ==========
     'slab_111': -983.050573,
        }
# 直接使用预计算好原子能量的数值，单位：eV
ref_dict = {}
ref_dict['O'] = -4.37
ref_dict['C'] = -9.28
ref_dict['H'] = -1.11
ref_dict['111'] = abinitio_energies['slab_111']


def get_formation_energies(energy_dict, ref_dict):
    formation_energies = {}
    for key in energy_dict.keys():
        E0 = energy_dict[key]
        # 先按照第一个下划线拆分
        parts = key.split('_', 1)
        if len(parts) == 2:
            name, site = parts
        else:
            name = parts[0]
            site = '_'.join(parts[1:])

        if 'slab' not in name:
            if site == '111':
                E0 -= ref_dict[site]
            formula = name.replace('-', '')
            composition = string2symbols(formula)
            for atom in composition:
                E0 -= ref_dict[atom]
            E0 = round(E0, 3)
            formation_energies[key] = E0
    return formation_energies

formation_energies = get_formation_energies(abinitio_energies,ref_dict)

for key in formation_energies:
    print(str(key) + ' ' + str(formation_energies[key]))


frequency_dict = {
                    'H2O2_gas':[],
                    'C6H12_gas':[],
                    'C6H12O_gas':[],
                    'H2O_gas':[],
                    'H2O2_111':[],
                    'Ha_111':[],
                    'Hb_111':[],
                    'OOH_111':[],
                    'C6H12_111':[],
                    'C6H12O_111':[],
                    'O_111':[],
                    'OH_111':[],
                    'H2O_111':[],
                    'OOH-C6H12_111':[],
                    'Hb-O_111':[],
                    'slab_111':[],
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
                surface = 'tife'
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

file_name = 'tife_energies.txt'
make_input_file(file_name,formation_energies,frequency_dict)

#Test that input is parsed correctly
from catmap.model import ReactionModel
from catmap.parsers import TableParser
rxm = ReactionModel()
#The following lines are normally assigned by the setup_file
#and are thus not usually necessary.
rxm.surface_names = ['tife']
rxm.adsorbate_names = ('H2O2','Ha','Hb','OOH','C6H12','C6H12O','O','OH','H2O')
rxm.transition_state_names = ('OOH-C6H12','Hb-O')
rxm.gas_names = ('H2O2_g','C6H12_g','C6H12O_g','H2O_g')
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
