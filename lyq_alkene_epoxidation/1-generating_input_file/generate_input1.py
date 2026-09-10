# 将 data_only_not_executable.py 的体系和数值套用到本脚本中，修复键名并添加 CSV 输出及容错解析。
from ase.symbols import string2symbols
import csv

abinitio_energies = {
        # ========== 气相物种 ==========
    'H2O2_gas': -18.122091,
    'C6H12_gas': -95.031223,
    'C6H12O_gas': -101.037771,
    'H2O_gas': -14.204457,
         # ========== 吸附物种 ==========
     'H2O2_111': -888.055199,
     'Ha_111': -873.9655615,  # H在位置1
     'Hb_111': -875.018362,  # H在位置2
     'OOH_111': -883.605234,
     'C6H12_111': -966.14116,
     'C6H12O_111': -972.2223845,
     'O_111': -874.6632055,
     'OH_111': -880.2182855,
     'H2O_111': -885.0369905,
         # ========== 过渡态 ==========
     'OOH-C6H12_111': -980.1874005,
     'Hb-O_111': -879.6075115,
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
    """
    计算形成能（基于给定的参考能 ref_dict）。会移除 transition-state 名称中的 '-' 后按元素逐个扣除参考能。
    """
    formation_energies = {}
    for key in energy_dict.keys():
        E0 = energy_dict[key]
        # 保护性处理：确保 key 中有下划线分隔 site
        if '_' not in key:
            # 如果没有 site 信息则跳过
            continue
        name, site = key.split('_', 1)
        if 'slab' not in name:
            # 对吸附态需减去 slab 能
            if site == '111':
                E0 -= ref_dict['111']
            # 移除转移态中的连字符，便于解析组成
            formula = name.replace('-', '')
            # 使用 ase.symbols.string2symbols 获取组成列表（会将字母序列解析为元素符号）
            try:
                composition = string2symbols(formula)
            except Exception:
                # 如果解析失败，跳过该项并打印警告
                print(f'Warning: 无法解析化学式 "{formula}"，跳过该项的形成能计算。原 key: {key}')
                continue
            for atom in composition:
                if atom in ref_dict:
                    E0 -= ref_dict[atom]
                else:
                    # 如果遇到未知原子种类，打印提示并跳过该原子
                    print(f'Warning: 未在 ref_dict 中找到原子 {atom}，在项 {key} 中跳过该原子扣除。')
            E0 = round(E0, 3)
            formation_energies[key] = E0
    return formation_energies


formation_energies = get_formation_energies(abinitio_energies, ref_dict)

# 打印形成能检查
for key in formation_energies:
    print(f'{key} {formation_energies[key]}')

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


def make_input_file(file_name, energy_dict, frequency_dict):
    """
    生成制表符分隔的 energies_origin.txt（与原始行为一致），并保证 frequency 字段存在。
    """
    header = '\t'.join(['surface_name', 'site_name', 'species_name', 'formation_energy', 'frequencies', 'reference'])
    lines = []
    for key in energy_dict.keys():
        E = energy_dict[key]
        if '_' not in key:
            continue
        name, site = key.split('_', 1)
        if 'slab' not in name:
            # 为了防止 KeyError，使用 dict.get 并在缺失时使用空列表
            frequency = frequency_dict.get(key, [])
            surface = None if site == 'gas' else 'tife'
            outline = [surface, site, name, E, frequency, 'Input File Tutorial.']
            line = '\t'.join([str(w) for w in outline])
            lines.append(line)
    lines.sort()
    lines = [header] + lines
    input_file = '\n'.join(lines)
    with open(file_name, 'w', encoding='utf-8') as f:
        f.write(input_file)
    print(f'Successfully created input file: {file_name}')


def make_csv_file(file_name, energy_dict, frequency_dict):
    """
    将数据写为 CSV：surface_name, site_name, species_name, formation_energy, frequencies, reference
    frequencies 会被转为字符串形式（例如 [3657, 1595] -> "3657;1595"）以便在表格软件中查看。
    """
    with open(file_name, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['surface_name', 'site_name', 'species_name', 'formation_energy', 'frequencies', 'reference']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for key in energy_dict.keys():
            E = energy_dict[key]
            if '_' not in key:
                continue
            name, site = key.split('_', 1)
            if 'slab' in name:
                continue
            freq_list = frequency_dict.get(key, [])
            # 将频率列表序列化为分号分隔的字符串，空列表写为空字符串
            freq_str = ';'.join(
                [str(int(f)) if isinstance(f, (int, float)) and float(f).is_integer() else str(f) for f in
                 freq_list]) if freq_list else ''
            surface = None if site == 'gas' else 'tife'
            writer.writerow({
                'surface_name': surface,
                'site_name': site,
                'species_name': name,
                'formation_energy': E,
                'frequencies': freq_str,
                'reference': 'Input File Tutorial.'
            })
    print(f'Successfully created CSV file: {file_name}')

# 生成输出文件
txt_file_name = 'tife_energies.txt'
csv_file_name = 'tife_energies.csv'
make_input_file(txt_file_name, formation_energies, frequency_dict)
make_csv_file(csv_file_name, formation_energies, frequency_dict)

# 以下为简单的解析测试（将解析过程包在 try/except 中，避免解析异常导致脚本退出）
try:
    from catmap.model import ReactionModel
    from catmap.parsers import TableParser

    rxm = ReactionModel()
    # 与 data_only_not_executable.py 保持一致的设置
    rxm.surface_names = ['tife']
    rxm.adsorbate_names = ('H2O2', 'Ha', 'Hb', 'OOH', 'C6H12', 'C6H12O', 'O', 'OH', 'H2O')
    rxm.transition_state_names = ('OOH-C6H12', 'Hb-O')
    rxm.gas_names = ('H2O2_g', 'C6H12_g', 'C6H12O_g', 'H2O_g')
    rxm.site_names = ('s',)
    rxm.species_definitions = {'s': {'site_names': ['111']}}
    parser = TableParser(rxm)
    parser.input_file = txt_file_name
    parser.parse()
    for key in rxm.species_definitions:
        print(f'{key} {rxm.species_definitions[key]}')
except Exception as e:
    print('解析测试出现异常（已捕获），解析步骤未能完成。异常信息如下：')
    print(e)
