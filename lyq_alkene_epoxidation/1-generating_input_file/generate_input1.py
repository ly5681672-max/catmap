# 生成 CatMAP 输入文件：计算形成能、输出 TXT/CSV，并验证解析结果。
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
    'Ha_111': -873.9655615,  # H 在位置 a
    'Hb_111': -875.018362,   # H 在位置 b
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

# 预计算原子参考能，单位：eV
ref_dict = {
    'O': -4.37,
    'C': -9.28,
    'H': -1.11,
    '111': abinitio_energies['slab_111'],
}

# 特殊物种名称 -> 实际化学组成。
# Ha/Hb 表示不同吸附位置上的 H，不是化学元素 Ha/Hb。
formula_map = {
    'Ha': 'H',
    'Hb': 'H',
    'Hb-O': 'HO',
}


def get_formation_energies(energy_dict, references):
    """根据原子参考能和裸板能计算形成能。"""
    formation_energies = {}

    for key, energy in energy_dict.items():
        if '_' not in key:
            raise ValueError(f'物种 key 缺少 site 信息: {key}')

        name, site = key.split('_', 1)

        if 'slab' in name:
            continue

        E0 = energy

        # 表面吸附态/过渡态先减去裸板能；气相不减裸板。
        if site == '111':
            E0 -= references['111']

        # 特殊名称显式映射，普通名称则移除 TS 连字符后按化学式解析。
        formula = formula_map.get(name, name.replace('-', ''))

        try:
            composition = string2symbols(formula)
        except Exception as exc:
            raise ValueError(
                f'无法解析物种 {name}，对应化学式 {formula}，原 key: {key}'
            ) from exc

        for atom in composition:
            if atom not in references:
                raise KeyError(f'物种 {key} 中的原子 {atom} 未在 ref_dict 中定义参考能')
            E0 -= references[atom]

        formation_energies[key] = round(E0, 3)

    return formation_energies


formation_energies = get_formation_energies(abinitio_energies, ref_dict)

# 打印形成能检查
for key, value in formation_energies.items():
    print(f'{key} {value}')

frequency_dict = {
    'H2O2_gas': [],
    'C6H12_gas': [],
    'C6H12O_gas': [],
    'H2O_gas': [],
    'H2O2_111': [],
    'Ha_111': [],
    'Hb_111': [],
    'OOH_111': [],
    'C6H12_111': [],
    'C6H12O_111': [],
    'O_111': [],
    'OH_111': [],
    'H2O_111': [],
    'OOH-C6H12_111': [],
    'Hb-O_111': [],
    'slab_111': [],
}


def make_input_file(file_name, energy_dict, frequencies):
    """生成 CatMAP TableParser 使用的制表符分隔输入文件。"""
    header = '\t'.join([
        'surface_name',
        'site_name',
        'species_name',
        'formation_energy',
        'frequencies',
        'reference',
    ])

    lines = []

    for key, energy in energy_dict.items():
        if '_' not in key:
            continue

        name, site = key.split('_', 1)

        if 'slab' in name:
            continue

        frequency = frequencies.get(key, [])
        surface = None if site == 'gas' else 'tife'
        outline = [surface, site, name, energy, frequency, 'Input File Tutorial.']
        lines.append('\t'.join(str(value) for value in outline))

    lines.sort()
    input_file = '\n'.join([header] + lines)

    with open(file_name, 'w', encoding='utf-8') as file:
        file.write(input_file)

    print(f'Successfully created input file: {file_name}')


def make_csv_file(file_name, energy_dict, frequencies):
    """将同一套形成能数据输出为 CSV，便于表格软件检查。"""
    fieldnames = [
        'surface_name',
        'site_name',
        'species_name',
        'formation_energy',
        'frequencies',
        'reference',
    ]

    with open(file_name, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for key, energy in energy_dict.items():
            if '_' not in key:
                continue

            name, site = key.split('_', 1)

            if 'slab' in name:
                continue

            freq_list = frequencies.get(key, [])
            freq_str = ';'.join(
                str(int(freq))
                if isinstance(freq, (int, float)) and float(freq).is_integer()
                else str(freq)
                for freq in freq_list
            ) if freq_list else ''

            writer.writerow({
                'surface_name': None if site == 'gas' else 'tife',
                'site_name': site,
                'species_name': name,
                'formation_energy': energy,
                'frequencies': freq_str,
                'reference': 'Input File Tutorial.',
            })

    print(f'Successfully created CSV file: {file_name}')


# 生成输出文件
txt_file_name = 'tife_energies.txt'
csv_file_name = 'tife_energies.csv'
make_input_file(txt_file_name, formation_energies, frequency_dict)
make_csv_file(csv_file_name, formation_energies, frequency_dict)

# CatMAP 解析测试
try:
    from catmap.model import ReactionModel
    from catmap.parsers import TableParser

    rxm = ReactionModel()
    rxm.surface_names = ['tife']
    rxm.adsorbate_names = (
        'H2O2', 'Ha', 'Hb', 'OOH', 'C6H12', 'C6H12O', 'O', 'OH', 'H2O'
    )
    rxm.transition_state_names = ('OOH-C6H12', 'Hb-O')
    rxm.gas_names = ('H2O2_g', 'C6H12_g', 'C6H12O_g', 'H2O_g')
    rxm.site_names = ('s',)

    # CatMAP 对 Ha/Hb/Hb-O 无法从名称自动推断元素组成，因此显式指定。
    rxm.species_definitions = {
        's': {'site_names': ['111']},
        'Ha': {'composition': {'H': 1}},
        'Hb': {'composition': {'H': 1}},
        'Hb-O': {'composition': {'H': 1, 'O': 1}},
    }

    parser = TableParser(rxm)
    parser.input_file = txt_file_name
    parser.parse()

    print('CatMAP parsing successful.')
    for key, value in rxm.species_definitions.items():
        print(f'{key} {value}')

except Exception as exc:
    print('CatMAP 解析测试失败：')
    print(exc)
    # 不再吞掉异常：解析失败时让进程返回非 0 退出码，便于定位问题。
    raise
