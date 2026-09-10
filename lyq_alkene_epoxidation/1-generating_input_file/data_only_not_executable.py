from ase.symbols import string2symbols

# abinitio_energies：
# 1. 绝对 DFT 能量（不考虑温度，后续可加入频率修正）
# 2. 或已经得到的吉布斯自由能（后续无需再加频率）
# 单位：eV
abinitio_energies = {
    # ========== 气相物种 ==========
    'H2O2_gas': -18.122091,
    'C6H12_gas': -95.031223,
    'C6H12O_gas': -101.037771,
    'H2O_gas': -14.204457,

    # ========== 吸附物种 ==========
    'H2O2_111': -888.055199,
    'Ha_111': -873.9655615,   # H 在位置 a
    'Hb_111': -875.018362,    # H 在位置 b
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

# 原子参考能，单位：eV
ref_dict = {
    'O': -4.37,
    'C': -9.28,
    'H': -1.11,
    '111': abinitio_energies['slab_111'],
}

# Ha/Hb 是不同吸附位置上的 H，不是标准化学元素符号。
# 因此显式指定这些特殊物种对应的真实化学组成。
formula_map = {
    'Ha': 'H',
    'Hb': 'H',
    'Hb-O': 'HO',
}


def get_formation_energies(energy_dict, references):
    """根据裸板能和原子参考能计算形成能。"""
    formation_energies = {}

    for key, energy in energy_dict.items():
        if '_' not in key:
            raise ValueError(f'物种 key 缺少 site 信息: {key}')

        # 只按第一个下划线拆分，避免物种名中出现其他字符时拆分过多。
        name, site = key.split('_', 1)

        if 'slab' in name:
            continue

        E0 = energy

        # 表面吸附态和过渡态减去裸板能，气相不减。
        if site == '111':
            E0 -= references['111']

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

        # 必须只拆分一次；不能使用 key.split('_')。
        name, site = key.split('_', 1)

        if 'slab' in name:
            continue

        frequency = frequencies.get(key, [])
        surface = None if site == 'gas' else 'tife'
        outline = [
            surface,
            site,
            name,
            energy,
            frequency,
            'Input File Tutorial.',
        ]
        lines.append('\t'.join(str(value) for value in outline))

    lines.sort()
    input_file = '\n'.join([header] + lines)

    with open(file_name, 'w', encoding='utf-8') as file:
        file.write(input_file)

    print(f'Successfully created input file: {file_name}')


file_name = 'tife_energies.txt'
make_input_file(file_name, formation_energies, frequency_dict)


# 测试 CatMAP 是否能够正确解析生成的输入文件。
from catmap.model import ReactionModel
from catmap.parsers import TableParser

rxm = ReactionModel()
rxm.surface_names = ['tife']

# CatMAP 内部表面物种 key 使用 *_s 形式；输入表中的 species_name 仍保持
# H2O2、Ha、Hb 等名称，由 TableParser 根据 site_name=111 匹配到 s 位点。
rxm.adsorbate_names = (
    'H2O2_s',
    'Ha_s',
    'Hb_s',
    'OOH_s',
    'C6H12_s',
    'C6H12O_s',
    'O_s',
    'OH_s',
    'H2O_s',
)
rxm.transition_state_names = ('OOH-C6H12_s', 'Hb-O_s')
rxm.gas_names = ('H2O2_g', 'C6H12_g', 'C6H12O_g', 'H2O_g')
rxm.site_names = ('s',)

# 非标准名称必须显式定义元素组成。
rxm.species_definitions = {
    's': {'site_names': ['111']},
    'Ha_s': {'composition': {'H': 1}},
    'Hb_s': {'composition': {'H': 1}},
    'Hb-O_s': {'composition': {'H': 1, 'O': 1}},
}

parser = TableParser(rxm)
parser.input_file = file_name
parser.parse()

print('CatMAP parsing successful.')
for key, value in rxm.species_definitions.items():
    print(f'{key} {value}')
