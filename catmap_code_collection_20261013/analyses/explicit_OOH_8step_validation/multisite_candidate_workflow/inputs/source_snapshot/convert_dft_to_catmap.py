"""只做 DFT Gibbs 能量换算，不运行 CatMAP，不修改原有 .mkm/energies.txt。

在 lyq_alkene_epoxidation/1-generating_input_file/ 运行：
    python convert_dft_to_catmap.py

输出 3 个文件（均位于本脚本所在目录）：
    DFT_formation_catmap.txt   CatMAP TableParser 形成能输入（独立构型）
    DFT_adsorption.csv         对照气相参照的吸附自由能
    DFT_TS_relative.csv        Excel 共吸附 NEB 的 IS/TS/FS 相对自由能

注意：输入 TXT 中 OOH-C6H12 和 Hb-O 来源于独立态数据表中的
      候选结构，并不等同于共吸附 NEB 的 TS 自由能；后者只见相对能量 CSV。
"""
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET
import csv
import re

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / '己烯环氧化.xlsx'
A = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
REL = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'
PKG = '{http://schemas.openxmlformats.org/package/2006/relationships}'
SURFACES = {'TiFe':'tife','TiMn':'timn','TiHf':'tihf','TiRe':'tire',
            'TiNb':'tinb','TiMo':'timo','TiV':'tiv','TiZr':'tizr',
            'TiCo':'tico','TiTi':'titi','TiW':'tiw','TiTa':'tita',
            'TiCr':'ticr','单核Ti':'ti','TiNi':'tini'}
MU = {'C':-9.28, 'H':-1.11, 'O':-4.37}  # 原始 DFT 原子参考能，eV
GAS_ROWS = {'H2O2':57,'C6H12':58,'H2O':59,'C6H12O':60,'OOH':61}
SINGLE_ROWS = {'H2O2':4,'Ha':8,'Hb':12,'OOH':16,'C6H12':20,
               'C6H12O':24,'O':29,'OH':34,'H2O':39,
               'OOH-C6H12':43,'Hb-O':47}
FORMULA = {'H2O2':{'H':2,'O':2},'Ha':{'H':1},'Hb':{'H':1},
           'OOH':{'H':1,'O':2},'C6H12':{'C':6,'H':12},
           'C6H12O':{'C':6,'H':12,'O':1},'O':{'O':1},
           'OH':{'O':1,'H':1},'H2O':{'H':2,'O':1},
           'OOH-C6H12':{'C':6,'H':13,'O':2},'Hb-O':{'H':1,'O':1}}

def load_sheets(path):
    """仅用 Python 标准库读取 XLSX 中已计算好的 G 单元格。"""
    with ZipFile(path) as z:
        workbook=ET.fromstring(z.read('xl/workbook.xml'))
        rels=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        targets={r.attrib['Id']:r.attrib['Target']
                 for r in rels.findall(PKG+'Relationship')}
        shared=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            xml=ET.fromstring(z.read('xl/sharedStrings.xml'))
            shared=[''.join(t.text or '' for t in node.iter(A+'t'))
                    for node in xml.findall(A+'si')]
        out={}
        for sh in workbook.find(A+'sheets').findall(A+'sheet'):
            title=sh.attrib['name']
            if title not in ('catmap数据集','吉布斯自由能汇总'):
                continue
            member=targets[sh.attrib[REL+'id']].lstrip('/')
            if not member.startswith('xl/'):
                member='xl/'+member
            xml=ET.fromstring(z.read(member))
            cells={}
            for cell in xml.iter(A+'c'):
                value=cell.find(A+'v')
                if value is None or value.text is None:
                    continue
                v=(shared[int(value.text)] if cell.attrib.get('t')=='s'
                   else float(value.text))
                cells[cell.attrib['r']]=v
            out[title]=cells
        if len(out)!=2:
            raise ValueError('找不到 catmap数据集 / 吉布斯自由能汇总 工作表')
        return out

def cell(col,row):
    letters=''
    while col:
        col,k=divmod(col-1,26)
        letters=chr(65+k)+letters
    return f'{letters}{row}'

def surfaces(table):
    """通过第一行催化剂名、第三行 G 表头定位每组的 Gibbs 列。"""
    positions=[]
    for addr,v in table.items():
        if re.fullmatch(r'[A-Z]+1',addr) and v in SURFACES:
            letters=addr[:-1];num=0
            for ch in letters:num=num*26+ord(ch)-64
            positions.append((num,SURFACES[v]))
    positions.sort()
    if len(positions)!=15:raise ValueError('催化剂数应为 15')
    out={}
    for index,(start,name) in enumerate(positions):
        stop=positions[index+1][0] if index+1<len(positions) else start+9
        gcols=[c for c in range(start,stop) if table.get(cell(c,3))=='G']
        ecols=[c for c in range(start,stop) if table.get(cell(c,3))=='E0']
        if len(gcols)!=1 or len(ecols)!=1:
            raise ValueError(f'{name} 的 E0/G 列无法唯一定位')
        out[name]=(ecols[0],gcols[0])
    return out

def ref_energy(composition):
    return sum(MU[atom]*count for atom,count in composition.items())

def save_csv(path,header,rows):
    with path.open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(header);w.writerows(rows)

def main():
    tables=load_sheets(SOURCE)
    singles=tables['catmap数据集'];summary=tables['吉布斯自由能汇总']
    cols1=surfaces(singles);cols2=surfaces(summary)
    get=lambda table,col,row:float(table[cell(col,row)])
    gas={species:get(summary,cols2['tife'][1],row)
         for species,row in GAS_ROWS.items()}
    gas_e0={species:get(summary,cols2['tife'][0],row)
            for species,row in GAS_ROWS.items()}
    if abs(gas['C6H12']+95.031223)>1e-5:
        raise ValueError('气相自由能与已审核源数据不一致，请检查 Excel')
    out=['surface_name\tsite_name\tspecies_name\tformation_energy\tfrequencies\treference']
    adsorption=[];neb=[]
    for species in ('H2O2','C6H12','C6H12O','H2O'):
        energy=gas[species]-ref_energy(FORMULA[species])
        out.append(f'None\tgas\t{species}\t{energy:.9f}\t[]\tDFT_Gibbs_gas')
    for surface in SURFACES.values():
        e1,c1=cols1[surface]
        e2,c2=cols2[surface]
        slab=get(singles,c1,51)
        slab_e0=get(singles,e1,51)
        for species,row in SINGLE_ROWS.items():
            g=get(singles,c1,row)
            formation=g-slab-ref_energy(FORMULA[species])
            provenance='DFT_Gibbs_TS_candidate' if species in ('OOH-C6H12','Hb-O') else 'DFT_Gibbs_single'
            out.append(f'{surface}\t111\t{species}\t{formation:.9f}\t[]\t{provenance}')
            if species in gas:
                e0=get(singles,e1,row)
                adsorption.append([surface,species,round(e0,9),round(slab_e0,9),
                                   round(gas_e0[species],9),round(e0-slab_e0-gas_e0[species],9),
                                   round(g,9),round(slab,9),round(gas[species],9),
                                   round(g-slab-gas[species],9)])
        g3,gts1,g4=(get(summary,c2,r) for r in (12,16,20))
        g5,gts2,g6=(get(summary,c2,r) for r in (24,28,32))
        # TiMo: row 20 is adsorbate-only with gas epoxide omitted;
        # row 21 is the original workbook's total including gaseous epoxide.
        fs1_total=get(summary,c2,21) if surface=='timo' else g4
        # TiMo: row 32 includes gaseous epoxide, unlike row 24/28.
        fs2_total=g6-gas['C6H12O'] if surface=='timo' else g6
        neb.append([surface,*[round(v,9) for v in (g3,gts1,fs1_total,g5,gts2,fs2_total)],
                    round(gts1-g3,9),round(gts1-fs1_total,9),
                    round(gts2-g5,9),round(gts2-fs2_total,9),
                    round(get(summary,e2,16)-get(summary,e2,12),9),
                    round(get(summary,e2,28)-get(summary,e2,24),9),
                    'TiMo: Excel totals / gas epoxide normalization' if surface=='timo'
                    else 'Original Excel Gibbs rows 12,16,20 / 24,28,32'])
    (HERE/'DFT_formation_catmap.txt').write_text('\n'.join(out)+'\n',encoding='ascii')
    save_csv(HERE/'DFT_adsorption.csv',
             ['surface','species','E0_ads_eV','E0_slab_eV','E0_gas_eV','Eads_eV',
              'G_ads_eV','G_slab_eV','G_gas_eV','dGads_eV'],adsorption)
    save_csv(HERE/'DFT_TS_relative.csv',
             ['surface','G_IS1','G_TS1','G_FS1_total','G_IS2','G_TS2','G_FS2_adjusted',
              'dG_TS1_minus_IS1','dG_TS1_minus_FS1','dG_TS2_minus_IS2',
              'dG_TS2_minus_FS2','dE0_TS1_minus_IS1','dE0_TS2_minus_IS2',
              'reference_note'],neb)
    print(f'转换完成：15组催化剂；{len(out)-1}条 CatMAP 形成能数据；'
          f'{len(adsorption)}条吸附能；{len(neb)}组 NEB 相对能量。')
    print('没有修改原 Excel、energies.txt、.mkm；没有运行微观动力学。')
    print('提示：单独构型的 TS 候选形成能与共吸附 NEB 的 TS-IS 能垒不是同一种数据。')

if __name__=='__main__':
    main()
