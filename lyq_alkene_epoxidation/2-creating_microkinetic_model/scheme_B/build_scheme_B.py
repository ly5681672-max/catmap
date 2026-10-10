"""Verify/regenerate coadsorbed Scheme B energies from the ORIGINAL XLSX.
Standard-library-only; original workbook and original A model never changed.
"""
import argparse
import csv
import io
import re
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET

HERE = Path(__file__).resolve().parent
DEFAULT_XLSX = HERE.parent.parent / '己烯环氧化.xlsx'
N='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
R='{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'
P='{http://schemas.openxmlformats.org/package/2006/relationships}'
LABELS={'TiFe':'tife','TiMn':'timn','TiHf':'tihf','TiRe':'tire',
'TiNb':'tinb','TiMo':'timo','TiV':'tiv','TiZr':'tizr','TiCo':'tico',
'TiTi':'titi','TiW':'tiw','TiTa':'tita','TiCr':'ticr','单核Ti':'ti','TiNi':'tini'}
MU={'C':-9.28,'H':-1.11,'O':-4.37}
NAMES=('S1','S2','S3','TS1','S4','S5','TS2','S6','S7')
F_H2O2={'H':2,'O':2}
F_HEX={'C':6,'H':14,'O':2}
F_WATER={'H':2,'O':1}
GAS_FORM={'H2O2':F_H2O2,'C6H12':{'C':6,'H':12},
'C6H12O':{'C':6,'H':12,'O':1},'H2O':F_WATER}

def colno(a):
    n=0
    for c in re.match('[A-Z]+',a).group():
        n=n*26+ord(c)-64
    return n

def addr(col,row):
    s=''
    while col:
        col,k=divmod(col-1,26)
        s=chr(k+65)+s
    return s+str(row)

def load(path):
    with ZipFile(path) as z:
        wb=ET.fromstring(z.read('xl/workbook.xml'))
        rel=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        paths={r.attrib['Id']:r.attrib['Target'] for r in rel.findall(P+'Relationship')}
        shared=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            root=ET.fromstring(z.read('xl/sharedStrings.xml'))
            shared=[''.join(t.text or '' for t in si.iter(N+'t')) for si in root.findall(N+'si')]
        sheets={}
        for sh in wb.find(N+'sheets').findall(N+'sheet'):
            name=sh.attrib['name']
            if name not in ('吉布斯自由能汇总','catmap数据集'):
                continue
            f=paths[sh.attrib[R+'id']].lstrip('/')
            if not f.startswith('xl/'): f='xl/'+f
            table={}
            for _,c in ET.iterparse(io.BytesIO(z.read(f)),events=('end',)):
                if c.tag!=N+'c': continue
                v=c.find(N+'v')
                if v is not None and v.text is not None:
                    if c.attrib.get('t')=='e': raise ValueError(name+':'+c.attrib['r'])
                    table[c.attrib['r']]=(shared[int(v.text)] if c.attrib.get('t')=='s' else float(v.text))
                c.clear()
            sheets[name]=table
        if len(sheets)!=2: raise ValueError('Two required XLSX sheets not found')
        return sheets

def find_cols(table):
    starts=sorted((colno(a),v) for a,v in table.items()
                  if re.fullmatch('[A-Z]+1',a) and v in LABELS)
    if len(starts)!=15: raise ValueError('Expected 15 surface columns')
    out={}
    for idx,(start,label) in enumerate(starts):
        stop=starts[idx+1][0] if idx+1<len(starts) else start+9
        choices=[i for i in range(start,stop) if table.get(addr(i,3))=='G']
        if len(choices)!=1: raise ValueError('Bad G header: '+label)
        out[label]=choices[0]
    return out

def gf(G,slab,comp):
    return G-slab-sum(MU[a]*n for a,n in comp.items())

def create(path):
    tables=load(path);a=tables['catmap数据集'];b=tables['吉布斯自由能汇总']
    ia,ib=find_cols(a),find_cols(b)
    get=lambda tab,k,r:tab[addr(k,r)]
    gasraw={x:get(b,ib['TiFe'],r) for x,r in
            [('H2O2',57),('C6H12',58),('H2O',59),('C6H12O',60)]}
    gas={g:round(gf(raw,0,GAS_FORM[g]),3) for g,raw in gasraw.items()}
    assert gas=={'H2O2':-7.162,'C6H12':-26.031,'H2O':-7.614,'C6H12O':-27.668},gas
    rec=[]
    for label,surface in LABELS.items():
        ka,kb=ia[label],ib[label]; sa=get(a,ka,51); sb=get(b,kb,4)
        v={
            'S1':gf(get(b,kb,8),sb,F_H2O2),
            'S2':gf(get(a,ka,8),sa,{'H':1})+gf(get(a,ka,16),sa,{'H':1,'O':2}),
            'S3':gf(get(b,kb,12),sb,F_HEX),
            'TS1':gf(get(b,kb,16),sb,F_HEX),
            'S4':gf(get(b,kb,21 if surface=='timo' else 20),sb,F_HEX),
            'S5':gf(get(b,kb,24),sb,F_WATER),
            'TS2':gf(get(b,kb,28),sb,F_WATER),
            'S6':gf(get(b,kb,32)-(gasraw['C6H12O'] if surface=='timo' else 0),sb,F_WATER),
            'S7':gf(get(a,ka,39),sa,F_WATER),
        }
        steps=[v['S1']-gas['H2O2'],v['S2']-v['S1'],
        v['S3']-v['S2']-gas['C6H12'],v['S4']-v['S3'],
        v['S5']+gas['C6H12O']-v['S4'],v['S6']-v['S5'],
        v['S7']-v['S6'],gas['H2O']-v['S7']]
        err=sum(steps)-(gas['C6H12O']+gas['H2O']-gas['H2O2']-gas['C6H12'])
        if v['TS1']-v['S4']<-1e-8 or abs(err)>1e-7:
            raise ValueError('TS1 below FS or cycle not closed: '+surface)
        rec.append((surface,v,steps,err))
    return rec,gas

def rows(rec,gas):
    out=['surface_name\tsite_name\tspecies_name\tformation_energy\tfrequencies\treference']
    for g in ('H2O2','C6H12','C6H12O','H2O'):
        out.append(f'None\tgas\t{g}\t{gas[g]:.3f}\t[]\tExcel:GibbsSummary/gas')
    for surf,v,_,_ in rec:
        for s in NAMES:
            source=('Excel:SingleAds/Ha+OOH_proxy' if s=='S2' else
                'Excel:SingleAds/H2O' if s=='S7' else
                'Excel:GibbsSummary/TiMo_gas_adjusted' if surf=='timo' and s in ('S4','S6') else
                'Excel:GibbsSummary/coadsorbed_G')
            out.append(f'{surf}\t111\t{s}\t{v[s]:.9f}\t[]\t{source}')
    return '\n'.join(out)+'\n'

def audit(rec):
    out=['surface,S1,S2,S3,TS1,S4,S5,TS2,S6,S7,TS1_minus_IS,TS1_minus_FS,TS2_minus_IS,TS2_minus_FS,cycle_deltaG,cycle_error,note']
    for surf,v,steps,err in rec:
        out.append(','.join([surf]+[f'{v[s]:.9f}' for s in NAMES]+[
          f'{v["TS1"]-v["S3"]:.9f}',f'{v["TS1"]-v["S4"]:.9f}',
          f'{v["TS2"]-v["S5"]:.9f}',f'{v["TS2"]-v["S6"]:.9f}',
          f'{sum(steps):.9f}',f'{err:.9f}',
          'TiMo_gas_normalized' if surf=='timo' else 'verified']))
    return '\n'.join(out)+'\n'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--xlsx',type=Path,default=DEFAULT_XLSX)
    p.add_argument('--write',action='store_true',help='Rebuild ONLY scheme_B files')
    args=p.parse_args()
    rec,gas=create(args.xlsx)
    expected=rows(rec,gas)
    # Windows CatMAP TableParser may open energies.txt using the GBK locale.
    # Keep parser-facing files ASCII-only; the source XLSX remains untouched.
    expected.encode('ascii')
    if args.write:
        (HERE/'energies.txt').write_text(expected,encoding='utf-8',newline='\n')
        (HERE/'energy_audit.csv').write_text(audit(rec),encoding='utf-8',newline='\n')
        print('Generated energies.txt and energy_audit.csv in scheme_B')
    else:
        current=(HERE/'energies.txt').read_text(encoding='utf-8-sig').replace('\r\n','\n')
        if current!=expected: raise ValueError('energies.txt differs from original Excel')
        print('MATCH energies.txt / original Excel')
        current_audit=(HERE/'energy_audit.csv').read_text(encoding='utf-8-sig').replace('\r\n','\n')
        if current_audit.replace('-0.000000000','0.000000000')!=audit(rec).replace('-0.000000000','0.000000000'): raise ValueError('energy_audit.csv differs from original Excel')
        print('MATCH energy_audit.csv / original Excel')
    print('15 surfaces; TS1 below FS = 0; max cycle error =',max(abs(e) for _,_,_,e in rec))
    for surf,v,_,_ in rec:
        flags=[f'{n}={v[ts]-v[state]:+.3f} eV' for n,ts,state in
          [('TS1-IS','TS1','S3'),('TS2-IS','TS2','S5'),('TS2-FS','TS2','S6')]
          if v[ts]-v[state]<-1e-7]
        if flags: print('FLAG',surf,', '.join(flags))

if __name__=='__main__': main()
