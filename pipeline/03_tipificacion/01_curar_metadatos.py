#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Curación de metadatos v2 - Proyecto Acinetobacter baumannii Sudamérica
Entrada:  metadatos/maestra_cruda.tsv
Salida:   metadatos/maestra_curada.tsv + reporte en pantalla
"""
import csv, re, unicodedata
from collections import Counter

ENTRADA = "maestra_cruda.tsv"
SALIDA = "maestra_curada.tsv"

# ---------------------------------------------------------------- utilidades
def sin_tildes(s):
    return ''.join(c for c in unicodedata.normalize('NFD', s)
                   if unicodedata.category(c) != 'Mn')

def titulo(s):
    menores = {'de', 'do', 'da', 'dos', 'das', 'del', 'la', 'las', 'los', 'e', 'y'}
    palabras = s.lower().split()
    out = []
    for i, w in enumerate(palabras):
        out.append(w if (i > 0 and w in menores) else w.capitalize())
    return ' '.join(out)

# ---------------------------------------------------------------- geografía
PAIS_STD = {
    'argentina': 'Argentina', 'bolivia': 'Bolivia', 'brazil': 'Brasil',
    'chile': 'Chile', 'colombia': 'Colombia', 'ecuador': 'Ecuador',
    'paraguay': 'Paraguay', 'peru': 'Peru', 'venezuela': 'Venezuela',
    'uruguay': 'Uruguay',
}

CIUDAD_FIX = {
    'buenos aeres': 'Buenos Aires',
    'buenos aires city': 'Buenos Aires',
    'ciudad autonoma buenos aires': 'Buenos Aires',
    'buenos aires, san andres': 'Buenos Aires',
    'rosario, santa fe': 'Rosario',
    'belem-pa': 'Belem',
    'uberlandia, mg': 'Uberlandia',
    'sao paulo, sao paulo': 'Sao Paulo',
    'pariquera-au': 'Pariquera-Acu',
    'mato grosso, pantanal': 'Pantanal',
    'cave furna do fim do morro do parafuso, paripiranga, bahia': 'Paripiranga',
}

REGIONES = {'Amapa', 'Roraima', 'Norte de Santander', 'Ayacucho', 'Amazonas',
            'Pantanal', 'Fernando de Noronha'}

INSTITUCION = re.compile(r'hospital|clinic|instituto|institute|universidad|'
                         r'university|laborator', re.I)

def parse_location(loc, pais_busqueda):
    pais = PAIS_STD.get(pais_busqueda.lower(), titulo(pais_busqueda))
    if not loc or loc.strip() == '':
        return pais, 'NA', 'NA', 'pais'

    partes = loc.split(':', 1)
    resto = partes[1].strip() if len(partes) > 1 else ''
    if resto == '':
        return pais, 'NA', 'NA', 'pais'

    sitio = 'NA'
    fragmentos = [f.strip() for f in resto.split(',')]
    lugares = []
    for f in fragmentos:
        if INSTITUCION.search(f):
            sitio = titulo(sin_tildes(f))
        else:
            lugares.append(f)
    lugar = ', '.join(lugares).strip()

    if lugar == '':
        return pais, 'NA', sitio, 'pais'

    clave = sin_tildes(lugar).lower().strip()
    if clave in CIUDAD_FIX:
        ciudad = CIUDAD_FIX[clave]
    else:
        lugar = re.sub(r'[-,]\s*[A-Z]{2}$', '', lugar)
        ciudad = titulo(sin_tildes(lugar))
        pp = [p.strip() for p in ciudad.split(',')]
        if len(pp) == 2 and pp[0].lower() == pp[1].lower():
            ciudad = pp[0]

    nivel = 'region' if ciudad in REGIONES else 'ciudad'
    return pais, ciudad, sitio, nivel

# ---------------------------------------------------------------- origen
AMBIENTAL = re.compile(
    r'cave|lake water|soil|petroleum|oil reservoir|root|bird|swab of a migratory|'
    r'limestone|pantanal|river|stream|^water$|water sample|'
    r'pharmaceutical industry', re.I)
AMB_HOSP = re.compile(
    r'high-touch|sink|surface|bed rail|equipment|ventilator circuit|'
    r'patient area|bathroom|keyboard|inanimate', re.I)

def clasificar_origen(iso_type, iso_source):
    s = (iso_source or '').strip()
    t = (iso_type or '').strip().lower()
    if AMB_HOSP.search(s):
        return 'ambiente_hospitalario'
    if AMBIENTAL.search(s):
        return 'ambiental'
    if t == 'clinical':
        return 'clinico'
    if t.startswith('environmental'):
        return 'ambiental'
    if s == '':
        return 'no_determinado'
    return 'clinico'

# ---------------------------------------------------------------- muestra
VACIOS = {'hospital', 'other', 'others', 'not collected', 'missing',
          'not applicable', 'unknown', 'na', 'nd', '-', 'not available',
          'other clinical sources', 'unspecified diagnostic sample',
          'clinical', 'clinical sample', 'human'}

def clasificar_muestra(src, origen):
    s = sin_tildes((src or '').strip().lower())
    if origen in ('ambiental', 'ambiente_hospitalario'):
        return 'no_aplica'
    if s == '' or s in VACIOS:
        return 'no_especificado'
    if re.search(r'blood|hemocultiv|sangre', s):
        return 'sangre'
    if re.search(r'cerebrospinal|csf|liquido cefalo|liquor', s):
        return 'lcr'
    if re.search(r'urine|urina|orina', s):
        return 'orina'
    if re.search(r'rectal|surveillance|screening|perianal|colonization|carrier', s):
        return 'vigilancia'
    if re.search(r'nasopharyng|throat|oropharyng|nasal', s):
        return 'respiratorio_superior'
    if re.search(r'trachea|endotracheal|bronch|sputum|lung|pulmonary|'
                 r'lower respiratory|inferior respiratory|\bbal\b', s):
        return 'respiratorio_inferior'
    if re.search(r'^respiratory|^aspirate|airway', s):
        return 'respiratorio_ne'
    if re.search(r'cath?eth?er|device|cateter|probe|\btube\b', s):
        return 'dispositivo'
    if re.search(r'peritone|pleural|ascit|synovial|bile|drain|abdominal fluid|'
                 r'joint|periprosthetic|lavage|punctate|bodily', s):
        return 'otro_esteril'
    if re.search(r'wound|tissue|ulcer|abscess|purulent|secretion|skin|'
                 r'soft tissue|burn|herida|bone|eschar|diabetic foot|'
                 r'exudate|tejido|biopsy|^swab|^hand$', s):
        return 'herida_tejido'
    return 'otro'

# ---------------------------------------------------------------- fecha
def parse_anio(f):
    m = re.match(r'\s*(\d{4})', f or '')
    if m:
        a = int(m.group(1))
        if 1950 <= a <= 2030:
            return str(a)
    return 'NA'

# ---------------------------------------------------------------- proceso
def main():
    with open(ENTRADA, newline='', encoding='utf-8') as fh:
        r = csv.DictReader(fh, delimiter='\t')
        filas = list(r)
        campos = r.fieldnames

    nuevos = ['pais', 'ciudad', 'sitio', 'nivel_geo', 'anio',
              'origen', 'tipo_muestra', 'tiene_ensamblado', 'tiene_reads']
    salida = campos + nuevos

    cont = {k: Counter() for k in
            ['pais', 'origen', 'tipo_muestra', 'nivel_geo', 'anio']}

    with open(SALIDA, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=salida, delimiter='\t')
        w.writeheader()
        for row in filas:
            pais, ciudad, sitio, nivel = parse_location(
                row.get('Location', ''), row.get('pais_busqueda', ''))
            row['pais'] = pais
            row['ciudad'] = ciudad
            row['sitio'] = sitio
            row['nivel_geo'] = nivel
            row['anio'] = parse_anio(row.get('Collection date', ''))
            row['origen'] = clasificar_origen(
                row.get('Isolation type', ''), row.get('Isolation source', ''))
            row['tipo_muestra'] = clasificar_muestra(
                row.get('Isolation source', ''), row['origen'])
            a = row.get('Assembly', '').strip()
            rn = row.get('Run', '').strip()
            row['tiene_ensamblado'] = 'si' if a not in ('', 'NULL') else 'no'
            row['tiene_reads'] = 'si' if rn not in ('', 'NULL') else 'no'
            for k in cont:
                cont[k][row[k]] += 1
            w.writerow(row)

    print(f"\nGenomas procesados: {len(filas)}")
    print(f"Salida: {SALIDA}\n")

    print("--- PAIS ---")
    for k, v in cont['pais'].most_common():
        print(f"  {k:<12} {v:5d}")

    print("\n--- ORIGEN ---")
    for k, v in cont['origen'].most_common():
        print(f"  {k:<24} {v:5d}")

    print("\n--- TIPO DE MUESTRA ---")
    for k, v in cont['tipo_muestra'].most_common():
        print(f"  {k:<24} {v:5d}")

    print("\n--- NIVEL GEOGRAFICO ---")
    for k, v in cont['nivel_geo'].most_common():
        print(f"  {k:<12} {v:5d}")

    sin_anio = cont['anio'].get('NA', 0)
    print(f"\n--- FECHA ---\n  con anio  {len(filas)-sin_anio:5d}"
          f"\n  sin anio  {sin_anio:5d}")

if __name__ == '__main__':
    main()
