#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figura 7: los tres mecanismos de ausencia aparente."""
import os
B = os.path.expanduser("~/abaumannii")
W,H = 940, 700
def t(x,y,txt,sz=10,w="400",anc="start",col="#222",fam="Helvetica,Arial"):
    return (f'<text x="{x}" y="{y}" font-family="{fam}" font-size="{sz}" font-weight="{w}" '
            f'text-anchor="{anc}" fill="{col}">{txt}</text>')
s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
   f'<rect width="{W}" height="{H}" fill="white"/>']
X0=50; PW=840
BLOQ=[
 ("A","Fragmentacion del ensamblado","#8c6bb1",
  "El gen y el elemento con el que se quiere evaluar su contexto quedan en contigs distintos.",
  "4,5 % de evaluabilidad", "944 de 988 copias de bla-OXA-51-like no evaluables para el contexto de ISAba1"),
 ("B","Linealizacion de elementos circulares","#d95f02",
  "El punto de corte del circulo cae dentro del gen, que queda partido entre los dos extremos del contig.",
  "+9 portadores (14 %)", "64 &#8594; 73 portadores de bla-OXA-72 tras el rastreo dirigido; ningun falso positivo"),
 ("C","Divergencia de secuencia","#1b9e77",
  "El gen esta presente pero su identidad queda por debajo del umbral de notificacion de la herramienta.",
  "5 de 900 genomas", "3 genomas ST2742 con delecion en marco de 102 pb que elimina el motivo STFK"),
]
y=48
for letra,tit,c,desc,cifra,pie in BLOQ:
    s.append(f'<rect x="{X0}" y="{y}" width="{PW}" height="185" fill="#fbfbfb" stroke="#ddd" rx="4"/>')
    s.append(f'<rect x="{X0}" y="{y}" width="5" height="185" fill="{c}"/>')
    s.append(t(X0+18,y+24,letra,15,"700",col=c))
    s.append(t(X0+42,y+24,tit,13,"700"))
    s.append(t(X0+42,y+42,desc,9.5,col="#555"))
    gy=y+62; gx=X0+42; GW=520
    if letra=="A":
        s.append(f'<rect x="{gx}" y="{gy}" width="200" height="20" fill="#e8e0f0" stroke="#8c6bb1"/>')
        s.append(f'<rect x="{gx+52}" y="{gy+3}" width="42" height="14" fill="#8c6bb1" rx="2"/>')
        s.append(t(gx+73,gy+13,"gen",8,"600",anc="middle",col="white"))
        s.append(t(gx+100,gy+34,"contig 1",8.5,anc="middle",col="#666"))
        s.append(f'<rect x="{gx+240}" y="{gy}" width="200" height="20" fill="#e8e0f0" stroke="#8c6bb1"/>')
        s.append(f'<rect x="{gx+330}" y="{gy+3}" width="52" height="14" fill="#666" rx="2"/>')
        s.append(t(gx+356,gy+13,"ISAba1",7.5,"600",anc="middle",col="white"))
        s.append(t(gx+340,gy+34,"contig 2",8.5,anc="middle",col="#666"))
        s.append(t(gx+218,gy+14,"?",14,"700",anc="middle",col="#b03a20"))
        s.append(t(gx,gy+58,"La distancia entre ambos no es determinable.",9,col="#555"))
    elif letra=="B":
        s.append(f'<circle cx="{gx+70}" cy="{gy+22}" r="30" fill="none" stroke="#d95f02" stroke-width="2.5"/>')
        s.append(f'<path d="M {gx+70} {gy-8} A 30 30 0 0 1 {gx+96} {gy+37}" fill="none" stroke="#d95f02" stroke-width="9"/>')
        s.append(t(gx+70,gy+56,"elemento circular",8.5,anc="middle",col="#666"))
        s.append(t(gx+150,gy+26,"&#8594;",16,anc="middle",col="#999"))
        s.append(f'<rect x="{gx+180}" y="{gy+12}" width="300" height="20" fill="#fbe7d8" stroke="#d95f02"/>')
        s.append(f'<rect x="{gx+180}" y="{gy+15}" width="46" height="14" fill="#d95f02"/>')
        s.append(f'<rect x="{gx+434}" y="{gy+15}" width="46" height="14" fill="#d95f02"/>')
        s.append(t(gx+330,gy+54,"el gen queda partido en los dos extremos",9,anc="middle",col="#555"))
    else:
        s.append(f'<rect x="{gx}" y="{gy+10}" width="440" height="20" fill="#e0efe8" stroke="#1b9e77"/>')
        s.append(f'<rect x="{gx+150}" y="{gy+10}" width="70" height="20" fill="#b03a20" opacity=".18"/>')
        s.append(f'<line x1="{gx+150}" y1="{gy+6}" x2="{gx+220}" y2="{gy+6}" stroke="#b03a20" stroke-width="2"/>')
        s.append(t(gx+185,gy+2,"102 pb ausentes",8.5,"600",anc="middle",col="#b03a20"))
        s.append(t(gx+185,gy+24,"STFK",8.5,"700",anc="middle",col="#b03a20"))
        s.append(t(gx+220,gy+52,"240 de 274 residuos identicos a la referencia; sitio activo eliminado",9,anc="middle",col="#555"))
    s.append(f'<rect x="{X0+600}" y="{y+62}" width="220" height="52" fill="{c}" opacity=".10" rx="3"/>')
    s.append(t(X0+710,y+94,cifra,15,"700",anc="middle",col=c))
    s.append(t(X0+42,y+168,pie,9,col="#444"))
    y+=200
s.append(t(X0,y+4,"Los tres mecanismos producen ausencias aparentes indistinguibles de ausencias reales "
                  "sin verificacion dirigida.",10,"600",col="#333"))
s.append('</svg>')
open(f"{B}/resultados/figuras/Fig7_mecanismos.svg","w").write("\n".join(s))
print("escrita: resultados/figuras/Fig7_mecanismos.svg")
