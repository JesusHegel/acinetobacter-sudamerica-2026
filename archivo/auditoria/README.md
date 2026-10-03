# Auditoría

Comprobaciones realizadas sobre los resultados a medida que se obtenían. No
producen ninguna cifra del manuscrito: verifican que las que sí lo hacen son
consistentes entre sí.

| Script | Qué comprueba |
|---|---|
| `auditoria_aritmetica.py` | Que los totales y porcentajes del texto cuadran con las tablas |
| `auditoria_bioprojects.py` | Proyectos que abarcan más de un país, e independencia de las ausencias replicadas |
| `auditoria_42.py` | Inventario de la sección de resultados: tipos capsulares, replicones y eventos |
| `E1_verificaciones.py` | Cuatro comprobaciones puntuales sobre el conjunto |
| `B2_B3.py` | Replicón del segundo aislado ecuatoriano ST108 y causas de los genomas sin ST asignado |
| `cierre_49_50.py` | Distribución de r3-T18 y regla de asignación de genes a eventos |
| `diag_contig_reps.py` | Distribución de replicones por contig y diámetro interno de los agrupamientos |
| `diag_rep_contiguidad.py` | Si la detección de genes *rep* depende de la contigüidad del ensamblado |
| `resolver_contradiccion.py` | Dos contigs de 7851 pb, uno con blaOXA-72 y otro sin él |
| `oxa51_delecion.py`, `oxa51_union.py` | Caracterización del gen intrínseco en los aislados ST2742 |

Los dos últimos documentan una variante de `blaOXA-51-like` con una deleción
que el catálogo de referencia no notifica a nivel de variante. El hallazgo no
entró al manuscrito pero explica parte de las decisiones de recuento descritas
en métodos.
