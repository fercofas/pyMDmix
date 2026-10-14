# Construcción de cajas de solvente mixto con agua OPC, paso a paso

Este documento describe cómo se construye, valida y registra en pyMDMix una caja de solvente mixto (cosolvente orgánico en agua OPC). Aplica a los cosolventes de la biblioteca (ETA, ANT, MAM, ISO, ISO5, MOH, PYR1 e ION) y a cualquier cosolvente nuevo cuyos parámetros ya estén en un archivo `.off`. Los archivos están en `examples/solvent_boxes/`.

> Estado de validación. La caja de etanol (ETAO) se construyó y validó en una instalación real: densidad de 0.9678 ± 0.0052 g/cm³ (experimental aproximado para 17 % p/p y 300 K: 0.97 g/cm³), sin deriva entre mitades de la trayectoria y 1672 sitios virtuales EPW. Los perfiles de los demás cosolventes siguen el mismo procedimiento y se están validando uno a uno. Cada caja debe pasar la validación de la sección 5 antes de usarse.

## 1. Fundamento

1. **Los cosolventes no se reparametrizan.** Las cargas RESP y los tipos atómicos del cosolvente se derivaron en fase gas o con GAFF y no dependen del modelo de agua. Lo que cambia es el equilibrio solvente-agua, por eso la caja se reconstruye, se reequilibra y se valida (densidad y agregación).
2. **OPC es un modelo de cuatro puntos.** La molécula de agua incluye un sitio virtual `EPW`. Por eso hay que cargar `leaprc.water.opc`, escribir los PDB con `include_ep` y comprobar que el número de `EPW` es igual al número de aguas.
3. **ff19SB y OPC se diseñaron en conjunto.** El script carga `leaprc.protein.ff19SB` y `leaprc.gaff` (los cosolventes orgánicos usan tipos GAFF).
4. **Tamaño de caja.** Las cajas originales de MDMix miden unos 20 Å. Con un corte de 9 Å, `pmemd` exige que el radio de la esfera máxima de la celda sea mayor que corte + 2 Å, y `pmemd.cuda` rechaza cajas pequeñas ("Small box detected"). El script replica la caja original `REP` veces por lado (por defecto, el menor entero que da al menos 36 Å). La fracción de mezcla no cambia.

## 2. Requisitos

```bash
conda activate pymdmix3
source $AMBERHOME/amber.sh
which tleap cpptraj packmol pmemd pmemd.cuda
ls $AMBERHOME/dat/leap/cmd | grep -E "^leaprc\.(protein.ff19SB|gaff|water.opc)$"
```

## 3. Decidir la composición con `boxcalc.py`

`boxcalc.py` reemplaza a la hoja `molecule_calculator.ods`. Calcula el número de moléculas a partir de volúmenes de componentes puros (V = M / (ρ N_A)). No considera el volumen de exceso de la mezcla, por eso la densidad final se obtiene del NPT.

```bash
cd examples/solvent_boxes
python boxcalc.py table                                   # cosolventes disponibles
python boxcalc.py library                                 # concentración real de las cajas TIP3P incluidas
python boxcalc.py mix ETA --vv 20 --side 20               # moléculas para 20 % v/v en 20 Å
python boxcalc.py check ETA --nco 136 --nwat 1672 --side 39.67
python boxcalc.py pure ETA --side 23                      # caja de cosolvente puro
python boxcalc.py molar UREA --molarity 7 --volume 0.1    # solución de molaridad dada
```

Concentración real de las cajas TIP3P de la biblioteca, calculada con `boxcalc.py library`:

| Caja | Cosolvente + agua | % v/v | % p/p | Fracción molar | Molaridad |
|---|---|---|---|---|---|
| ETAWAT20 | 17 + 209 | 20.9 | 17.2 | 0.0752 | 3.5 M |
| ANTWAT20 | 17 + 209 | 19.1 | 15.6 | 0.0752 | 3.6 M |
| MAMWAT20 | 18 + 209 | 19.6 | 22.0 | 0.0793 | 3.8 M |
| ISOWAT20 | 14 + 209 | 22.1 | 18.3 | 0.0628 | 2.9 M |
| ISOWAT5 | 4 + 310 | 5.2 | 4.1 | 0.0127 | 0.7 M |
| MOHWAT20 | 56 + 470 | 21.1 | 17.5 | 0.1065 | 5.2 M |
| PYRWAT1 | 1 + 438 | 1.0 | 1.0 | 0.0023 | 0.08 M (ver nota) |

Observaciones sobre la hoja de cálculo original:

- En la hoja `ethanol_20%`, el volumen molecular del etanol (`ETHOH_vol`, celda B16) se calcula con la masa molar de la urea (60.06 g/mol, celda B6) y no con la del etanol (46.07 g/mol). El resultado correcto es 97.0 Å³ por molécula y no 126.4 Å³. Con el valor correcto, 17 ETA + 209 WAT corresponde a 20.9 % v/v, consistente con el nombre "20".
- Las densidades de ANT, MAM, ISO, MOH y PYR que `boxcalc.py` usa no estaban en la hoja. Son valores de literatura (marcados `lit` en `boxcalc.py table`); confírmalos con tu proveedor antes de publicar.
- La caja PYRWAT1 guardada en la biblioteca tiene 438 aguas en unos 28 Å de lado (densidad de partida cercana a 0.60 g/mL), es decir, es una caja con mucho vacío. `build_box.sh` no usa ese lado: para PYR1 parte de 23.65 Å, que corresponde a la composición real.

## 4. Construir la caja

Una carpeta de trabajo vacía por solvente:

```bash
mkdir -p ~/cajas/ANT && cd ~/cajas/ANT
nohup bash /ruta/a/pyMDmix/examples/solvent_boxes/build_box.sh ANT > build.log 2>&1 &
tail -f build.log
```

Perfiles incluidos: `ETA`, `ANT`, `MAM`, `ISO`, `ISO5`, `MOH`, `PYR1`, `ION`. Cada perfil toma el `.off` TIP3P de `pyMDMix/data/solventlib/`, el número de moléculas y el lado de la caja original.

Etapas:

| Etapa | Qué hace | Salida principal |
|---|---|---|
| 1 | Extrae las unidades del cosolvente del `.off` TIP3P | `RES.pdb`, `RES.lib` |
| 2 | Obtiene una molécula de agua OPC y comprueba el átomo `EPW` | `wat_opc.pdb` |
| 3 | Empaqueta con packmol y construye la topología con tleap | `mix.parm7`, `mix.rst7` |
| 4 | Minimiza y calienta (50 ps) | `heat.rst7` |
| 5 | NPT de 10 ns en GPU | `npt.nc`, `npt.out` |
| 6 | Calcula densidad y tamaño de caja; elige el cuadro representativo | `6_analysis.log`, `best_frame.txt` |
| 7 | Escribe `best.pdb` (con `include_ep`) y el `.off` final | `UNIT.off` |
| 8 | Genera el `.config` del solvente y los RDF de agregación | `SOLV.config`, `rdf_*.dat` |

Para reanudar desde una etapa: `bash build_box.sh ANT 5`. Variables que se pueden sobrescribir desde el entorno: `REP`, `NPT_NS`, `DISCARD_NS`, `TEMP`, `FFS`, `HEAT_CMD`, `NPT_CMD`, `SRC_OFF`, `SOLV`, `UNIT`.

Ejemplos:

```bash
REP=3 bash build_box.sh MOH            # caja más grande
NPT_NS=20 bash build_box.sh ISO5       # NPT más largo
REP=1 HEAT_CMD=pmemd NPT_CMD=pmemd bash build_box.sh PYR1   # CPU, solo si la caja ya es grande
```

### 4.1 Cosolvente nuevo

Si ya tienes un `.off` con la unidad del cosolvente parametrizada, define el perfil por variables:

```bash
python boxcalc.py mix TFE --vv 20 --side 21 --rep 2    # da NCO y NWAT
COSOLV=TFE SRC_OFF=/ruta/TFE.off COMP="TFE:NN" NW0=NNN SIDE_EQ=21 \
SOLV=TFEO UNIT=TFEOPC20 CFG_SRC=/ruta/tfe.config RDFA=":TFE@O1" RDFB=":TFE@O1" \
bash build_box.sh TFE
```

Aquí `COMP` y `NW0` describen la caja original sin replicar (el script aplica `REP`). El `.off` debe contener una unidad con el nombre de cada residuo de `COMP`. Si el cosolvente no está parametrizado, primero hay que derivar cargas y parámetros (antechamber, RESP), cosa que este procedimiento no hace.

## 5. Validar la caja

Condiciones mínimas antes de registrar:

1. `6_analysis.log`: la deriva de densidad entre la primera y la segunda mitad debe ser cercana a cero (|deriva| menor que la desviación estándar, típicamente < 0.002 g/cm³), y la densidad media debe estar dentro de pocos por ciento del valor experimental de la mezcla a esa fracción en peso y temperatura.
2. La etapa 7 debe terminar con "sitios virtuales EPW en best.pdb" igual al valor esperado, y el `.off` debe listar las unidades de los cosolventes, `WAT` y la unidad de la caja.
3. La etapa 8 imprime el máximo del primer pico y el valor de cola de g(r) entre 9 y 12 Å (debe ser cercano a 1). Un g(r) cosolvente-cosolvente con un primer máximo alto (> 3) a menos de 5 Å indica agregación. Compara con la caja TIP3P si dudas. Las cajas con pocas moléculas de cosolvente (ISO5, PYR1) dan RDF ruidosos.
4. Revisa visualmente la última trayectoria (sin moléculas atravesando la caja de forma anómala).

Para obtener la densidad experimental de referencia de una mezcla etanol-agua, por ejemplo, usa las tablas de densidad por % p/p y temperatura (ver la columna "% p/p" de `boxcalc.py check`).

## 6. Registrar el solvente en pyMDMix

```bash
mkdir ~/mdmix_cajas && cd ~/mdmix_cajas
cp ~/cajas/ANT/ANTOPC20.off ~/cajas/ANT/ANTO.config .
grep -n "name\|objectfile\|boxunit\|watermodel" ANTO.config    # comprobar
mdmix create solvents -f ANTO.config
mdmix info solvents | grep ANTO
```

El `.config` generado por la etapa 8 reutiliza las secciones `[PROBES]` y `[TYPES]` del `.config` original. Si el comando falla con `InvalidPath`, el mensaje indica qué archivo no se encontró (normalmente `objectfile`).

Nombres de solvente que usa esta construcción: ETAO, ANTO, MAMO, ISOO, ISO5O, MOHO, PYR1O e IONO. Mantienen separadas las cajas TIP3P (ETA, ANT, ...) de las OPC.

## 7. Probar la caja con un sistema

```ini
[SYSTEM]
NAME = PEP
OFF = pep.off
# Con agua OPC es obligatorio cargar su campo de fuerza
EXTRAFF = leaprc.water.opc

[MDSETTINGS]
SOLVENTS = ETAO
NREPL = 1
NANOS = 1
RESTR = FREE
```

```bash
mdmix create project -n PEP_OPC -f pep_opc.cfg
```

En la réplica creada, `tleap` debe cargar `leaprc.water.opc` y la topología debe conservar los sitios `EPW`.

## 8. Problemas frecuentes

| Síntoma | Causa y solución |
|---|---|
| `max pairlist cutoff must be less than unit cell max sphere radius` o `Small box detected` | Caja demasiado pequeña. Sube `REP`. |
| tleap no encuentra un tipo atómico o parámetro del cosolvente | Falta cargar el campo de fuerza del tipo (por ejemplo `leaprc.gaff`) o un `frcmod`. Ajusta `FFS`. |
| `saveOff` rechaza la lista de unidades | Revisa `7_off.log`. La sintaxis `saveOff { ... }` se verificó con la caja de etanol (ETAOPC20.off contiene ETA, WAT y la unidad de la caja). |
| El número de `EPW` no coincide | El PDB se escribió sin `include_ep`, o el residuo de agua no es `WAT`. |
| Aviso de `tleap` con "Errors = N" distinto de 0 | Revisa `3_build.log` antes de seguir. |
| La densidad no converge | Alarga el NPT con `NPT_NS=20` y aumenta `DISCARD_NS`. |
