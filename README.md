## EnsambladorAutonomoCloroplastos.
[![PyPI version](https://img.shields.io/pypi/v/numpy)](https://pypi.org/project/numpy)

## Índice

- [Descripción](#descripcion)
- [Instalación](#instalacion)
- [Uso general](#uso-general)
- [Casos de uso](#casos-de-uso)
- [Características avanzadas](#caracteristicas-avanzadas)
- [Limitaciones](#limitaciones)

## Descripcion

Esta herramienta lo que hace es juntar una serie de herramientas para realizar la tarea de obtencion del genoma 
del cloroplasto circular y alinearlo con un cloroplasto de referencia.

Para ello se usan las siguientes herramientas: reaper, Trimmomatic, novowrap, samtools, minimap2, y mafft. 

El programa a demas selecciona el fasta generado de novowrap con una mejor puntuacion con respecto al 
cloroplasto de referencia, descarta los ensamblajes que se consideran erroneos dejando una lista con 
estos para poder revisarlos manualmente y gracias a samtools se generara un tsv con los nucleotidos que pueden
generar error, y corrige los fallos de estos.

Por ultimo antes de alinear cambia el nombre al nombre de la especie + nombre del archivo + numero de bases del
coting.

## Instalación  <!-- id="instalacion" -->

La herramienta actualmente no se encuentra fuera de github, por lo que a la hora de instalar solo se puede 
descargando este repositorio.

## Uso general

Para ejecutar el programa de manera general se tiene que ejecutar el script cargadordearchivos.py este se encarga 
de instalar las dependencias necesarias.

Porteriormente te pide que le pases una pareja de archivos .fastq.gz de manera manual o que le pases un directorio,
si usas esta segunda opcion los archivos tienen que acabar en _1.fastq.gz _2.fasq.gz para que los localize 
automaticamente y los empareje

## Casos de uso  <!-- id="casos-de-uso" -->

El programa funciona actualmente solo con fastq.gz pareados, en linux y para el ensamblaje de cloroplastos. Los
comandos de las herramientas originales estan limitados, asi que el principal caso de uso es el de tener que hacer
un ensamblaje de cloroplastos con un gran volumen de archivos en el formato indicado. De esta manera puedes dejar
correr el programa durante horas o dias asta tener los archivos resultantes correctos, y una lista con los archivos
fallidos para usar herramientas mas especificas.

## Características avanzadas  <!-- id="caracteristicas-avanzadas" -->

En la actualidad usa un cloroplasto de hurdeum vulgare como referencia para reconstruir el cloroplasto. Para cambiar
esto debes descargar otro cloroplasto de referencia y guardarlo en la carpeta cloroplasto_referencia en formato gb

# Argumentos de scripts.
  **instaladordependecias**
  
    ``` 
    # Reinstalar todas las dependencias.
    ./instaladordependencias.py --force



## Limitaciones  <!-- id="limitaciones" -->

No dispone de una herramienta de argumentos completa, y para usar los argumentos tienes que ejecutar los diferentes
scripts por separado. Tampoco dispone de portabilidad total a windows.

Si una herramienta falla con una pareja de datos, no existen caminos secundarios para corregir el problema.

Otra limitacion es a la hora de nombrar los fasta, aun no esta añadido un argumento para cambiar el nombre de la 
especie por otra diferente.

