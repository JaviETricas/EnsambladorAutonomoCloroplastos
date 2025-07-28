EnsambladorAutonomoCloroplastos.

Esta herramienta lo que hace es juntar una serie de herramientas para realizar la tarea de obtencion del genoma 
del cloroplasto circular y alinearlo con un cloroplasto de referencia.
Para ello se usan las siguientes herramientas: reaper, Trimmomatic, novowrap, samtools, minimap2, y mafft. 

El programa a demas selecciona el fasta generado de novowrap con una mejor puntuacion con respecto al 
cloroplasto de referencia, descarta los ensamblajes que se consideran erroneos dejando una lista con 
estos para poder revisarlos manualmente y gracias a samtools se generara un tsv con los nucleotidos que pueden
generar error, y corrige los fallos de estos.

Por ultimo antes de alinear cambia el nombre al nombre de la especie + nombre del archivo + numero de bases del
coting.

Actualmente la herramienta esta en una version muy alpha por lo que aun no tengo un listado de los comandos que
se pueden usar, pero para poder hacerlo funcionar solo se requiere ejecutar el script/instaladordependencias.py 
y posteriormente ejecutar el cargadordearchivos.py y pasarle manualmente la direccion de los pareados en formato
fastq.gz o pasarle una carpeta con los pareados nombrados igual pero acabados en _1.fastq.gz o _2.fastq.gz para 
que el programa los detecte automaticamente y los procese.

## Índice
- [Descripción](#descripcion)
- [Instalación](#instalacion)
- [Uso general](#uso-general)
- [Casos de uso](#casos-de-uso)
- [Características avanzadas](#caracteristicas-avanzadas)
- [Limitaciones](#limitaciones)

## Descripción  <!-- id="descripcion" -->
…

## Instalación  <!-- id="instalacion" -->
…

## Uso general  <!-- id="uso-general" -->
…

## Casos de uso  <!-- id="casos-de-uso" -->
…

## Características avanzadas  <!-- id="caracteristicas-avanzadas" -->
…

## Limitaciones  <!-- id="limitaciones" -->
