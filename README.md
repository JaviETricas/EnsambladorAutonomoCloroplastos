## EnsambladorAutonomoCloroplastos.
[![PyPI version](https://img.shields.io/pypi/v/numpy)](https://pypi.org/project/numpy) 


Para una comprobacion rapida, use el argumento --test y descargara una carpeta llamada test con unos archivos 
fastq.gz pareados que podra usar como ejemplo.

El proceso total puede costar varias horas.


## Índice

- [Descripción](#descripcion)

- [Requisitos](#requisitos)
- [Instalación](#instalacion)
- [Uso general](#uso-general)
- [Inicio rapido](#inicio-rapido)
- [Casos de uso](#casos-de-uso)
- [Características avanzadas](#caracteristicas-avanzadas)
  - [Argumentos de scripts](#argumentos-de-scripts)
- [Limitaciones](#limitaciones)
- [Bibliografia](#bibliografia)

## Descripcion

Esta herramienta lo que hace es juntar una serie de herramientas para realizar la tarea de obtencion del genoma 
del cloroplasto circular y alinearlo con un cloroplasto de referencia.

Para ello se usan las siguientes herramientas: reaper, Trimmomatic, novowrap, samtools, minimap2, y mafft. 
El programa a demas selecciona el fasta generado de novowrap con una mejor puntuacion con respecto al 
cloroplasto de referencia, descarta los ensamblajes que se consideran erroneos dejando una lista con 
estos para poder revisarlos manualmente y gracias a samtools se generara un tsv con los nucleotidos que pueden
generar error, y corrige los fallos de estos.

Por ultimo antes de alinear cambia el nombre al nombre de la especie + nombre del archivo + numero de bases del
coting, y corrige los errores basicos de nucleotidos por el proceso de mayoria, si alguna base esta en duda, se 
selecciona la mayoria.

Una vez Alineado todo, el programa revisara en secciones de 400 en 400pb que el k-mers sea correcto y si este 
baja del umbral del 60% y la seccion es superior a 2000pb probara a invertir la seccion y volver a ensamblar para
posteriormente elegir la cadena con un resultado mejor y guardar la otra en temporaldocs.

<img width="1175" height="316" alt="Captura desde 2025-08-28 12-57-57" src="https://github.com/user-attachments/assets/7712b317-53b5-4d40-9a6d-035f78696869" />

<img width="1114" height="316" alt="Captura desde 2025-08-28 12-58-57" src="https://github.com/user-attachments/assets/59a10905-e9d6-4927-b36b-2890be0f9570" />

Esto se realiza por que las secciones IRa, IRb, y SSC se suelen voltear, a demas algunas veces el ensamblador voltea
estas secciones sin querer y asi evitamos los errores que nos podria causar esto.

## Requisitos

Para que el programa funcione correctamente en la actualidad se necesita:

 - Sistema operativo linux.
 - Conda instalado
 - Bioconda
 - Comandos de Github

## Instalacion

La herramienta actualmente no se encuentra fuera de github, por lo que a la hora de instalar solo se puede 
descargando este repositorio.

## Uso general

Para ejecutar el programa de manera general se tiene que ejecutar el script cargadordearchivos.py este se encarga 
de instalar las dependencias necesarias.

Porteriormente te pide que le pases una pareja de archivos .fastq.gz de manera manual o que le pases un directorio,
si usas esta segunda opcion los archivos tienen que acabar en _1.fastq.gz _2.fasq.gz para que los localize 
automaticamente y los empareje

## Casos de uso

El programa funciona actualmente solo con fastq.gz pareados, en linux y para el ensamblaje de cloroplastos. Los
comandos de las herramientas originales estan limitados, asi que el principal caso de uso es el de tener que hacer
un ensamblaje de cloroplastos con un gran volumen de archivos en el formato indicado. De esta manera puedes dejar
correr el programa durante horas o dias asta tener los archivos resultantes correctos, y una lista con los archivos
fallidos para usar herramientas mas especificas.

## Inicio rapido

Para comprobar que todo funciona correctamente sigue los siguientes pasos.

    ``` 
    # Si no tienes archivos fastq.gz para la prueba usa
    ./cargadordearchivos.py --test

    ``` 
    # Si ya tienes su archivo usa
    ./cargadordearchivos.py

Despues le aparecera una ventana como esta:

    ```   
    ¿Quieres introducir manualmente dos archivos .fastq.gz? [y/n]:
  
Puede introducir manualmente la direccion de cada fastq.gz o darle a n, en ese caso vera que aparece lo siguiente:

    ``` 
    Ruta de la carpeta con .fastq.gz:
  
Si le pasas la ruta a una carpeta busca manualmente los archivos fastq.gz que acaben en 1 y 2 y se llamen igual,
de esa forma los empareja y ya no tienes que hacer nada mas asta obtener los resultados. Si no detecta todos, pasalos
manualmente.
  

## Caracteristicas avanzadas

En la actualidad usa un cloroplasto de hurdeum vulgare como referencia para reconstruir el cloroplasto. Para cambiar
esto debes descargar otro cloroplasto de referencia y guardarlo en la carpeta cloroplasto_referencia en formato gb

### Argumentos de scripts

  **cargadordearchivos**

    ``` 
    # Nombre cientifico de la especie para el fasta
    ./cargadordearchivos.py --species "name sp"


    ``` 
    # Borra los bam al ser archivos pesados
    ./cargadordearchivos.py --dellbam

    ``` 
    # vacia la papelera de los documentos borrados
    ./cargadordearchivos.py --dell

    ``` 
    # Fuerza la reinstalacion de las librerias
    ./cargadordearchivos.py --force


  **instaladordependecias**
  
    ``` 
    # Reinstalar todas las dependencias.
    ./instaladordependencias.py --force


  **SeleccionNovowrap**
  
    ```
    # Directorio raiz
    ./SeleccionNovowrap.py -r
    # Documento de errores
    ./SeleccionNovowrap.py -o

    
## Limitaciones

No dispone de una herramienta de argumentos completa, y para usar los argumentos tienes que ejecutar los diferentes
scripts por separado. Tampoco dispone de portabilidad total a windows.

Si una herramienta falla con una pareja de datos, no existen caminos secundarios para corregir el problema.

Otra limitacion es a la hora de nombrar los fasta, aun no esta añadido un argumento para cambiar el nombre de la 
especie por otra diferente.

El programa busca inversiones y las alinea con el consenso, y te guarda la no invertida para la filogenia, pero no
detecta a que genes pertenece ni que seccion o si es un fallo del ensamblador o una mutacion del cloroplasto.

## Bibliografia

Sancho, R., Cantalapiedra, C. P., López‐Alvarez, D., Gordon, S. P., Vogel, J. P., Catalán, P., & Contreras‐Moreira, B. (2017). Comparative plastome genomics and phylogenomics of Brachypodium: flowering time signatures, introgression and recombination in recently diverged ecotypes. New Phytologist, 218(4), 1631-1644. https://doi.org/10.1111/nph.14926

Gong, C., Huang, Y., Liu, M., Zhou, Y., Xu, Y., Mohammed, N., Qiao, X., Zuccolo, A., Xie, W., Wing, R. A., Zhang, J., Zhou, F., & Lin, Y. (2025). Continuous infiltration and evolutionary trajectory of nuclear organelle DNA inOryza. Genome Research. https://doi.org/10.1101/gr.279609.124

Wu P, Xu C, Chen H, Yang J, Zhang X, Zhou S. NOVOWrap: An automated solution for plastid genome assembly and structure standardization. Mol Ecol Resour. 2021 Aug;21(6):2177-2186. doi: 10.1111/1755-0998.13410. Epub 2021 May 25. PMID: 33934526.

Greiner S, Lehwark P, Bock R. OrganellarGenomeDRAW (OGDRAW) version 1.3.1: expanded toolkit for the graphical visualization of organellar genomes. Nucleic Acids Res. 2019 Jul 2;47(W1):W59-W64. doi: 10.1093/nar/gkz238. PMID: 30949694; PMCID: PMC6602502.


