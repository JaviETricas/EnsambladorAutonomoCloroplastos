Ensamblador Autonomo Cloroplastos

Actualmente esta en una fase con errores!!!

El programa lo que realiza es la concatenacion de una serie de herramientas de manera autonoma,
con el objetivo de ensamblar los cloroplastos de una serie de muestras pareadas elegidas. Para 
lograr este objetivo las muestras tienen que introducirse manualmente o estar en una carpeta la
cual cada pareja termine en _1.fastq.gz, _2.fastq.gz para la autoseleccion.
Una vez el programa elige las multiples parejas, elimina duplicados y cebadores usando las erra-
mienta tally y Trimmomatic ya portables con el programa, despues ensambla el cloroplasto con no-
vowrap y genera una lista con las muestras que no se han ensamblado correctamente. Las que si 
que lo han logrado se procesan y se comparan con una muestra elegida a travez de mafti.

Actualmente el proceso funciona bien asta la seleccion de novowrap.




