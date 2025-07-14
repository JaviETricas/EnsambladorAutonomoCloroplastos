Ensamblador Autonomo Cloroplastos

Actualmente esta en una fase con errores!!!

El programa lo que realiza es la concatenacion de una serie de herramientas de manera autonoma,con el objetivo de ensamblar los cloroplastos de una serie de muestras pareadas elegidas. Para lograr este objetivo las muestras tienen que introducirse manualmente o estar en una carpeta la cual cada pareja termine en _1.fastq.gz, _2.fastq.gz para la autoseleccion.

Una vez el programa elige las multiples parejas, elimina duplicados y cebadores usando las erramienta tally y Trimmomatic ya portables con el programa, despues ensambla el cloroplasto con novowrap y genera una lista con las muestras que no se han ensamblado correctamente. Despues de esta seleccion se generara un bam y un informe con samtools y minimap2 los cuales nos permiten pasar de minimap2 → SAM y de samtools sort → BAM para posteriormente con la herramienta de samtools mpileup obtener un tsv el cual permitira en el futuro procesar las cadenas de DNA antes de la siguiente fase. 
Por ultimo se modifica el nombre que los programas han alterado, se posiciona todo el ADN en la misma linea, se comparan con una muestra elegida a travez de mafti y se genera un archivo que permita visualizar todo.

Actualmente el proceso funciona bien asta la seleccion de novowrap.

No dispone de manera de modificar los parametros de las erramientas actualmente sin entrar en el
codigo. 




