## Dépendances :

Pour installer les packages nécessaire au fonctionnement, créezun environnement conda à partir du fichier projet_long.yml : *conda env create --file projet_long.yml*
Il est également nécessaird d'avoir installé **chromosight** disponible sur le dépôt Github : https://github.com/koszullab/chromosight

## Arguments : 

Le script à lancer est le script pipeline_final.sh. Si ce dernier n'est pas exécutable, utilisez chmod +x pipeline_final.sh.
Les arguments à entrer sont les suivants : 
* Un fichier .txt (ex : hic.txt) contenant un SRR (numéro d'accession Short Read Archive) d'une expérience de HIC, par ligne
* Un fichier .txt (ex: chip_seq.txt) contenant un SRR d'une expérience de ChIP-Seq (avec immunoprécipitation), par ligne
* Un fichier .txt (ex: input.txt) contenant un SRR d'une expérience de ChIP-Seq (sans immunoprécipitation), par ligne
* Un répértoire où l'arborescence des analyses sera construite
* Un répértoire où se situe les fichiers fasta des chromosomes de levures

### Remarques :
* **Attention** pour les deux derniers agruments, à la fin du chemin absolu, ou relatif ne pas mettre de '/'
* **Attention** pour le dernier argument, les fichiers dans le répértoire doivent avoir l'extension .fasta, et non pas .fa ou .faa 

## Exemple de commande :

./pipeline_final.sh hic.txt chipseq.txt input.txt /home/anthony data

## Sortie :

* HIC/results/SRRXXXXXXX/ : cartes de contacts pour chaque chromosome
* HIC/results/SRRXXXXXXX/ : barplot d'enrichissement final

