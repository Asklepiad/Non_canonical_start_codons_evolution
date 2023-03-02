# !/bin/bash

proj="C_psittaci"
identity=75
proteinortho6.pl --project=${proj} ./orto_rows/*.fasta --cpus=7 --debug --selfblast --singles --identity=${identity}
