ldak=~/project-pasaniuc/software/ldak/ldak5.linux
bfile=./out/10k

$ldak --cut-weights ./out/sections --bfile $bfile --no-thin YES

$ldak --calc-weights-all ./out/sections --bfile $bfile

