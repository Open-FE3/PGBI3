# Coded on 2017.12.25

DIR=$(cd $(dirname $0);pwd)
CDIR=`echo $DIR | awk -F"/" '{print($NF)}'`
echo $CDIR
RDIR=tminami@eic.eri.u-tokyo.ac.jp:/home/tminami/tsunami/3D/ana_comp/src/solver/

#RDIR=minami@sakaki01.eri.u-tokyo.ac.jp:/home/minami/volcano/FEM_edge/FEM_edge_Nakadake/${CDIR}/
echo $RDIR
rsync -avz --delete -e ssh ./ $RDIR
#rsync -avz  -e ssh ./ $RDIR
