wget https://github.com/alecjacobson/common-3d-test-models/raw/master/data/stanford-bunny.obj
bin/ACVD stanford-bunny.obj 3000 0
bin/ACVD stanford-bunny.obj 3000 1.5
wget https://github.com/alecjacobson/common-3d-test-models/raw/master/data/fandisk.obj
bin/ACVDQ fandisk.obj 3000 0
wget https://github.com/alecjacobson/common-3d-test-models/raw/master/data/horse.obj
bin/AnisotropicRemeshingQ horse.obj 1000 1.5
wget http://graphics.stanford.edu/data/3Dscanrep/xyzrgb/xyzrgb_statuette.ply.gz
gunzip xyzrgb_statuette.ply.gz
bin/ACVDQ xyzrgb_statuette.ply 100000 1.5
bin/ACVDQP xyzrgb_statuette.ply 100000 1.5
bin/ACVDQP xyzrgb_statuette.ply 100000 1.5 -np 3
