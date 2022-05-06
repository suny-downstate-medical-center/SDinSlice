#!/usr/bin/bash 
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_standard/ --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/pad_standard.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/thick100/ --Lz=100 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/thick100.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/thick200/ --Lz=200 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/thick200.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/thick300/ --Lz=300 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/thick300.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/thick500/ --Lz=500 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/thick500.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/thick600/ --Lz=600 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/thick600.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/thick700/ --Lz=700 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/thick700.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/thick800/ --Lz=800 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/thick800.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/d45000/ --density=45000 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/d45000.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/d67500/ --density=67500 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/d67500.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/d112500/ --density=112500 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/d112500.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/d120000/ --density=120000 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/d120000.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/d45000_constBeta/ --constBeta=True --density=45000 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/d45000_constBeta.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/d67500_constBeta/ --constBeta=True --density=67500 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/d67500_constBeta.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/d112500_constBeta/ --constBeta=True --density=112500 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/d112500_constBeta.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/d120000_constBeta/ --constBeta=True --density=120000 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/d120000_constBeta.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/s2v02/ --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/s2v02.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/s2v1/ --sa2v=1.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/s2v1.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/s2v2/ --sa2v=2.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/s2v2.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/s2v4/ --sa2v=4.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/s2v4.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/s2v5/ --sa2v=5.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/s2v5.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/s2v6/ --sa2v=6.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/s2v6.json
python3 makeSbatch.py thick100
python3 makeSbatch.py thick200
python3 makeSbatch.py thick300
python3 makeSbatch.py thick500
python3 makeSbatch.py thick600
python3 makeSbatch.py thick700
python3 makeSbatch.py thick800
python3 makeSbatch.py d45000
python3 makeSbatch.py d67500
python3 makeSbatch.py d112500
python3 makeSbatch.py d120000
python3 makeSbatch.py d45000_constBeta
python3 makeSbatch.py d67500_constBeta
python3 makeSbatch.py d112500_constBeta
python3 makeSbatch.py d120000_constBeta
python3 makeSbatch.py s2v02
python3 makeSbatch.py s2v1
python3 makeSbatch.py s2v1
python3 makeSbatch.py s2v1
python3 makeSbatch.py s2v1
python3 makeSbatch.py s2v1
python3 makeSbatch.py pad_standard