#!/usr/bin/bash 
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/pad_standard/ --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/pad_standard.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/thick100/ --Lz=100 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/thick100.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/thick200/ --Lz=200 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/thick200.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/thick300/ --Lz=300 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/thick300.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/thick500/ --Lz=500 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/thick500.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/thick600/ --Lz=600 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/thick600.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/thick700/ --Lz=700 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/thick700.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/thick800/ --Lz=800 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/thick800.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/d45000/ --density=45000 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/d45000.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/d67500/ --density=67500 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/d67500.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/d112500/ --density=112500 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/d112500.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/d120000/ --density=120000 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/d120000.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/d45000_constBeta/ --constBeta=True --density=45000 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/d45000_constBeta.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/d67500_constBeta/ --constBeta=True --density=67500 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/d67500_constBeta.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/d112500_constBeta/ --constBeta=True --density=112500 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/d112500_constBeta.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/d120000_constBeta/ --constBeta=True --density=120000 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/d120000_constBeta.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/s2v02/ --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/s2v02.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/s2v1/ --sa2v=1.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/s2v1.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/s2v2/ --sa2v=2.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/s2v2.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/s2v4/ --sa2v=4.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/s2v4.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/s2v5/ --sa2v=5.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/s2v5.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/s2v6/ --sa2v=6.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/s2v6.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/betaNrn_165/ --alphaNrn=0.165 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/betaNrn_165.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/betaNrn_315/ --alphaNrn=0.315 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/betaNrn_315.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/betaNrn_39/ --alphaNrn=0.39 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/betaNrn_39.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/betaNrn_465/ --alphaNrn=0.465 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/betaNrn_465.json
python3 genCfgs.py --tstop=10000 --dir=/expanse/lustre/scratch/ckelley/temp_project/pad_data/betaNrn_54/ --alphaNrn=0.54 --sa2v=3.0 --r0=100 --k0=70 --ox=pad --nthreads=48 --nrec=48 --O2consume=True --pas=-70 --gpas=0.0001  --uniformRec=True --o2bath=0.01 --o2bc=0.01 cfgs/betaNrn_54.json
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
python3 makeSbatch.py s2v2
python3 makeSbatch.py s2v4
python3 makeSbatch.py s2v5
python3 makeSbatch.py s2v6
python3 makeSbatch.py pad_standard
python3 makeSbatch.py betaNrn_165
python3 makeSbatch.py betaNrn_315
python3 makeSbatch.py betaNrn_39
python3 makeSbatch.py betaNrn_465
python3 makeSbatch.py betaNrn_54
