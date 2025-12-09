import subprocess
import random
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil

# === PARAMÈTRES ===
natoms = 64 # Nombre de particules
parasites = 9        # Lignes parasites par timestep
dumpfreq = 2500      # Fréquence des dumps
time = 18000000
Nn = 8
numrep = 10
ncore = 2

def filtrer_fichier(input_filepath, natoms, parasites, output_filepath):
    lignesn = []
    with open(input_filepath, 'r') as f:
        lignes = f.readlines()

    for i in range(len(lignes)):
        if i % (natoms + parasites) > parasites - 1:
            lignesn.append(lignes[i])
    dumpsteps = len(lignesn) // natoms

    with open(output_filepath, "w") as fa:
        fa.writelines(lignesn)

    return output_filepath, dumpsteps



def fichier_to_matrice(natoms, dumpsteps, input_filepath):
    datamat = np.loadtxt(input_filepath)
    numelements = datamat.shape[1]
    data_reshaped = datamat.reshape((dumpsteps, natoms, numelements)) 
    return data_reshaped

def rwd(data_reshaped, dumpsteps):
    Weights = np.zeros((dumpsteps, natoms, Nn))
    Qreward = np.zeros((dumpsteps, natoms))
    Positions = np.zeros((dumpsteps, natoms, 2))
    Speeds = np.zeros((dumpsteps, natoms, 2))
    for t in range(dumpsteps):
        for i in range(natoms):
            c = 1
            Qreward[t,i] = data_reshaped[t,i,c]
            c +=1 
            for k in range(Nn):
                Weights[t, i, k] = data_reshaped[t, i, c]  
                c +=1 
            Positions[t,i,0] = data_reshaped[t,i,c]
            c +=1
            Positions[t,i,1] = data_reshaped[t,i,c]
            c+=1
            Speeds[t,i,0] = data_reshaped[t,i,c]
            c +=1
            Speeds[t,i,1] = data_reshaped[t,i,c]


    return Qreward, Weights, Positions, Speeds


input_file = "evo.in"

dp = [1e-5]  # Valeurs de Dp à tester
nump = len(dp)
alphaq = [0.03] 
numq = len(alphaq)
Seed = np.zeros((nump, numq, numrep))

for p in range(nump):
    Dp = dp[p]
    for i in range(numq):
        alpha_q = alphaq[i]
        for rep in range(numrep):
            seed = random.randint(10000, 99999)
            Seed[p,i,rep] = seed
            log_file = f"log_seed_{seed}_{Dp}_{alpha_q}.lmp"

            print(f"Lancement de la simulation avec seed = {seed}, Dp = {Dp}, alpha_q = {alpha_q}, {rep+1}/{numrep}")
    
            # Exécute LAMMPS avec la seed variable
            subprocess.run(["mpirun", "--oversubscribe","-np",str(ncore),
                "/home/johann1/Desktop/lammps/src/lmp_mpi",
                "-in", input_file,
                "-var", "Seed", str(seed),
                "-var", "Dp", str(Dp),
                "-var", "alphaq", str(alpha_q)
            ], stdout=open(log_file, "w"))

            
            
            # enlever les lignes parasites
            monfichier, dumpsteps = filtrer_fichier(f'dumpnn_{Dp}_{alpha_q}_{seed}.txt', natoms, parasites, f'filtered_{Dp}_{alpha_q}_{seed}.txt')
            # fichier brut vers disque pour visualiser ovito
            shutil.move(f'dumpnn_{Dp}_{alpha_q}_{seed}.txt', f'/media/johann1/DATA/Johann/END_to_END_cfg_6/{natoms}/dumpnn_{Dp}_{alpha_q}_{seed}.txt')

            
            simu = fichier_to_matrice(natoms, dumpsteps, f'filtered_{Dp}_{alpha_q}_{seed}.txt')
            score, Poids, Pos, Speeds = rwd(simu, dumpsteps)
            os.remove(f'filtered_{Dp}_{alpha_q}_{seed}.txt')


            np.save(f'Weights_{Dp}_{alpha_q}_{seed}.npy', Poids)
            shutil.move(f'Weights_{Dp}_{alpha_q}_{seed}.npy', f'/media/johann1/DATA/Johann/END_to_END_cfg_6/{natoms}/Weights_{Dp}_{alpha_q}_{seed}.npy')
            np.save(f'Rewards_{Dp}_{alpha_q}_{seed}.npy', score)
            shutil.move(f'Rewards_{Dp}_{alpha_q}_{seed}.npy', f'/media/johann1/DATA/Johann/END_to_END_cfg_6/{natoms}/Rewards_{Dp}_{alpha_q}_{seed}.npy')
            np.save(f'Positions_{Dp}_{alpha_q}_{seed}.npy', Pos)
            shutil.move(f'Positions_{Dp}_{alpha_q}_{seed}.npy', f'/media/johann1/DATA/Johann/END_to_END_cfg_6/{natoms}/Positions_{Dp}_{alpha_q}_{seed}.npy')
            np.save(f'Speeds_{Dp}_{alpha_q}_{seed}.npy', Speeds)
            shutil.move(f'Speeds_{Dp}_{alpha_q}_{seed}.npy', f'/media/johann1/DATA/Johann/END_to_END_cfg_6/{natoms}/Speeds_{Dp}_{alpha_q}_{seed}.npy')
            
            
            
np.save('Seed.npy',Seed)