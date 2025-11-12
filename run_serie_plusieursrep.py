import subprocess
import random
import numpy as np
import matplotlib.pyplot as plt
import os

# === PARAMÈTRES ===
natoms = 64    # Nombre de particules
parasites = 9        # Lignes parasites par timestep
dumpfreq = 1000      # Fréquence des dumps
time = 14400000
Nn = 8
numrep = 10

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

def filter_x(input_filepath, natoms, parasites, output_filepath):
    lignesn = []
    with open(input_filepath, 'r') as f:
        lignes = f.readlines()

    for i in range(len(lignes)):
        if i % (natoms + parasites) > parasites - 1:
            lignesn.append(lignes[i])
    with open(output_filepath, "w") as fa:
        fa.writelines(lignesn)

    return output_filepath


def fichier_to_matrice(natoms, dumpsteps, input_filepath):
    datamat = np.loadtxt(input_filepath)
    numelements = datamat.shape[1]
    data_reshaped = datamat.reshape((dumpsteps, natoms, numelements)) 
    return data_reshaped


def rwd(data_reshaped, dumpsteps):
    Weights = np.zeros((dumpsteps, natoms, Nn))
    Qreward = np.zeros((dumpsteps, natoms))

    for t in range(dumpsteps):
        for i in range(natoms):
            Qreward[t,i] = data_reshaped[t,i,1]
            for k in range(Nn):
                Weights[t, i, k] = data_reshaped[t, i, k + 2]  # Poids NN
    return Qreward, Weights

def pos(data_reshaped, dumpsteps):
    Positions = np.zeros((dumpsteps, natoms, 2))

    for t in range(dumpsteps):
        for i in range(natoms):
            Positions[t,i,0] = data_reshaped[t,i,1]
            Positions[t,i,1] = data_reshaped[t,i,2]
    return Positions

# Chemin vers ton fichier input
input_file = "evo.in"

dp = [1e-5]  # Valeurs de Dp à tester
nump = len(dp)
alphaq = [0.0083] 
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
            subprocess.run([
                "/home/atef/Documents/Johann_Cours_ENS/stage_S2/lammps/src/lmp_serial",
                "-in", input_file,
                "-var", "Seed", str(seed),
                "-var", "Dp", str(Dp),
                "-var", "alphaq", str(alpha_q)
            ], stdout=open(log_file, "w"))

            monfichier, dumpsteps = filtrer_fichier(f'dumpnn_{Dp}_{alpha_q}_{seed}.txt', natoms, parasites, f'filtered_{Dp}_{alpha_q}_{seed}.txt')
            ficher_pos = filter_x(f'Evolopti_{natoms}_4_{alpha_q}_{Dp}_{seed}.cfg', natoms, parasites, f'pos_{Dp}_{alpha_q}_{seed}.txt')
            simu_p = fichier_to_matrice(natoms, dumpsteps, f'filtered_{Dp}_{alpha_q}_{seed}.txt')
            simu_pos = fichier_to_matrice(natoms,dumpsteps,f'pos_{Dp}_{alpha_q}_{seed}.txt')
            score, Poids = rwd(simu_p, dumpsteps)
            Pos = pos(simu_pos, dumpsteps)

            np.save(f'Weights_{Dp}_{alpha_q}_{seed}.npy', Poids)
            np.save(f'Rewards_{Dp}_{alpha_q}_{seed}.npy', score)
            np.save(f'Positions_{Dp}_{alpha_q}_{seed}.npy', Pos)
            os.remove(f'filtered_{Dp}_{alpha_q}_{seed}.txt')
            os.remove(f'pos_{Dp}_{alpha_q}_{seed}.txt')
            os.remove(f'dumpnn_{Dp}_{alpha_q}_{seed}.txt')
np.save('Seed.npy',Seed)